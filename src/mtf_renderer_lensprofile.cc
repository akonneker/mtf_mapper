/*
Copyright 2011 Frans van den Bergh. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY Frans van den Bergh ''AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Frans van den Bergh OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of the Council for Scientific and Industrial Research (CSIR).
*/

#include "include/mtf_renderer_lensprofile.h"
#include "include/sampling_rate.h"
#include <opencv2/imgcodecs/imgcodecs.hpp>

void Mtf_renderer_lensprofile::render(const vector<Block>& blocks) {
    Point2d centr(img.cols/2, img.rows/2);
    
    vector<double> resolution;
    for (size_t i=0; i < min(size_t(3), in_resolution.size()); i++) {
        resolution.push_back(in_resolution[i] / pixel_size);
        logger.debug("resolution = %lf (in= %lf, ps=%lf)\n", resolution.back(), in_resolution[i], pixel_size);
    }
    
    vector< vector<Ordered_point> > sagittal(resolution.size());
    vector< vector<Ordered_point> > meridional(resolution.size());
    for (size_t i=0; i < blocks.size(); i++) {
    
        const double angle_thresh = 15.0;
        
        for (size_t k=0; k < 4; k++) {
            if (blocks[i].get_mtf50_value(k) >= 1.0 || blocks[i].get_quality(k) <= 0.3) continue;
        
            Point2d ec = blocks[i].get_edge_centroid(k);
            
            Point2d udir = ec - centr;
            double radial_len = norm(udir);
            Point2d dir = udir * (1.0/radial_len);

            Point2d norm = blocks[i].get_normal(k);
            double delta = dir.x*norm.x + dir.y*norm.y;
            
            const vector<double>& sfr = blocks[i].get_sfr(k);
            for (size_t j=0; j < resolution.size(); j++) {
                double res = resolution[j] * NYQUIST_FREQ*2;
                int lidx = min((int)floor(res), NYQUIST_FREQ*2-2);
                double frac = res - lidx;
            
                double contrast = (1 - frac)*sfr[lidx] + frac*sfr[lidx+1]; 
                
                //contrast /= sin(resolution[j]*M_PI)/(resolution[j]*M_PI); // TODO: optionally remove photosite aperture / AA MTF
                
                if (acos(fabs(delta))/M_PI*180.0 < angle_thresh) { // edge perp to tangent
                    meridional[j].push_back(Ordered_point(radial_len, contrast));
                } 
                if (acos(fabs(delta))/M_PI*180 > (90 - angle_thresh)) { // edge perp to radial : TODO: check math
                    sagittal[j].push_back(Ordered_point(radial_len, contrast));
                }
            }
        }
    }    
    
    FILE* fout = fopen((wdir + prname).c_str(), "wt");
    
    if (sagittal[0].size() == 0 || meridional.size() == 0) {
        logger.error("%s\n", "Fatal error: lens profile requested, but insufficient edges detected to generate profile.\nSkipping.");
        return;
    }
    
    logger.debug("got %d sagittal / %d meridional samples\n", (int)sagittal[0].size(), (int)meridional[0].size());
    
    vector< vector<Ordered_point> > s_fitted(resolution.size());
    vector< vector<Ordered_point> > s_spread(resolution.size());
    
    vector< vector<Ordered_point> > m_fitted(resolution.size());
    vector< vector<Ordered_point> > m_spread(resolution.size());
    
    
    fprintf(fout, "# ");
    for (size_t j=0; j < resolution.size(); j++) {
    
        sort(sagittal[j].begin(), sagittal[j].end());
        sort(meridional[j].begin(), meridional[j].end());
        
        lsfit(sagittal[j], s_fitted[j], s_spread[j]);
        lsfit(meridional[j], m_fitted[j], m_spread[j]);
        
        fprintf(fout, "distance  contrast(%.1flp/mm)  ", resolution[j] * pixel_size);
    }
    fprintf(fout, "\n");
    
    double scale = 1.0/pixel_size;
    
    /*
    FILE* fraw = fopen((wdir + "lp_raw.txt").c_str(), "wt");
    const int ridx = 1;
    for (size_t i=0; i < sagittal[ridx].size(); i++) {
        fprintf(fraw, "%lf %lf\n", scale*sagittal[ridx][i].first, sagittal[ridx][i].second);
    }
    fprintf(fraw, "\n\n");
    for (size_t i=0; i < meridional[0].size(); i++) {
        fprintf(fraw, "%lf %lf\n", scale*meridional[ridx][i].first, meridional[ridx][i].second);
    }
    fclose(fraw);
    */
    
    double lower_limit = 5;
    double upper_limit = -5;
    
    fprintf(fout, "#sagittal curve\n");
    for (size_t i=0; i < s_fitted[0].size(); i++) {
        for (size_t j=0; j < resolution.size(); j++) {
            fprintf(fout, "%lf %lf ", scale*s_fitted[j][i].first, s_fitted[j][i].second);
            lower_limit = min(lower_limit, s_fitted[j][i].second);
            upper_limit = max(upper_limit, s_fitted[j][i].second);
        }
        fprintf(fout, "\n");
    }
    
    fprintf(fout, "\n\n#meridional curve\n");
    for (size_t i=0; i < m_fitted[0].size(); i++) {
        for (size_t j=0; j < resolution.size(); j++) {
            fprintf(fout, "%lf %lf ", scale*m_fitted[j][i].first, m_fitted[j][i].second);
            lower_limit = min(lower_limit, m_fitted[j][i].second);
            upper_limit = max(upper_limit, m_fitted[j][i].second);
        }
        fprintf(fout, "\n");
    }
    
    
    fprintf(fout, "\n\n#sagittal bounds\n");
    for (size_t i=0; i < s_spread[0].size(); i++) {
        for (size_t j=0; j < resolution.size(); j++) {
            fprintf(fout, "%lf %lf ", scale*s_spread[j][i].first, s_spread[j][i].second);
            lower_limit = min(lower_limit, s_spread[j][i].second);
            upper_limit = max(upper_limit, s_spread[j][i].second);
        }
        fprintf(fout, "\n");
    }
    
    fprintf(fout, "\n\n#meridional bounds\n");
    for (size_t i=0; i < m_spread[0].size(); i++) {
        for (size_t j=0; j < resolution.size(); j++) {
            fprintf(fout, "%lf %lf ", scale*m_spread[j][i].first, m_spread[j][i].second);
            lower_limit = min(lower_limit, m_spread[j][i].second);
            upper_limit = max(upper_limit, m_spread[j][i].second);
        }
        fprintf(fout, "\n");
    }
    fclose(fout);
    
    vector<string> linecolor;
    linecolor.push_back("#e00000");
    linecolor.push_back("web-blue");
    linecolor.push_back("web-green");
    
    vector<string> shadecolor;
    shadecolor.push_back("#ffe0e0");
    shadecolor.push_back("#f0d0d0");
    shadecolor.push_back("#d0d0f7");
    shadecolor.push_back("#d0d0f0");
    shadecolor.push_back("#d0f7d0");
    shadecolor.push_back("#d0f0d0");
    
    string resmode = lpmm_mode ? "lp/mm" : "c/p";
    
    FILE* gpf = fopen( (wdir + string("lensprofile.gnuplot")).c_str(), "wt");
    fprintf(gpf, "set xlab \"distance (%s)\"\n", lpmm_mode ? "mm" : "pixels");
    fprintf(gpf, "set ylab \"contrast\"\n");
    fprintf(gpf, "set key left bottom\n");
    fprintf(gpf, "set ytics 0.1\n");
    fprintf(gpf, "set style line 11 lc rgb \"#f0f0f0\" lt 1 lw 1\n");
    fprintf(gpf, "set grid xtics ytics ls 11\n");
    double ar = 768.0/1024.0;
    int fontsize = lrint(12.0*gnuplot_width/1024.0);
    int title_fontsize = lrint(14.0*gnuplot_width/1024.0);
    int linewidth = lrint(2*gnuplot_width/1024.0);
    #ifdef _WIN32
    fprintf(gpf, "set term pngcairo dashed transparent enhanced size %d, %d font 'Verdana,%d'\n", gnuplot_width, int(gnuplot_width*ar), fontsize);
    #else
    fprintf(gpf, "set term pngcairo dashed transparent enhanced size %d, %d font 'Arial,%d'\n", gnuplot_width, int(gnuplot_width*ar), fontsize);
    #endif
    
    fprintf(gpf, "if (GPVAL_VERSION >= 5.0) set linetype 12 dashtype 2;\n"); // force dashed lines in meridional plot
    fprintf(gpf, "set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb \"white\" behind\n");
    fprintf(gpf, "set output \"%slensprofile.png\"\n", wdir.c_str());
    if (img_filename.length() > 0) {
        fprintf(gpf, "set title \"%s\" font \",%d\"\n", img_filename.c_str(), title_fontsize);
    }
    
    lower_limit = min(-0.05, lower_limit);
    upper_limit = max(1.0, upper_limit);
    if (fixed_size) {
        double max_radius = sqrt(img.cols*img.cols/4 + img.rows*img.rows/4);
        max_radius /= pixel_size;
        fprintf(gpf, "plot [0:%.3lf][%.3lf:%.3lf] ", max_radius, lower_limit, upper_limit);
    } else {
        fprintf(gpf, "plot [][%.3lf:%.3lf] ", lower_limit, upper_limit);
    }
    for (size_t j=0; j < resolution.size(); j++) {
        double res = lpmm_mode ? resolution[j]*pixel_size : resolution[j] ;
        fprintf(gpf,   
            "\"%s\" index 2 u 1:%d w filledcurve fs transparent solid 0.5 lc rgb \"%s\" lt 16 notitle,"
            "\"%s\" index 3 u 1:%d w filledcurve fs transparent solid 0.5 lc rgb \"%s\" lt 16 notitle,"
            "\"%s\" index 0 u 1:%d w l lc rgb \"%s\" lw %d lt 16 t \"S %.2lf %s\","
            "\"%s\" index 1 u 1:%d w l lc rgb \"%s\" lt 12 lw %d t \"M %.2lf %s\"",
            (wdir+prname).c_str(), int(2*j)+2, shadecolor[2*j].c_str(),
            (wdir+prname).c_str(), int(2*j)+2, shadecolor[2*j+1].c_str(),
            (wdir+prname).c_str(), int(2*j)+2, linecolor[j].c_str(), linewidth, res, resmode.c_str(),
            (wdir+prname).c_str(), int(2*j)+2, linecolor[j].c_str(), linewidth, res, resmode.c_str()
        );
        fprintf(gpf, "%c", (j < resolution.size() - 1) ? ',' : '\n');
    }    
    fclose(gpf);
    
    char* buffer = new char[1024];
    #ifdef _WIN32
    sprintf(buffer, "\"\"%s\" \"%slensprofile.gnuplot\"\"", gnuplot_binary.c_str(), wdir.c_str());
    #else
    sprintf(buffer, "\"%s\" \"%slensprofile.gnuplot\"", gnuplot_binary.c_str(), wdir.c_str());
    #endif
    int rval = system(buffer);
    if (rval != 0) {
        logger.error("Failed to execute gnuplot (error code %d)\n", rval);
        logger.info("You can try to execute [%s] to render the plots manually\n", buffer);
        gnuplot_failure = true;
    } else {
        logger.debug("%s\n", "Gnuplot plot completed successfully. Look for lensprofile.png");
        // Hack: trim the png file to remove the black border that gnuplot 4.4 leaves
        cv::Mat src_img;
        src_img = cv::imread(wdir.c_str() + string("/lensprofile.png"));
        if (src_img.data) {
            src_img.adjustROI(-1, -1, -1, -1);
            cv::imwrite(wdir.c_str() + string("/lensprofile.png"), src_img);
        }
    }
    
    delete [] buffer;
}

void Mtf_renderer_lensprofile::lsfit(const vector<Ordered_point>& in_data, 
    vector<Ordered_point>& recon, 
    vector<Ordered_point>& spread, int recon_samples) {
    
    recon  = vector<Ordered_point> (recon_samples);
    spread  = vector<Ordered_point> (recon_samples, Ordered_point(0,-1.0));
    
    vector<Ordered_point> data(in_data);
    
    double x_span = data.back().first - data.front().first;
    
    for (int i=0; i < 5; i++) {
        data.push_back(Ordered_point(in_data[i].first - 0.01*x_span, in_data[i].second));
        data.push_back(Ordered_point(in_data[i].first - 0.015*x_span, in_data[i+1].second));
        data.push_back(Ordered_point(in_data[in_data.size()-1-i].first + 0.01*x_span, in_data[in_data.size()-1-i].second));
        data.push_back(Ordered_point(in_data[in_data.size()-1-i].first + 0.015*x_span, in_data[in_data.size()-1-i-1].second));
    }
    sort(data.begin(), data.end());
    
    int tdim = 2;
    double scale = 5.0;
    double h = x_span/scale; 
    double bin_width = x_span/double(recon.size()-1);
    
    size_t upper_idx = 0;
    size_t lower_idx = 0;
    double prev_scale = scale;
    for (size_t q=0; q < recon.size(); q++) {
        prev_scale = scale;
        if (q < recon.size()*0.2 || q >= recon.size()*0.8 || data.size() < 80) {
            tdim = 3;
            scale = 5.0;
        } else {
            tdim = 4;
            scale = 5.0;
        }
        h = x_span/scale;
        if (scale != prev_scale) {
            lower_idx = 0;
            upper_idx = 0;
        }
        
        double xmid = (double(q)-0.5)*x_span/double(recon.size()-2) + in_data.front().first;
        
        for (; lower_idx < data.size() && data[lower_idx].first < (xmid-h); lower_idx++);
        for (; upper_idx < data.size() && data[upper_idx].first < (xmid+h); upper_idx++);
        
        if (sparse_chart) {
            int cur_samples = upper_idx - lower_idx;
            while (data.size() > 20 && cur_samples < 20) {
                if (lower_idx > 0) {
                    lower_idx--;
                    cur_samples++;
                }
                if (upper_idx < data.size()-1) {
                    upper_idx++;
                    cur_samples++;
                }
            }
        }
        
        vector< vector<double> > cov(tdim, vector<double>(tdim, 0.0));
        vector<double> b(tdim, 0.0);
        vector<double> a(tdim);
        
        for (size_t r=lower_idx; r < upper_idx; r++) { 
            double sh = 1.0/scale;
            double x = (data[r].first - xmid)/x_span; 
            if (fabs(x) > 1) x /= fabs(x);
            double w = fabs(x)*fabs(x)*fabs(x)/(sh*sh*sh);
            w = fabs(x) < sh ? (1 - w)*(1 - w)*(1 - w) : 0;
            w += sparse_chart ? 0.25 : 0; // if we have very few samples, we must try to use them to bridge the gaps
            
            a[0] = w;
            double px = x;
            for (int jj=1; jj < tdim; jj++) {
                a[jj] = w*px;
                px *= x;
            }
            for (int col=0; col < tdim; col++) {
                for (int icol=0; icol < tdim; icol++) { // build covariance of design matrix : A'*A
                    cov[col][icol] += a[col]*a[icol];
                }
                b[col] += a[col]*data[r].second*w; // build rhs of system : A'*b
            }
        }
        
        // hardcode cholesky decomposition
        bool singular = false;
        vector<double> ldiag(tdim, 0.0);
        for (int i=0; i < tdim && !singular; i++) {
            for (int j=i; j < tdim && !singular; j++) {
                double sum = cov[i][j];
                for (int k=i-1; k >= 0; k--) {
                    sum -= cov[i][k]*cov[j][k];
                }
                if (i == j) {
                    if (sum <= 0.0) {
                        singular = true;
                    } else {
                        ldiag[i] = sqrt(sum);
                    }
                } else {
                    cov[j][i] = sum/ldiag[i];
                }
            }
        }
        
        if (!singular) {
            // hardcode backsubstitution
            vector<double> x(tdim, 0.0);
            for (int i=0; i < tdim; i++) {
                double sum = b[i];
                for (int k=i-1; k >= 0; k--) {
                    sum -= cov[i][k]*x[k];
                }
                x[i] = sum/ldiag[i];
            }
            for (int i=tdim-1; i >= 0; i--) {
                double sum = x[i];
                for (int k=i+1; k < tdim; k++) {
                    sum -= cov[k][i]*x[k];
                }
                x[i] = sum/ldiag[i];
            }
        
            
            recon[q].first = xmid;
            recon[q].second = x[0]; // x is centred on xmid, so parabola constant is midpoint extimate
            
            vector<double> residuals;
            for (size_t r=lower_idx; r < upper_idx; r++) {
                // only compute residuals on really close points?
                double cent_x = data[r].first - xmid;
                
                if (fabs(cent_x) < 2*bin_width) {
                    double delta = data[r].second - polyeval(cent_x/x_span, x);
                    residuals.push_back(fabs(delta));
                }
            }
            if (residuals.size() < 10) {
                for (size_t r=lower_idx; r < upper_idx; r++) {
                    double cent_x = data[r].first - xmid;
                    
                    if (fabs(cent_x) < 4*bin_width) {
                        double delta = data[r].second - polyeval(cent_x/x_span, x);
                        residuals.push_back(fabs(delta));
                    }
                }
            }
            if (residuals.size() > 2) {
                sort(residuals.begin(), residuals.end());
                double p90 = residuals[lrint((residuals.size() - 1)*0.9)];
                spread[q].first = spread[q].second = p90;
            }
        }
    }
    
    size_t nz_idx = 0;
    while (nz_idx < spread.size() && spread[nz_idx].second == -1) nz_idx++;
    if (nz_idx == spread.size()) { // no nonzero points at all!
        nz_idx = 0;
        spread.front() = spread.back() = Ordered_point(1,1);
    } else {
        for (size_t q=0; q < nz_idx; q++) { // pad out leading zeros with first nonzero
            spread[q] = spread[nz_idx];
        }
    }
    size_t lnz_idx = spread.size() - 1;
    while (lnz_idx > 0 && spread[lnz_idx].second == -1) lnz_idx--;
    for (size_t q=lnz_idx; q < spread.size(); q++) { // pad out trailing zeros with last nonzero
        spread[q] = spread[lnz_idx];
    }
    // now we have valid endpoints, so we can interpolate
    size_t prev_nz = nz_idx;
    for (size_t q=nz_idx+1; q < lnz_idx; q++) {
        if (spread[q].second == -1) {
            // find next nonzero index
            size_t next_nz = q+1;
            while (next_nz < lnz_idx && spread[next_nz].second == -1) next_nz++;
            double slope = (spread[next_nz].second - spread[prev_nz].second) / (recon[next_nz].first - recon[prev_nz].first);
            double offset = spread[prev_nz].second - slope*recon[prev_nz].first;
            while (q < next_nz) {
                spread[q].first = spread[q].second = recon[q].first * slope + offset;
                q++;
            }
        } 
        prev_nz = q;
    }
    
    // now hit recon with SG smoothing
    const int sgh = 7;
    vector<double> recon_smoothed(recon.size(), 0);
    const double* sgw = 0;
    for (int q=0; q < (int)recon.size(); q++) {
        if (q < sgh || q > (int(spread.size()) - 1 - sgh)) {
            recon_smoothed[q] = recon[q].second;
        } else {
            //const int stride = 3;
            int filter_order = 5; //min(5, (q-5)/stride);
            sgw = savitsky_golay[filter_order];
            for (int x=-sgh; x <= sgh; x++) { 
                // since magnitude has extra samples at the end, we can safely go past the end
                recon_smoothed[q] += recon[q+x].second * sgw[x+sgh];
            }
        }
    }
    for (int q=0; q < (int)recon.size(); q++) {
        recon[q].second = recon_smoothed[q];
    }
    
    // apply some exponential smoothing to spread
    vector<double> smoothed(spread.size(), 0);
    double fsmooth = spread.front().first;
    double bsmooth = spread.back().first;
    double alpha = 0.35;
    vector<double> fs(spread.size(), 0);
    vector<double> bs(spread.size(), 0);
    for (size_t q=0; q < spread.size(); q++) {
        fsmooth = (1 - alpha)*fsmooth + alpha * spread[q].first;
        bsmooth = (1 - alpha)*bsmooth + alpha * spread[spread.size()-1 - q].first;
        fs[q] = fsmooth;
        bs[spread.size()-1 - q] = bsmooth;
    }
    for (size_t q=0; q < spread.size(); q++) {
        smoothed[q] = 0.5*fs[q] + 0.5*bs[q];
        if (q < 10) {
            smoothed[q] = bs[q];
        }
        if (q > (spread.size()-1-10)) {
            smoothed[q] = fs[q];
        }
    }
    
    spread = vector<Ordered_point>(spread.size()*2);
    // every element of spread is now interpolated, so add spread to recon
    for (size_t q=0; q < smoothed.size(); q++) {
        spread[q].first = recon[q].first;
        spread[q].second = recon[q].second + smoothed[q];
        
        spread[spread.size()-1-q].first = recon[q].first;
        spread[spread.size()-1-q].second = recon[q].second - smoothed[q];
    }
}
    

