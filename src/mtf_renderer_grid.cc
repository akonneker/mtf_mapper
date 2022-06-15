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
#include "include/mtf_renderer_grid.h"
#include <opencv2/imgproc/imgproc.hpp>

class Grid_functor_mtf : public Grid_functor {
  public:
    Grid_functor_mtf(bool sparse_mode = false) : sparse_mode(sparse_mode) {}

    virtual double value(const Block& block, size_t edge) const {
        return  block.get_quality(edge) > very_poor_quality ? block.get_mtf50_value(edge) : nodata();
    }
    
    virtual bool in_range(double val) const {
        return val > 0 && val < 1.0; 
    }
    
    virtual double clamp(double value, double upper, double /*lower*/) const {
        return std::max(0.0, std::min(value, upper));
    }
    
    virtual double nodata(void) const {
        return 2.0;
    }

    virtual double smoothing_factor(void) const {
        return sparse_mode ? 1e-1 : 1e-3;
    }

    virtual int pruning_threshold(void) const {
        return sparse_mode ? 2 : 1;
    }

  private:
    bool sparse_mode = false;
};

Mtf_renderer_grid::Mtf_renderer_grid(
    const std::string& img_filename,
    const std::string& wdir, const std::string& fname, 
    const std::string& gnuplot_binary, 
    const cv::Mat& img, int gnuplot_width,
    bool lpmm_mode, double pixel_size,
    double in_zscale, 
    double surface_max,
    int mtf_contrast)
    :  Mtf_renderer(img_filename),
        wdir(wdir), fname(fname), 
        gnuplot_binary(gnuplot_binary),
        img(img), lpmm_mode(lpmm_mode), pixel_size(pixel_size),
        gnuplot_failure(false), gnuplot_warning(true),
        m_lower(0), m_upper(0), gnuplot_width(gnuplot_width),
        surface_max(surface_max), mtf_contrast(mtf_contrast) {

    const int coarse_grid_size = 40;
    const int fine_grid_size = 200;
    if (img.rows > img.cols) {
        grid_y_coarse = coarse_grid_size;
        grid_x_coarse = coarse_grid_size * img.cols / img.rows;
        grid_y_fine = fine_grid_size;
        grid_x_fine = fine_grid_size * img.cols / img.rows;
    } else {
        grid_x_coarse = coarse_grid_size;
        grid_y_coarse = coarse_grid_size * img.rows / img.cols;
        grid_x_fine = fine_grid_size;
        grid_y_fine = fine_grid_size * img.rows / img.cols;
    }
    
    zscale = std::max(0.0, std::min(in_zscale, 1.0));
}
    
void Mtf_renderer_grid::render(const vector<Block>& blocks) {

    // first check if we have enough good data to generate a plot
    vector<double> allvals;
    for (size_t i=0; i < blocks.size(); i++) {
        for (size_t k=0; k < 4; k++) {
            double val = blocks[i].get_mtf50_value(k);
            if (val > 0 && val < 1.0 && blocks[i].get_quality(k) > very_poor_quality) {
                allvals.push_back(val);
            }
        }
    }

    if (allvals.size() < 20) {
        logger.error("%s\n", "Too few valid edges found. No surface can be generated.");
        return;
    }

    sort(allvals.begin(), allvals.end());
    m_upper = allvals[97*allvals.size()/100];
    m_lower = allvals[5*allvals.size()/100];
    
    cv::Mat grid_mer_coarse(grid_y_coarse, grid_x_coarse, CV_32FC1, 0.0);
    cv::Mat grid_mer_fine(grid_y_fine, grid_x_fine, CV_32FC1, 0.0);
    cv::Mat grid_sag_coarse(grid_y_coarse, grid_x_coarse, CV_32FC1, 0.0);
    cv::Mat grid_sag_fine(grid_y_fine, grid_x_fine, CV_32FC1, 0.0);
    
    Grid_functor_mtf mtf_ftor(sparse_chart);
    cv::Size img_dims(img.cols, img.rows);
    interpolate_grid(mtf_ftor, MERIDIONAL, grid_mer_coarse, grid_mer_fine, img_dims, blocks, m_upper, mtf_ftor.smoothing_factor(), mtf_ftor.pruning_threshold());
    interpolate_grid(mtf_ftor, SAGITTAL, grid_sag_coarse, grid_sag_fine, img_dims, blocks, m_upper, mtf_ftor.smoothing_factor(), mtf_ftor.pruning_threshold());
    
    double zmax = 0;
    double zmin = 1e30;
    FILE* file = fopen((wdir+fname).c_str(), "wt");
    fprintf(file, "#coarse meridional grid\n");
    for (int y=0; y < grid_mer_coarse.rows; y++) {
        for (int x=0; x < grid_mer_coarse.cols; x++) {
            fprintf(file, "%lf %lf %.5lf\n", 
                x*img.cols/grid_mer_coarse.cols/pixel_size, 
                y*img.rows/grid_mer_coarse.rows/pixel_size, 
                grid_mer_coarse.at<float>(y,x)*pixel_size
            );
            zmax = std::max(zmax, (double)grid_mer_coarse.at<float>(y, x) * pixel_size);
            zmin = std::min(zmin, (double)grid_mer_coarse.at<float>(y, x) * pixel_size);
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n\n");
    fprintf(file, "#coarse sagittal grid\n");
    for (int y=0; y < grid_sag_coarse.rows; y++) {
        for (int x=0; x < grid_sag_coarse.cols; x++) {
            fprintf(file, "%lf %lf %.5lf\n", 
                x*img.cols/grid_sag_coarse.cols/pixel_size, 
                y*img.rows/grid_sag_coarse.rows/pixel_size, 
                grid_sag_coarse.at<float>(y,x)*pixel_size
            );
            zmax = std::max(zmax, (double)grid_sag_coarse.at<float>(y, x) * pixel_size);
            zmin = std::min(zmin, (double)grid_sag_coarse.at<float>(y, x) * pixel_size);
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n\n");
    
#if 0
    // we no longer use the fine grid, but rather let gnuplot interpolate by itself
    fprintf(file, "#fine meridional grid\n");
    for (int y=0; y < grid_mer_fine.rows; y++) {
        for (int x=0; x < grid_mer_fine.cols; x++) {
            fprintf(file, "%lf %lf %.5lf\n", 
                x*img.cols/grid_mer_fine.cols/pixel_size, 
                y*img.rows/grid_mer_fine.rows/pixel_size, 
                grid_mer_fine.at<float>(y,x)*pixel_size
            );
            zmax = std::max(zmax, (double)grid_mer_fine.at<float>(y,x)*pixel_size);
            zmin = std::min(zmin, (double)grid_mer_fine.at<float>(y,x)*pixel_size);
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n\n");
    
    fprintf(file, "#fine sagittal grid\n");
    for (int y=0; y < grid_sag_fine.rows; y++) {
        for (int x=0; x < grid_sag_fine.cols; x++) {
            fprintf(file, "%lf %lf %.5lf\n", 
                x*img.cols/grid_sag_fine.cols/pixel_size, 
                y*img.rows/grid_sag_fine.rows/pixel_size, 
                grid_sag_fine.at<float>(y,x)*pixel_size
            );
            zmax = std::max(zmax, (double)grid_sag_fine.at<float>(y,x)*pixel_size);
            zmin = std::min(zmin, (double)grid_sag_fine.at<float>(y,x)*pixel_size);
        }
        fprintf(file, "\n");
    }
#endif
    
    fclose(file);
    
    if (surface_max > 0 && surface_max > zmax) { // valid user override on maximum value specified
        zmax = surface_max;
    }
    
    double span = zmax - zmin;
    zmax = zmax + 0.05*span;
    zmin = zmin - 0.05*span;
    zmin = std::max(zmin, 0.0);
    
    // adapt zmin according to zscale setting
    // zscale = 0 -> use zmin=0
    // zscale = 1 -> use measured zmin
    zmin *= zscale;
    
    const int width_in_pixels = int(gnuplot_width*600.0/1024);
    const int height_in_pixels_3d = int(gnuplot_width*1200.0/1024);
    int fontsize = std::max(long(10), lrint(10.0*gnuplot_width/1024.0));
    int title_fontsize = fontsize + 2;
    int fontsize3d = std::max(long(9), lrint(9.0*gnuplot_width/1024.0));
    double linewidth = std::max(double(0.5), double(0.5)*gnuplot_width/1024.0);
    
    FILE* gpf = fopen((wdir + std::string("grid.gnuplot")).c_str(), "wt");
    
    fprintf(gpf, "%s\n", diverging_palette.c_str());
    fprintf(gpf, "set cbrange [%lf:%lf]\n", zmin, zmax);
    fprintf(gpf, "set xlab \"column (%s)\"\n", lpmm_mode ? "mm" : "pixels");
    fprintf(gpf, "set ylab \"row (%s)\"\n",  lpmm_mode ? "mm" : "pixels");
    fprintf(gpf, "set xtics out nomirror\n");
    fprintf(gpf, "set xtics scale 0.3\n");
    fprintf(gpf, "set ytics out nomirror\n");
    fprintf(gpf, "set ytics scale 0.3\n");
    if (!lpmm_mode) {
        fprintf(gpf, "set ytics rotate by 45\n");
    }
    fprintf(gpf, "set pm3d map impl interpolate 3,3\n");
    fprintf(gpf, "set hidden3d\n");
    fprintf(gpf, "set cntrlabel onecolor\n");
    fprintf(gpf, "set contour surface\n");
    fprintf(gpf, "set cntrparam order 8\n");
    fprintf(gpf, "set cntrparam bspline\n");
    fprintf(gpf, "set autoscale xfix\n");
    fprintf(gpf, "set autoscale yfix\n");
    fprintf(gpf, "set term pngcairo dashed transparent enhanced size %d, %d font '%s,%d'  background rgb \"white\"\n",
        width_in_pixels, 
        (int)lrint(width_in_pixels*2.3*grid_mer_fine.rows/double(grid_mer_fine.cols)), 
        #ifdef _WIN32
        "Verdana",
        #else
        "Arial",
        #endif
        fontsize
    );
    
    fprintf(gpf, "set output \"%sgrid_image.png\"\n", wdir.c_str());
    if (img_filename.length() > 0) {
        fprintf(gpf, "set multiplot title \"%s\" font \",%d\"\n", img_filename.c_str(), title_fontsize);
        fprintf(gpf, "set tmargin 4\n");
    } else {
        fprintf(gpf, "set multiplot\n");
    }
    fprintf(gpf, "set size 1,0.5\n");
    fprintf(gpf, "set origin 0.0,0.5\n");
    fprintf(gpf, "set title \"Meridional MTF%2d (%s)\"\n", mtf_contrast, lpmm_mode ? "lp/mm" : "c/p");
    fprintf(gpf, "set yrange [*:*] reverse\n");
    fprintf(gpf, "splot \"%s\" i 0 notitle w l lc rgb \"#77303030\"\n", (wdir+fname).c_str());
    fprintf(gpf, "set origin 0.0,0.0\n");
    fprintf(gpf, "set title \"Sagittal MTF%2d (%s)\"\n", mtf_contrast, lpmm_mode ? "lp/mm" : "c/p");
    fprintf(gpf, "splot \"%s\" i 1 notitle w l lc rgb \"#77303030\"\n", (wdir+fname).c_str());
    fprintf(gpf, "unset multiplot\n");
    fprintf(gpf, "reset\n");
    fprintf(gpf, "%s\n", diverging_palette.c_str());
    fprintf(gpf, "set cbrange [%lf:%lf]\n", zmin, zmax);
    fprintf(gpf, "set zrange [%lf:%lf]\n", zmin, zmax);
    fprintf(gpf, "set yrange [*:*] reverse\n");
    fprintf(gpf, "unset label 11\n");
    fprintf(gpf, "set autoscale xfix\n");
    fprintf(gpf, "set autoscale yfix\n");
    fprintf(gpf, "set pm3d hidden3d 8 corners2color median\n");
    fprintf(gpf, "set term pngcairo dashed transparent enhanced size %d, %d font '%s,%d'  background rgb \"white\"\n",
        (int)lrint(width_in_pixels*2*grid_mer_fine.rows/double(grid_mer_fine.cols)),
        height_in_pixels_3d,
        #ifdef _WIN32
        "Verdana",
        #else
        "Arial",
        #endif
        fontsize3d
    );
    fprintf(gpf, "set output \"%sgrid_surface.png\"\n", wdir.c_str());
    fprintf(gpf, "unset xlab\n");
    fprintf(gpf, "unset ylab\n");
    if (img_filename.length() > 0) {
        fprintf(gpf, "set multiplot title \"%s\" font \",%d\"\n", img_filename.c_str(), fontsize);
        fprintf(gpf, "set tmargin 5\n");
    } else {
        fprintf(gpf, "set multiplot\n");
    }
    fprintf(gpf, "set ticslevel %lf\n", 0.0);
    fprintf(gpf, "set view 25, 350\n");
    fprintf(gpf, "set title \"Meridional MTF%2d (%s)\"\n", mtf_contrast, lpmm_mode ? "lp/mm" : "c/p");
    fprintf(gpf, "set size 1,0.5\n");   
    fprintf(gpf, "set origin 0.0,0.5\n");
    fprintf(gpf, "splot \"%s\" i 0 w pm3d lc rgb \"black\" lw %lg notitle\n",  (wdir+fname).c_str(), linewidth);
    fprintf(gpf, "set view 25, 350\n");
    fprintf(gpf, "set title \"Sagittal MTF%2d (%s)\"\n", mtf_contrast, lpmm_mode ? "lp/mm" : "c/p");
    fprintf(gpf, "set origin 0.0,0.0\n");
    fprintf(gpf, "splot \"%s\" i 1 w pm3d lc rgb \"black\" lw %lg notitle\n", (wdir+fname).c_str(), linewidth);
    fprintf(gpf, "unset multiplot\n");
    
    fclose(gpf);
    
    char* buffer = new char[1024];
    #ifdef _WIN32
    sprintf(buffer, "\"\"%s\" \"%sgrid.gnuplot\"\"", gnuplot_binary.c_str(), wdir.c_str());
    #else
    sprintf(buffer, "\"%s\" \"%sgrid.gnuplot\"", gnuplot_binary.c_str(), wdir.c_str());
    #endif
    int rval = system(buffer);
    if (rval != 0) {
        logger.error("Failed to execute gnuplot (error code %d)\n", rval);
        logger.info("You can try to execute [%s] to render the plots manually\n", buffer);
        gnuplot_failure = true;
    } else {
        logger.debug("%s\n", "Gnuplot plot completed successfully. Look for grid_image.png and grid_surface.png");
    }
    
    delete [] buffer;
    
}
    

    

