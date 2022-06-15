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

#include "include/mtf_renderer_focus.h"

void Mtf_renderer_focus::render(const vector<Mtf_profile_sample>& samples, Bayer::bayer_t bayer, 
    vector<Ellipse_detector>* ellipses, cv::Rect* dimension_correction) {
    Point2d centroid(0,0);
    
    cv::Mat merged = draw.rimg;
    initial_rows = draw.initial_rows;
    
    if (samples.size() < 100) {
        draw.fail_with_message(wdir + prname, string("Insufficient targets detected"));
        return;
    }
    
    if (!distance_scale.fiducials_found) {
        draw.fail_with_message(wdir + prname, string("Fiducials not found, probably incorrect chart type."));
        return;
    }
    
    double mean_y = 0;
    vector<Sample> data;
    double min_long = 1e50;
    double max_long = -1e50;
    double min_pix_long = 1e50;
    double max_pix_long = -1e50;
    for (size_t i=0; i < samples.size(); i++) {
        if (samples[i].quality > very_poor_quality) {
            Point2d d = samples[i].p - zero;
            Point2d coord(
                d.x*longitudinal.x + d.y*longitudinal.y,
                d.x*transverse.x + d.y*transverse.y
            );
            
            Point2d wc = distance_scale.estimate_world_coords(samples[i].p.x, samples[i].p.y);
            
            if (fabs(wc.y) < 20*psf && fabs(wc.x) < 180*psf) { 
            
                mean_y += coord.y;
                
                double depth = 0;
                distance_scale.estimate_depth_img_coords(samples[i].p.x, samples[i].p.y, depth);
                coord.x = depth;
                
                // coord is now projected into the surface
                data.push_back(Sample(coord.x, samples[i].mtf, 1.0, 1.0));
                min_long = min(coord.x, min_long);
                max_long = max(coord.x, max_long);
                min_pix_long = min(min_pix_long, coord.x);
                max_pix_long = max(max_pix_long, coord.x);
            }
            
            //merged.at<cv::Vec3b>(lrint(samples[i].p.y), lrint(samples[i].p.x)) = cv::Vec3b(200, 200, 0);
        }
    }
    mean_y /= double(data.size());
    vector<Point2d> img_cents;
    auto sliding(sliding_edges);
    for (size_t k=0; k < sliding.size(); k+= 2) {
        cv::Scalar edge_col(200,200,0);
        const std::pair<Point2d, Point2d>& ec = sliding[k];
        const std::pair<Point2d, Point2d>& en = sliding[k+1];
        Point2d cent = 0.5*(en.first + ec.first + ec.second);
        //cv::line(merged, cent - transverse, cent + transverse, edge_col, 2, cv::LINE_AA);
        img_cents.push_back(cent);
    }
    double min_ic_dist = 1e50;
    for (size_t i=0; i < img_cents.size() - 1; i++) {
        double dist = cv::norm(img_cents[i] - img_cents[i+1]);
        min_ic_dist = std::min(dist, min_ic_dist);
    }
    if (distance_scale.fiducial_code_set >= 0 && distance_scale.fiducial_code_set <= 4) {
        vector<vector<Point2d>> central_block_centroids {
            /*A0 */ {{411.913,983.093}, {411.844,969.094}, {411.773,954.867}, {411.702,940.407}, {411.629,925.707}, {411.555,910.761}, {411.479,895.564}, {411.403,880.109}, {411.325,864.389}, {411.246,848.398}, {411.165,832.128}, {411.083,815.572}, {411,798.722}, {410.914,781.571}, {410.828,764.11}, {410.74,746.331}, {410.65,728.225}, {410.559,709.784}, {410.466,690.996}, {410.371,671.854}, {410.274,652.346}, {410.176,632.462}, {410.075,612.192}, {409.973,591.522}, {409.868,570.443}, {409.762,548.941}, {409.653,527.004}, {409.542,504.618}, {409.429,481.77}, {409.313,458.444}, {409.195,434.626}, {409.075,410.3}, {408.952,385.45}, {408.826,360.058}, {408.697,334.107}, {408.566,307.578}, {408.431,280.451}, {408.294,252.706}, {408.153,224.321}, {408.009,195.275}, {407.862,165.543}, {407.711,135.102}, {407.556,103.925}, {407.398,71.9857}},
            /*A1 */ {{290.59,710.57}, {290.552,699.549}, {290.513,688.394}, {290.473,677.104}, {290.433,665.675}, {290.393,654.105}, {290.352,642.391}, {290.31,630.532}, {290.268,618.523}, {290.226,606.362}, {290.182,594.047}, {290.139,581.573}, {290.094,568.939}, {290.05,556.142}, {290.004,543.177}, {289.958,530.041}, {289.912,516.732}, {289.865,503.246}, {289.817,489.578}, {289.768,475.727}, {289.719,461.687}, {289.669,447.454}, {289.619,433.026}, {289.568,418.398}, {289.516,403.565}, {289.463,388.524}, {289.41,373.269}, {289.356,357.797}, {289.3,342.102}, {289.245,326.18}, {289.188,310.026}, {289.131,293.635}, {289.073,277.001}, {289.014,260.119}, {288.954,242.984}, {288.893,225.589}, {288.831,207.929}, {288.768,189.998}, {288.704,171.789}, {288.64,153.296}, {288.574,134.512}, {288.507,115.43}, {288.439,96.0432}, {288.37,76.3441}},
            /*A2 */ {{205.279,510.658}, {205.258,502.215}, {205.237,493.697}, {205.215,485.102}, {205.194,476.431}, {205.172,467.681}, {205.15,458.852}, {205.129,449.943}, {205.106,440.952}, {205.084,431.879}, {205.061,422.723}, {205.038,413.481}, {205.015,404.154}, {204.992,394.739}, {204.968,385.236}, {204.945,375.643}, {204.921,365.959}, {204.896,356.183}, {204.872,346.314}, {204.847,336.349}, {204.822,326.288}, {204.797,316.13}, {204.772,305.872}, {204.746,295.514}, {204.72,285.054}, {204.694,274.49}, {204.668,263.821}, {204.641,253.045}, {204.614,242.161}, {204.587,231.167}, {204.56,220.061}, {204.532,208.842}, {204.504,197.507}, {204.475,186.056}, {204.447,174.487}, {204.418,162.796}, {204.389,150.983}, {204.359,139.046}, {204.329,126.983}, {204.299,114.791}, {204.268,102.468}, {204.238,90.0135}, {204.207,77.4239}, {204.175,64.6974}},
            /*A3 */ {{145.06,365.679}, {145.049,359.347}, {145.038,352.975}, {145.026,346.561}, {145.015,340.106}, {145.004,333.608}, {144.992,327.068}, {144.981,320.485}, {144.969,313.859}, {144.958,307.189}, {144.946,300.474}, {144.934,293.715}, {144.922,286.911}, {144.91,280.061}, {144.898,273.165}, {144.886,266.222}, {144.874,259.233}, {144.861,252.196}, {144.849,245.111}, {144.836,237.977}, {144.824,230.795}, {144.811,223.563}, {144.798,216.282}, {144.786,208.949}, {144.773,201.566}, {144.76,194.131}, {144.747,186.644}, {144.733,179.105}, {144.72,171.512}, {144.707,163.866}, {144.693,156.165}, {144.68,148.409}, {144.666,140.598}, {144.652,132.731}, {144.638,124.806}, {144.624,116.825}, {144.61,108.786}, {144.596,100.688}, {144.582,92.5307}, {144.567,84.3136}, {144.553,76.0361}, {144.538,67.6973}, {144.524,59.2966}, {144.509,50.8334}},
            /*A4 */ {{102.514,260.994}, {102.508,256.323}, {102.503,251.63}, {102.497,246.914}, {102.491,242.177}, {102.485,237.418}, {102.479,232.636}, {102.473,227.831}, {102.467,223.004}, {102.461,218.154}, {102.455,213.281}, {102.449,208.385}, {102.443,203.465}, {102.437,198.522}, {102.431,193.555}, {102.424,188.564}, {102.418,183.549}, {102.412,178.51}, {102.406,173.447}, {102.4,168.359}, {102.393,163.246}, {102.387,158.108}, {102.38,152.945}, {102.374,147.757}, {102.368,142.544}, {102.361,137.304}, {102.355,132.039}, {102.348,126.748}, {102.341,121.43}, {102.335,116.086}, {102.328,110.716}, {102.322,105.318}, {102.315,99.8936}, {102.308,94.4418}, {102.301,88.9625}, {102.294,83.4556}, {102.288,77.9208}, {102.281,72.3579}, {102.274,66.7667}, {102.267,61.1469}, {102.26,55.4984}, {102.253,49.821}, {102.246,44.1143}, {102.239,38.3783}}
        };
        vector<Point2d>& block_cents = central_block_centroids[distance_scale.fiducial_code_set];
        for (size_t k=0; k < block_cents.size(); k++) {
            cv::Scalar edge_col(20,20,200);
            double world_y = block_cents[k].x - distance_scale.page_scale_factor*420/sqrt(2.0)*0.5;
            double world_x = block_cents[k].y - distance_scale.page_scale_factor*420*0.5;
            Point2d cent = distance_scale.world_to_image(world_x, world_y);
            
            double min_dist = 1e50;
            for (auto& ic: img_cents) {
                min_dist = std::min(cv::norm(ic - cent), min_dist);
            }
            const double cross_w = 0.5*min_ic_dist - std::max(0.1*min_ic_dist, 2.0);
            const double cross_l = min_ic_dist;
            if (min_dist > 1.01*min_ic_dist) {
                cv::line(merged, cent - cross_l*transverse - cross_w*longitudinal, cent + cross_l*transverse + cross_w*longitudinal, edge_col, 2, cv::LINE_AA);
                cv::line(merged, cent - cross_l*transverse + cross_w*longitudinal, cent + cross_l*transverse - cross_w*longitudinal, edge_col, 2, cv::LINE_AA);
            }
        }
    } else {
        logger.error("mtf_renderer_focus: encountered unknown fiducial code %d\n", distance_scale.fiducial_code_set);
    }
    
    const int sh = 4;
    const double sgw[] = {-21/231.0, 14/231.0, 39/231.0, 54/231.0, 59/231.0, 54/231.0, 39/231.0, 14/231.0, -21/231.0};
    sort(data.begin(), data.end());
    
    // just pretend our samples are equally spaced
    vector<Sample> ndata;
    for (size_t i=sh; i < data.size() - sh; i++) {
        double val = 0;
        for (int w=-sh; w <= sh; w++) {
            val += sgw[w+sh] * data[i+w].y;    
        }
        
        
        if (fabs(val - data[i].y)/val < 0.04) {
            ndata.push_back(data[i]);
        }
        
    }
    logger.debug("dropped %lu samples, %lu remain\n", data.size() - ndata.size(), ndata.size());
    data = ndata;
    
    Ratpoly_fit cf(data, 4, 2);
    VectorXd sol = rpfit(cf);
    
    // perform a few steps of IRLS
    double prev_err = 1e50;
    int dccount;
    for (int iter=0; iter < 50; iter++) {
        double errsum = 0;
        double wsum = 0;
        dccount = 0;
        for (size_t k=0; k < data.size(); k++) {
            double y = cf.rpeval(sol, cf.scale(data[k].x))/cf.ysf;
            double e = fabs(y - data[k].y);
            errsum += e * data[k].yweight;
            wsum += data[k].yweight;
            data[k].yweight = 1.0; 
            if (e/y > 0.05) { // kill really poor outliers
                data[k].yweight = 0;
                dccount++;
            }
            
        }
        errsum /= wsum;
        logger.debug("iter %d err: %lg dc=%lg\n", iter, errsum, dccount / double(data.size()));
        cf.order_m = 2; // reset the denominator order, it will be reduced if necessary
        sol = rpfit(cf);
        if (iter > 0 && (prev_err - errsum)/prev_err < 0.0001) {
            logger.debug("bailing out at iter %d\n", iter);
            break;
        }
        prev_err = errsum;
    }
    double errsum = 0;
    double wsum = 0;
    for (size_t k=0; k < data.size(); k++) {
        double y = cf.rpeval(sol, cf.scale(data[k].x))/cf.ysf;
        double e = fabs(y - data[k].y);
        errsum += e * data[k].yweight;
        wsum += data[k].yweight;
    }
    errsum /= wsum;
    logger.debug("final model fit error (weighted): %lg\n", errsum);
    
    // do a quick bootstrap to estimate some bounds on the focus peak
    vector<double> mc_peaks;
    const int n_mc = 30;
    for (int i=0; i < n_mc; i++) {
        vector<Sample> mc_data;
        for (size_t j=i; j < data.size(); j += 10) {
            mc_data.push_back(data[j]);
        }
        Ratpoly_fit mc_cf(mc_data, 4, 2, true);
        VectorXd mc_sol = rpfit(mc_cf);
        
        double mc_peak = mc_cf.peak(mc_sol);
        mc_peaks.push_back(mc_peak);
    }
    sort(mc_peaks.begin(), mc_peaks.end());
    double mc_p5 = mc_peaks[lrint(n_mc*0.05)];
    double mc_p95 = mc_peaks[lrint(n_mc*0.95)];
    
    string rp_name = wdir + ((bayer != Bayer::NONE) ? Bayer::to_string(bayer) + "_" : "") + "profile_points.txt";
    FILE* fraw = fopen(rp_name.c_str(), "wt");
    for (size_t i=0; i < data.size(); i++) {
        fprintf(fraw, "%lf %lf\n", data[i].x, data[i].y);
    }
    fclose(fraw);
    
    // now we can plot the reconstruction?
    double peak_y = 0;
    string cp_name = wdir + ((bayer != Bayer::NONE) ? Bayer::to_string(bayer) + "_" : "") + "profile_curve.txt";
    FILE* fout = fopen(cp_name.c_str(), "wt");
    for (double x=min_long; x < max_long; x += 1) {
        double y = cf.rpeval(sol, cf.scale(x))/cf.ysf;
        fprintf(fout, "%lf %lf\n", x, y);
        
        if (y > peak_y) {
            peak_y = y;
        }
    }
    fclose(fout);
    
    double focus_peak = cf.peak(sol); // works on pre-stretched depth
    logger.debug("focus_plane %lg\n", focus_peak);
    
    draw.chart_centre_marker();
    draw.camera_centre_marker(
        dimension_correction->width/2 - dimension_correction->x, 
        dimension_correction->height/2 - dimension_correction->y
    );
    
    if (ellipses) {
        cv::Scalar ellipse_col(200,200,0);
        cv::Scalar vec_col(50,50,200);
        constexpr double error_tolerance = 1.5; // allow roughly 1 pixel error in each axis, with a little bit of a tolerance
        for (auto e: *ellipses) {
            if (!e.valid) continue;
            
            int fid_idx = -1;
            for (size_t li=0; li < distance_scale.fid_img_points.size(); li++) {
                if (cv::norm(distance_scale.fid_img_points[li] - cv::Point2d(e.centroid_x, e.centroid_y)) < 0.001) {
                    fid_idx = li;
                }
            }
            Point2d fid_offset(0, 0);
            if (fid_idx >= 0) {
                cv::Point3d& p3d = distance_scale.fid_world_points[fid_idx];
                Point2d backproj_pt = distance_scale.world_to_image(p3d.x, p3d.y, p3d.z);
                fid_offset = backproj_pt - distance_scale.fid_img_points[fid_idx];
            }
            if (cv::norm(fid_offset) < error_tolerance) { 
                fid_idx = -1; // skip this fiducial since the error is small enough
            }
            
            Point2d prev(0,0);
            for (double theta=0; theta < 2*M_PI; theta += M_PI/64.0) {
                double synth_x = e.major_axis * cos(theta);
                double synth_y = e.minor_axis * sin(theta);
                double rot_x = cos(e.angle)*synth_x - sin(e.angle)*synth_y + e.centroid_x;
                double rot_y = sin(e.angle)*synth_x + cos(e.angle)*synth_y + e.centroid_y;

                // clip to image size, just in case
                rot_x = max(rot_x, 0.0);
                rot_x = min(rot_x, (double)(merged.cols-1));
                rot_y = max(rot_y, 0.0);
                rot_y = min(rot_y, (double)(merged.rows-1));
                
                Point2d current(rot_x, rot_y);
                if (theta > 0) {
                    cv::line(merged, prev, current, ellipse_col, 1, cv::LINE_AA);
                    
                    if (fid_idx >= 0) {
                        cv::line(merged, prev + fid_offset, current + fid_offset, vec_col, 1, cv::LINE_AA);
                    }
                }
                
                prev = current;
            }
        }
        
        for (size_t li=0; li < distance_scale.fid_img_points.size(); li++) {
            cv::Point3d& p3d = distance_scale.fid_world_points[li];
            cv::Point2d& p2d = distance_scale.fid_img_points[li];
            
            Point2d backproj_pt = distance_scale.world_to_image(p3d.x, p3d.y, p3d.z);
            
            if (cv::norm(p2d - backproj_pt) >= error_tolerance) { 
                cv::line(merged, p2d, backproj_pt, vec_col, 1, cv::LINE_AA);
            }
        }
    }
    
    vector<Point2d> curve;
    min_pix_long = -merged.cols/2;
    max_pix_long = merged.cols/2;
    double mtf_peak_value = 0;
    double peak_wx = 0;
    // main MTF profile curve
    for (double x=min_pix_long; x < max_pix_long; x += 1) {
        double px = longitudinal.x * x  + longitudinal.y * mean_y + zero.x;
        double py = transverse.x * x + transverse.y *  mean_y + zero.y;
        
        double depth = 0;
        distance_scale.estimate_depth_img_coords(px, py, depth);
        
        if (depth > min_long && depth < max_long) {
            double mtf = cf.rpeval(sol, cf.scale(depth))/cf.ysf;
            
            double world_y = psf*((mtf / peak_y) * 130 - 130);
                
            Point2d wc = distance_scale.estimate_world_coords(px, py);
            
            if (mtf > mtf_peak_value) {
                mtf_peak_value = mtf;
                peak_wx = wc.x;
            }
                
            Point2d proj_p = distance_scale.world_to_image(wc.x, world_y);
            curve.push_back(proj_p);
        }
    }
    if (cf.has_poles(sol)) {
        draw.fail_circle();
    } 
    draw.curve(curve, cv::Scalar(128, 128, 128), 6);
    draw.curve(curve, cv::Scalar(40, 90, 40), 3, cv::Scalar(40, 255, 40));
    
    curve.clear();
    for (double wy=-130*psf; wy < 130*psf; wy += 1) {
        Point2d proj_p = distance_scale.world_to_image(peak_wx, wy);
        curve.push_back(proj_p);
    }
    draw.curve(curve, cv::Scalar(255, 30, 30), 3);
    
    
    cv::Scalar axisdark = cv::Scalar(30, 140, 190);
    cv::Scalar axislight = cv::Scalar(40, 187, 255);
    
    curve.clear();
    for (double ystep=-135*psf; ystep <= 135*psf; ystep += 2) {
        curve.push_back(distance_scale.world_to_image(-180*psf, ystep));
    }
    draw.curve(curve, axisdark, 2, axisdark);
    curve.clear();
    for (double xstep=-180*psf; xstep <= 180*psf; xstep += 2) {
        curve.push_back(distance_scale.world_to_image(xstep, -135*psf));
    }
    draw.curve(curve, axisdark, 2, axislight);
    curve.clear();
    curve.push_back(distance_scale.world_to_image(-180*psf, -135*psf));
    
    // find a reasonable value for zmax that remains inside the image
    double zmax = 0;
    const double border = 35;
    for (; zmax > -200; zmax -= 2) {
        Point2d p = distance_scale.world_to_image(-180*psf, -135*psf, zmax*psf);
        if (p.x < border || p.y < border || p.x > img.cols - 1 - border || p.y > initial_rows - 1 - border) break;
    }
    
    curve.push_back(distance_scale.world_to_image(-180*psf, -135*psf, zmax*psf));
    draw.curve(curve, axisdark, 2, axisdark);
    
    int font = cv::FONT_HERSHEY_DUPLEX; 
    char tbuffer[1024];
    
    draw.text(-178*psf,  135*psf, 0, axisdark, "Y");
    draw.text( 180*psf, -133*psf, 0, axislight, "X");
    draw.text(-180*psf, -132*psf, zmax*psf, axisdark, "Z");
    
    cv::Scalar resultcolour(255, 255, 30);
    sprintf(tbuffer, "%.3lf (c/p)", mtf_peak_value);
    int tbaseline=0;
    cv::Size tsize = cv::getTextSize(tbuffer, font, 1.0, 2, &tbaseline);
    Point2d textpos = distance_scale.world_to_image(peak_wx, -20*psf);
    textpos -= Point2d(tsize.width/2.0, 0);
    cv::putText(merged, tbuffer, textpos, font, 1, CV_RGB(20, 20, 20), 2, cv::LINE_AA);
    cv::putText(merged, tbuffer, textpos, font, 1, resultcolour, 1, cv::LINE_AA);
    double prev_line_height = tsize.height;
    
    sprintf(tbuffer, "%.1lf mm", focus_peak);
    tsize = cv::getTextSize(tbuffer, font, 1.0, 2, &tbaseline);
    textpos = distance_scale.world_to_image(peak_wx, -20*psf);
    textpos -= Point2d(tsize.width/2.0, -prev_line_height*1.5);
    cv::putText(merged, tbuffer, textpos, font, 1, CV_RGB(20, 20, 20), 2, cv::LINE_AA);
    cv::putText(merged, tbuffer, textpos, font, 1, resultcolour, 1, cv::LINE_AA);
    
    // blank out the text region (again)
    rectangle(merged, Point2d(0, img.rows), Point2d(merged.cols, merged.rows), cv::Scalar::all(255), cv::FILLED);
    
    cv::Scalar red(30, 30, 200);
    cv::Scalar yellow(40, 187, 255);
    cv::Scalar black(0, 0, 0);
    cv::Scalar green(30, 200, 30);
    
    sprintf(tbuffer, "Focus peak at depth %.1lf mm [%.1lf,%.1lf] relative to chart origin.",focus_peak, mc_p5, mc_p95);
    int baseline = 0;
    cv::Size ts = cv::getTextSize(tbuffer, font, 1, 1, &baseline);
    draw.text(Point2d(50, initial_rows + ts.height*1.75), black, "%s", tbuffer);
    
    draw.text(Point2d(50, initial_rows + ts.height*2*1.75), black, "Estimated chart distance=%.2lf mm.", distance_scale.centre_depth);
    
    
    if (distance_scale.bundle_rmse >= 1.5) {
        sprintf(tbuffer, "Test chart not flat? Large geom. calib. error of %.1lf pixels", distance_scale.bundle_rmse);
    } else {
        sprintf(tbuffer, "Bundle adjustment RMSE=%.3lf pixels (ideal < 1)", distance_scale.bundle_rmse);
    }
    draw.text(Point2d(50, initial_rows + ts.height*3*1.75), black, "%s", tbuffer);
    
    ts = cv::getTextSize(tbuffer, font, 1, 1, &baseline);
    double col2 = ts.width + 150;
    cv::Scalar rmse_col = green;
    if (distance_scale.bundle_rmse < 1.5) {
        if (distance_scale.bundle_rmse > 0.7) {
            rmse_col = yellow;
        }
        draw.checkmark(Point2d(25, initial_rows + ts.height*3*1.75), rmse_col);
    } else {
        draw.crossmark(Point2d(35, initial_rows + ts.height*2.75*1.75), red);
    }
    
    double rpy = initial_rows + ts.height*4*1.75;
    double rpx = 50;
    draw.text(Point2d(rpx, rpy), black, "Chart z-angle=%.1lf degrees (ideal = 45)", distance_scale.get_normal_angle_z());
    
    cv::Scalar zang_col = green;
    if (fabs(distance_scale.get_normal_angle_z() - 45) < 10) {
        if (fabs(distance_scale.get_normal_angle_z() - 45) > 5) {
            zang_col = yellow;
        }
        draw.checkmark(Point2d(25, rpy), zang_col);
    } else {
        draw.crossmark(Point2d(35, rpy - 0.25*1.75*ts.height), red);
    }
    
    rpy = initial_rows + ts.height*5*1.75;
    draw.text(Point2d(rpx, rpy), black, "Chart y-angle=%.1lf degrees (ideal = 0)", distance_scale.get_normal_angle_y());
    
    cv::Scalar yang_col = green;
    if (fabs(distance_scale.get_normal_angle_y()) < 2) {
        if (fabs(distance_scale.get_normal_angle_y()) > 1) {
            yang_col = yellow;
        }
        draw.checkmark(Point2d(25, rpy), yang_col);
    } else {
        draw.crossmark(Point2d(35, rpy - 0.25*1.75*ts.height), red);
    }
    
    
    double white_clip = 0;
    double black_clip = 0;
    double overexposure = 0;
    exposure_checks(Point2d(180*psf, 135*psf), white_clip, black_clip, overexposure);
    
    rpx = col2;
    rpy = initial_rows + ts.height*3*1.75;
    if (white_clip > 5){
        sprintf(tbuffer, "Probable highlight clipping, severity=%.1lf%%", white_clip);
    } else {
        sprintf(tbuffer, "No highlight clipping detected");
    }
    draw.text(Point2d(rpx, rpy), black, "%s", tbuffer);
    
    cv::Scalar wc_col = green;
    if (white_clip < 10) {
        draw.checkmark(Point2d(rpx-25, rpy), wc_col);
    } else {
        draw.crossmark(Point2d(rpx-15, rpy - 0.25*1.75*ts.height), red);
    }
    
    rpy = initial_rows + ts.height*4*1.75;
    if (black_clip > 5){
        sprintf(tbuffer, "Probable shadow clipping, severity=%.1lf%%", black_clip);
    } else {
        sprintf(tbuffer, "No shadow clipping detected");
    }
    draw.text(Point2d(rpx, rpy), black, "%s", tbuffer);
    
    cv::Scalar bc_col = green;
    if (black_clip < 10) {
        draw.checkmark(Point2d(rpx-25, rpy), bc_col);
    } else {
        draw.crossmark(Point2d(rpx-15, rpy - 0.25*1.75*ts.height), red);
    }
    
    
    rpy = initial_rows + ts.height*5*1.75;
    if (overexposure > 0){
        draw.text(Point2d(rpx, rpy), black, "Probable overexposure, severity=%.1lf%%", overexposure);
        draw.crossmark(Point2d(rpx-15, rpy - 0.25*1.75*ts.height), red);
    } 
    
    
    imwrite(wdir + prname, merged);
    
}

VectorXd Mtf_renderer_focus::rpfit(Ratpoly_fit& cf, bool scale, bool refine) {
    const vector<Sample>& pts_row = cf.get_data();
    
    if (scale) {
        double xmin = 1e50;
        double xmax = -1e50;
        double ysf=0;
        for (size_t i=0; i < pts_row.size(); i++) {
            xmin = std::min(xmin, pts_row[i].x);
            xmax = std::max(xmax, pts_row[i].x);
            ysf = max(ysf, fabs(pts_row[i].y));
        }
        cf.xs_min = 0.5*(xmin + xmax);
        cf.xs_scale = 2.0/(xmax - xmin);
        cf.ysf = ysf = 1.0/ysf;
    }
    
    VectorXd sol;
    bool done = false;
    
    while (!done) {
        int tdim = cf.dimension();
        MatrixXd cov = MatrixXd::Zero(tdim, tdim);
        VectorXd b = VectorXd::Zero(tdim);
        VectorXd a = VectorXd::Zero(tdim);
        
        
        
        for (int iter=0; iter < 1; iter++) {
            cov.setZero();
            b.setZero();
            a.setZero();
            
            for (size_t i=0; i < pts_row.size(); i++) {
                const Sample& sp = pts_row[i];
                double w = sp.weight * sp.yweight;
                a[0] = 1*w;
                double prod = cf.scale(sp.x); // top poly
                for (int j=1; j <= cf.order_n; j++) {
                    a[j] = w*cf.cheb(j, prod);
                }
                // bottom poly
                for (int j=1; j <= cf.order_m; j++) {
                    a[j+cf.order_n] = cf.cheb(j, prod)*w*sp.y*cf.ysf;
                }
                
                for (int col=0; col < tdim; col++) { 
                    for (int icol=0; icol < tdim; icol++) { // build covariance of design matrix : A'*A
                        cov(col, icol) += a[col]*a[icol];
                    }
                    b[col] += cf.base_value*a[col]*sp.y*cf.ysf*w; // build rhs of system : A'*b
                }
            }
            
            for (int col=cf.order_n+1; col < cov.cols(); col++) {
                cov.col(col) = -cov.col(col);
            }
            
            sol = cov.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
        }
        
        // now perform non-linear optimization
        
        if (refine) {
            sol = cf.gauss_newton_armijo(sol);
        }
        
        if (cf.has_poles(sol)) {
            cf.order_m--;
            if (!cf.silent) {
                logger.debug("reduced order: (%d, %d)\n", cf.order_n, cf.order_m);
            }
            done = false;
        } else {
            done = true;
        }
    }
    
    return sol;
}
    
void Mtf_renderer_focus::exposure_checks(const Point2d& dims, double& white_clip, double& black_clip, double& overexposure) {
    
    vector<Point2d> corners { {-dims.x, -dims.y}, {dims.x, -dims.y}, {dims.x, dims.y}, {-dims.x, dims.y}};
    vector<cv::Point> roi;
    for (const auto& p: corners) {
        Point2d v = distance_scale.world_to_image(p.x, p.y);
        roi.push_back(cv::Point(v.x, v.y));
    }
    
    cv::Mat mask(img.rows, img.cols, CV_8UC1, cv::Scalar::all(0));
    cv::fillConvexPoly(mask, (const cv::Point*)roi.data(), roi.size(), cv::Scalar::all(255));
    
    // we could (in theory) mask out all the objects we used (ellipses and blocks)
    // that should take care of most of the gradients?
    
    cv::Rect bounds=cv::boundingRect(roi);
    
    bounds.x = max(0, bounds.x);
    bounds.y = max(0, bounds.y);
    bounds.width = min((img.cols - 1 - bounds.x), bounds.width);
    bounds.height = min((img.rows - 1 - bounds.y), bounds.height);
    
    //cv::Mat masked;
    //img.copyTo(masked, mask);
    
    // calculate histogram
    vector<uint32_t> histo(65536 + 1, 0);
    for (int row=bounds.y; row < (bounds.y + bounds.height); row++) {
        for (int col=bounds.x; col < (bounds.x + bounds.width); col++) {
            int val = img.at<uint16_t>(row, col);
            if (mask.at<uint8_t>(row, col) > 0) {
                histo[val]++;
            }
        }
    }
    int last_nz = histo.size()-1;
    while (last_nz > 0 && histo[last_nz] == 0) last_nz--;
    
    #ifdef MDEBUG
    FILE* hfile = fopen("focus_histo.txt", "wt");
    size_t ch = 0;
    for (size_t j=0; j <= size_t(last_nz); j++) {
        if (histo[j] > 0) {
            ch += histo[j];
            fprintf(hfile, "%lu %u %lu\n", j, histo[j], ch);
        }
    }
    fclose(hfile);
    #endif
    
    vector<double> chisto(histo.size(), 0.0);
    chisto[0] = histo[0];
    for (size_t j=1; j < chisto.size(); j++) {
        chisto[j] = chisto[j-1] + histo[j];
    }
    
    // catch a histogram with an abnormally large count in the last bin
    white_clip = 0;
    if (chisto.back() - chisto[last_nz-1] > 0.05*chisto.back()) {
        white_clip = (chisto[last_nz] - chisto[last_nz-1])/chisto.back() * 100;
    }
    
    // catch a histogram with a bin within 5% of end with an abnormally large count
    int wc_sentinel = max(2.0, last_nz * 0.95);
    int wc_peak = last_nz;
    while (wc_peak > wc_sentinel) {
        if (chisto[wc_peak] - chisto[wc_peak-1] > 0.05*chisto.back()) {
            white_clip = (chisto[wc_peak] - chisto[wc_peak-1])/chisto.back() * 100;
        }
        wc_peak--;
    }
    
    overexposure = 0;
    // catch a histogram that has a > 70% weight in last 10% ....
    if ( (chisto.back() - chisto[65535*0.9]) > 0.7*histo.back() ) {
        overexposure = (chisto.back() - chisto[65535*0.9]) / (0.5*chisto.back()) * 100;
        overexposure = min(overexposure, 100.0);
    }
    
    int first_nz = 0;
    while (first_nz < (int)histo.size()-1 && histo[first_nz] == 0) first_nz++;
    
    vector<double> rhisto(histo.size(), 0.0);
    rhisto.back() = histo.back();
    for (size_t j=histo.size()-2; j > 0; j--) {
        rhisto[j] = rhisto[j+1] + histo[j];
    }
    rhisto[0] = rhisto[1] + histo[0];
    
    int p5 = first_nz+1;
    while (p5 < (int)histo.size()-1 && rhisto[p5] > 0.05*rhisto.front()) p5++;
    
    black_clip = 0;
    if (rhisto[first_nz] - rhisto[first_nz+1] > 0.0001*rhisto.front()) {
        black_clip = (rhisto[first_nz] - rhisto[first_nz+1])/(0.05*rhisto.front()) * 100;
        black_clip = min(black_clip, 100.0);
    }
    
    // compute integral image
    //cv::Mat temp;
    //masked.convertTo(temp, CV_32F);
    //integral(temp, masked, CV_32F);
    
    // homomorphic filtering tells us that we can subtract the log 
    // of the blurred version to obtain the "flat" lighting version
    // .... so the log of the blur is the slow-varying part
    // we want to see how uniform that is ...?
    
    // edges/nodata will be a problem (but we can construct a second mask?)
    
    
    // maybe we can simply blur a lot, and look at the
    // horizontal / vertical projections ? These should be flat?
}    