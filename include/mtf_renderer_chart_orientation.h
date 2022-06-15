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
#ifndef MTF_RENDERER_CHART_ORIENTATION_H
#define MTF_RENDERER_CHART_ORIENTATION_H

#include "include/logger.h"
#include "mtf_renderer.h"
#include "common_types.h"
#include "distance_scale.h"
#include "include/camera_draw.h"

class Mtf_renderer_chart_orientation : public Mtf_renderer {
  public:
    Mtf_renderer_chart_orientation(
        const std::string& img_filename,
        const std::string& wdir, 
        const std::string& co_fname, 
        [[maybe_unused]] const cv::Mat& img,
        [[maybe_unused]] int gnuplot_width,
        Distance_scale& distance_scale,
        cv::Rect* dimension_correction = NULL)
      :  Mtf_renderer(img_filename),
         wdir(wdir), co_fname(co_fname), 
         dimension_correction(dimension_correction),
         draw(img, distance_scale, (gnuplot_width / double(img.cols)), 50) {
    
      
    }
    
    void render(const vector<Block>& ) {
    
        Distance_scale& distance_scale = draw.distance_scale;
        
        if (!distance_scale.fiducials_found) {
            draw.fail_with_message(wdir + '/' + co_fname, string("Fiducials not found, probably incorrect chart type."));
            return;
        }
        
        draw.chart_centre_marker();        
        
        // draw centre-of-camera marker
        draw.camera_centre_marker(
            draw.scaling_factor*(dimension_correction->width/2 - dimension_correction->x), 
            draw.scaling_factor*(dimension_correction->height/2 - dimension_correction->y)
        );
        
        vector<Point2d> curve;
        
        cv::Scalar c_dark(20, 20, 20);
        cv::Scalar c_red(30, 30, 255);
        cv::Scalar c_green(30, 255, 30);
        cv::Scalar c_blue(255, 30, 30);
        
        cv::Scalar c_lred(80, 80, 255);
        cv::Scalar c_lgreen(130, 255, 130);
        cv::Scalar c_lblue(255, 182, 0);
        
        curve.clear();
        const double border = 35;
        double ymax = 0;
        for (; ymax > -100; ymax -= 2) {
            Point2d p = distance_scale.world_to_image(-draw.psf, ymax*draw.psf, 0);
            if (p.x < border || p.y < border || p.x > draw.rimg.cols - 1 - border || p.y > draw.rimg.rows - 1 - border) break;
        }
        
        double xmax = 0;
        for (; xmax > -100; xmax -= 2) {
            Point2d p = distance_scale.world_to_image(xmax*draw.psf, ymax*draw.psf, 0);
            if (p.x < border || p.y < border || p.x > draw.rimg.cols - 1 - border || p.y > draw.rimg.rows - 1 - border) break;
        }
        
        double zmax = 0;
        for (; zmax > -100; zmax -= 2) {
            Point2d p = distance_scale.world_to_image(xmax*draw.psf, ymax*draw.psf, zmax*draw.psf);
            if (p.x < border || p.y < border || p.x > draw.rimg.cols - 1 - border || p.y > draw.rimg.rows - 1 - border) break;
        }
        
        for (double ystep=ymax; ystep <= 0; ystep += 2) {
            curve.push_back(distance_scale.world_to_image(xmax*draw.psf, ystep*draw.psf));
        }
        draw.curve(curve, c_dark, 3, c_dark);
        draw.curve(curve, c_green, 2, c_green);
        
        // in landscape orientation, dy > dy on the y axis
        bool landscape = fabs(curve.front().y - curve.back().y) > fabs(curve.front().x - curve.back().x);
        bool flipped = false;
        if (landscape) {
            flipped = curve.front().y < curve.back().y;
        } else {
            flipped = curve.front().x > curve.back().x;
        }
        
        curve.clear();
        for (double xstep=xmax; xstep <= 0; xstep += 2) {
            curve.push_back(distance_scale.world_to_image(xstep*draw.psf, ymax*draw.psf));
        }
        draw.curve(curve, c_dark, 3, c_dark);
        draw.curve(curve, c_red, 2, c_red);
        
        curve.clear();
        
        curve.push_back(distance_scale.world_to_image(xmax*draw.psf, ymax*draw.psf));
        curve.push_back(distance_scale.world_to_image(xmax*draw.psf, ymax*draw.psf, zmax*draw.psf));
        draw.curve(curve, c_dark, 3, c_dark);
        draw.curve(curve, c_lblue, 2, c_lblue);
        
        // if roll angle is close to 90 degrees, then assume we are in portrait mode
        // and reduce angles appropriately
        
        double e_roll = distance_scale.roll_angle/M_PI*180;
        double e_pitch = distance_scale.pitch_angle/M_PI*180;
        double e_yaw = distance_scale.yaw_angle/M_PI*180;
        if (fabs(e_roll - 90) < fabs(e_roll)) {
            e_roll -= 90;
            e_pitch = distance_scale.yaw_angle/M_PI*180;
            e_yaw   = -distance_scale.pitch_angle/M_PI*180;
        } else {
            if (fabs(e_roll + 90) < fabs(e_roll)) {
                e_roll += 90;
                e_pitch = distance_scale.yaw_angle/M_PI*180;
                e_yaw   = -distance_scale.pitch_angle/M_PI*180;
            }
        }
        
        draw.arc_with_arrow(0, 10, fabs(xmax), c_lred,   (e_pitch < 0) ^ flipped);
        draw.arc_with_arrow(1, 10, fabs(ymax), c_lgreen, (e_yaw < 0) ^ flipped);
        draw.arc_with_arrow(2, 10, fabs(zmax), c_lblue,  (e_roll < 0) ^ flipped);

        
        if (flipped) {
            if (landscape) {
                draw.text_block_ra(xmax*draw.psf, ymax*draw.psf, zmax*draw.psf, c_lblue, "Roll=%.2lf", fabs(e_roll));
                draw.text_block_ra(xmax*draw.psf, 15*draw.psf, 0, c_lgreen, "Yaw=%.2lf", fabs(e_yaw));
                draw.text_block(10*draw.psf, (ymax+20)*draw.psf, 0, c_lred, "Pitch=%.2lf", fabs(e_pitch));
            } else {
                draw.text_block_ra(xmax*draw.psf, (ymax+20)*draw.psf, zmax*draw.psf, c_lblue, "Roll=%.2lf", fabs(e_roll));
                draw.text_block((xmax+15)*draw.psf, 10*draw.psf, 0, c_lgreen, "Yaw=%.2lf", fabs(e_yaw));
                draw.text_block_ra(10*draw.psf, ymax*draw.psf, 0, c_lred, "Pitch=%.2lf", fabs(e_pitch));
            }
        } else {
            if (landscape) {
                draw.text_block((xmax+20)*draw.psf, ymax*draw.psf, zmax*draw.psf, c_lblue, "Roll=%.2lf", fabs(e_roll));
                draw.text_block(xmax*draw.psf, 10*draw.psf, 0, c_lgreen, "Yaw=%.2lf", fabs(e_yaw));
                draw.text_block_ra(10*draw.psf, (ymax+20)*draw.psf, 0, c_lred, "Pitch=%.2lf", fabs(e_pitch));
            } else {
                draw.text_block(xmax*draw.psf, (ymax+20)*draw.psf, zmax*draw.psf, c_lblue, "Roll=%.2lf", fabs(e_roll));
                draw.text_block_ra((xmax+20)*draw.psf, 10*draw.psf, 0, c_lgreen, "Yaw=%.2lf", fabs(e_yaw));
                draw.text_block(10*draw.psf, ymax*draw.psf, 0, c_lred, "Pitch=%.2lf", fabs(e_pitch));
            }
        }
        
        cv::Scalar black(0, 0, 0);
        
        // ugly, but we have to obtain the text height somehow
        char tbuffer[] = "dummy";
        int baseline;
        const int font = cv::FONT_HERSHEY_DUPLEX; 
        cv::Size ts = cv::getTextSize(tbuffer, font, 1, 1, &baseline);
        
        if (distance_scale.user_provided_focal_ratio > 0) {
            draw.checkmark(Point2d(25, draw.initial_rows + ts.height*1*1.75), c_green);
            draw.text(Point2d(50, draw.initial_rows + ts.height*1*1.75), black, "EXIF or user-provided focal ratio: %.2lf", distance_scale.user_provided_focal_ratio);
        } else {
            draw.crossmark(Point2d(35, draw.initial_rows + ts.height*0.75*1.75), c_red);
            draw.text(Point2d(50, draw.initial_rows + ts.height*1*1.75), black, "Estimated focal ratio: %.2lf", distance_scale.focal_length);
        }
        
        
        imwrite(wdir + '/' + co_fname, draw.rimg);
    }
    
  private:
    string wdir;
    string co_fname; 
    
    cv::Rect* dimension_correction;
    
    Camera_draw draw;
};

#endif
