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
#ifndef MTF_RENDERER_FOCUS_H
#define MTF_RENDERER_FOCUS_H

#include "include/logger.h"
#include "mtf_renderer.h"
#include "common_types.h"
#include "loess_fit.h"
#include "mtf_profile_sample.h"

#include <stdlib.h>
#include <vector>
using std::vector;
using std::pair;
using std::make_pair;

#include "focus_surface.h"
#include "distance_scale.h"
#include "mtf50_edge_quality_rating.h"
#include "bayer.h"
#include "ellipse.h"
#include "camera_draw.h"

class Mtf_renderer_focus : public Mtf_renderer {
  public:
    Mtf_renderer_focus(Distance_scale& distance_scale,
        const vector<std::pair<Point2d, Point2d>>& sliding_edges,
        const std::string& wdir, const std::string& prof_fname, 
        const cv::Mat& img, 
        [[maybe_unused]] bool lpmm_mode=false, [[maybe_unused]] double pixel_size=1.0) 
      :  zero(distance_scale.zero), 
         transverse(distance_scale.transverse), 
         longitudinal(distance_scale.longitudinal),
         wdir(wdir), prname(prof_fname),
         img(img), 
         distance_scale(distance_scale),
         sliding_edges(sliding_edges),
         psf(distance_scale.page_scale_factor),
         draw(img, distance_scale, 1.0, 220) {
         
    }
    
    void render(const vector<Block>&) {
        logger.error("%s\n", "Fatal error. This function should not be used. Aborting");
        exit(1);
        return;
    }
    
    void render(const vector<Mtf_profile_sample>& samples, Bayer::bayer_t bayer = Bayer::NONE, 
        vector<Ellipse_detector>* ellipses = NULL, cv::Rect* dimension_correction = NULL);

  private:
    VectorXd rpfit(Ratpoly_fit& cf, bool scale=true, bool refine=true);
    void exposure_checks(const Point2d& dims, double& white_clip, double& black_clip, double& overexposure);
  
    Point2d& zero;
    Point2d& transverse;
    Point2d& longitudinal;
    std::string wdir;
    std::string prname;
    const cv::Mat& img;
    Distance_scale& distance_scale;
    const vector<std::pair<Point2d, Point2d>>& sliding_edges;
    int initial_rows;
    
    double psf;
    
    Camera_draw draw;
};

#endif
