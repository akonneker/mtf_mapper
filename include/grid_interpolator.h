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
#ifndef GRID_INTERPOLATOR_H
#define GRID_INTERPOLATOR_H
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "include/block.h"

typedef enum {MERIDIONAL, SAGITTAL, NEITHER} Edge_type;

class Grid_functor {
  public:
    virtual double value(const Block& /*block*/, size_t /*edge*/) const {
        return 0;
    }
    
    virtual bool in_range(double /*val*/) const {
        return true; 
    }
    
    virtual double clamp(double value, double /*upper*/, double /*lower*/) const {
        return value;
    }
    
    virtual double nodata(void) const {
        return 0;
    }

    virtual double smoothing_factor(void) const {
        return 1e-3;
    }

    virtual int pruning_threshold(void) const {
        return 1;
    }
};

void interpolate_grid(const Grid_functor& ftor, Edge_type target_edge_type, cv::Mat& grid_coarse, cv::Mat& grid_fine, 
    cv::Size img_dims, const vector<Block>& blocks, double upper, double smoothing_factor = 1e-3, int pruning_threshold = 1);

#endif
