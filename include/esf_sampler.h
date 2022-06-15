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
#ifndef ESF_SAMPLER_H
#define ESF_SAMPLER_H

#include "include/logger.h"
#include "include/common_types.h"
#include "include/ordered_point.h"
#include <opencv2/imgproc/imgproc.hpp>

#include "include/bayer.h"
#include "include/edge_model.h"
#include <map>
using std::map;

class Esf_sampler {
  public:
    
    // enumeration and string array (see .cc file) should be kept synchronized
    typedef enum {
        LINE=0,
        QUAD,
        PIECEWISE_QUAD,
        DEFERRED
    } esf_sampler_t;
    
    const static std::array<string, 4> esf_sampler_names;
    
    static esf_sampler_t from_string(const string& s) {
        for (size_t i=0; i < esf_sampler_names.size(); i++) {
            if (s.compare(esf_sampler_names[i]) == 0) {
                return static_cast<esf_sampler_t>(i);
            }
        }
        return LINE; // fallback
    }
  
    Esf_sampler(double max_dot, Bayer::cfa_mask_t cfa_mask=Bayer::ALL, double max_edge_length=1e6, double border_width=0) 
    : max_dot(max_dot), default_cfa_mask(cfa_mask), max_edge_length(max_edge_length), border_width(border_width) {
        
    }
    
    virtual void sample(Edge_model& edge_model, vector<Ordered_point>& local_ordered, 
        const map<int, scanline>& scanset, double& edge_length,
        const cv::Mat& geom_img, const cv::Mat& sampling_img, 
        Bayer::cfa_mask_t cfa_mask = Bayer::DEFAULT) = 0;
        
  protected:
    double max_dot;
    Bayer::cfa_mask_t default_cfa_mask;
    double max_edge_length;
    double border_width;
};

#endif
