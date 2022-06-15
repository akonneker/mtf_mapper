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

#ifndef UNDISTORT_H
#define UNDISTORT_H

#include "common_types.h"
#include "include/logger.h"
#include <opencv2/imgproc/imgproc.hpp>

#include <map>
using std::map;

#include <cmath>

class Undistort { 
  public:
    Undistort(const cv::Rect& r) : centre(r.width / 2, r.height / 2), offset(r.x, r.y), radmap(0), max_val(1, 1), last_padding(0, 0) {};
    
    cv::Point2i transform_pixel(int col, int row) {
        cv::Point2d tp = transform_point(double(col), double(row));
        return cv::Point2i(lrint(tp.x), lrint(tp.y));
    }
    
    cv::Point2d transform_point(double px, double py); // this function interpolates radmap
    
    virtual cv::Point2d slow_transform_point(double px, double py) = 0; // the "real" forward transformation, which could be slow
    virtual cv::Point2d inverse_transform_point(double px, double py) = 0;
    
    cv::Point2d transform_point(const cv::Point2d& p) {
        return transform_point(p.x, p.y);
    }
    
    virtual cv::Mat unmap(const cv::Mat& src, cv::Mat& rawimg) = 0;
    
    bool rectilinear_equivalent(void) const {
        return rectilinear;
    }
    
    void set_rectilinear_equivalent(bool b) {
        rectilinear = b;
    }
    
    void set_max_val(const Point2d& maxv) { 
        max_val = maxv;
    }

    void set_allow_crop(bool crop) {
        allow_crop = crop;
    }

    void apply_padding(vector<cv::Mat>& images);
    
    cv::Point2d centre;
    cv::Point2d offset;
    vector<double> radmap; // a vector of the transformed radius sampled at uniform spacing in the untransformed radius
    Point2d max_val;
    
    bool rectilinear = false;
    bool allow_crop = true;
    
    
  protected:  
    void build_radmap(void); // called from derived class constructor once parameters are known
    cv::Mat unmap_base(const cv::Mat& src, cv::Mat& rawimg, int pad_left, int pad_top);
    void estimate_padding(const cv::Mat& src, int& pad_left, int& pad_top);

    cv::Point2i last_padding;
};
    
#endif
    
