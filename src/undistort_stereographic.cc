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

#include "include/undistort_stereographic.h"
    
cv::Point2d Undistort_stereographic::slow_transform_point(double col, double row) {
    double px = (col + offset.x - centre.x)*pitch;
    double py = (row + offset.y - centre.y)*pitch;

    double rd = sqrt((px)*(px) + (py)*(py)); // radial distance in mm

    if (rd == 0) {
        return cv::Point2d(centre.x - offset.x, centre.y - offset.y);
    }

    double theta = atan(rd/f);
    double ru = tan(0.5*theta)*f*2.0;
    
    px = px*ru/rd / pitch + centre.x - offset.x;
    py = py*ru/rd / pitch + centre.y - offset.y;
    
    return cv::Point2d(px, py);
}

cv::Point2d Undistort_stereographic::inverse_transform_point(double col, double row) {
    double px = (col + offset.x - centre.x)*pitch;
    double py = (row + offset.y - centre.y)*pitch;
    
    double ru = sqrt((px)*(px) + (py)*(py)); // radial distance in mm

    if (ru == 0) {
        return cv::Point2d(centre.x - offset.x, centre.y - offset.y);
    }

    double theta = atan(0.5*ru/f);
    double rd = tan(2*theta)*f;
    
    px = px*rd/ru / pitch + centre.x - offset.x;
    py = py*rd/ru / pitch + centre.y - offset.y;
    
    return cv::Point2d(px, py);
}

// note: second parameter is the raw Bayer image, which must also be padded out
cv::Mat Undistort_stereographic::unmap(const cv::Mat& in_src, cv::Mat& rawimg) {
    
    int pad_left = 0;
    int pad_top = 0;
    estimate_padding(in_src, pad_left, pad_top);

    // TODO: perform some sanity checking ...
    pad_left = std::min(8000, pad_left);
    pad_top = std::min(8000, pad_top);

    return unmap_base(in_src, rawimg, pad_left, pad_top);
}

    
