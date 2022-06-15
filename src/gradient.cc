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
#include "include/gradient.h"
#include "include/common_types.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <opencv2/imgproc/imgproc.hpp>

//------------------------------------------------------------------------------
Gradient::Gradient(const cv::Mat& in_img)
 : _width(in_img.cols), _height(in_img.rows)
{
    
    cv::Mat in_float;
    double min_val = 0;
    double max_val = 0;
    minMaxLoc(in_img, &min_val, &max_val);
    in_img.convertTo(in_float, CV_32FC1, 1.0/max_val);
    
    cv::Mat smoothed;
    cv::GaussianBlur(in_float, smoothed, cv::Size(5,5), 1.2, 1.2);
    
    _compute_gradients(smoothed);
    smoothed.release();
}

//------------------------------------------------------------------------------
Gradient::~Gradient(void) {
}


//------------------------------------------------------------------------------
void Gradient::_compute_gradients(const cv::Mat& smoothed_im) {

    _gradient_x = cv::Mat(smoothed_im.rows, smoothed_im.cols, CV_32FC1);
    _gradient_y = cv::Mat(smoothed_im.rows, smoothed_im.cols, CV_32FC1);

    float* gxp = (float*)_gradient_x.data + 1;
    float* smp = (float*)smoothed_im.data + 1;
    for(int i=0; i < smoothed_im.rows*smoothed_im.cols - 2; i++) { // skip first and last pixels
        *gxp = *(smp+1) - *(smp-1);
        gxp++;
        smp++;
    }
    for(int r=0; r < smoothed_im.rows; r++){
        _gradient_x.ptr<float>(r)[0] = 0;
        _gradient_x.ptr<float>(r)[smoothed_im.cols-1] = 0;
    }
    float* gyp = (float*)_gradient_y.data + _gradient_y.cols;
    smp = (float*)smoothed_im.data + smoothed_im.cols;
    for(int i=0; i < smoothed_im.rows*smoothed_im.cols - 2*smoothed_im.cols; i++) { // skip first and last row
        *gyp = *(smp+smoothed_im.cols) - *(smp-smoothed_im.cols);
        gyp++;
        smp++;
    }
    for (int c=0; c < smoothed_im.cols; c++) {
        _gradient_y.ptr<float>(0)[c] = 0;
        _gradient_y.ptr<float>(smoothed_im.rows-1)[c] = 0;
    }
}

