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
#ifndef GRADIENT_H
#define GRADIENT_H

#include <cmath>
#include "include/common_types.h"

class Gradient {
public:
    Gradient(const cv::Mat& in_img);
    virtual ~Gradient(void);

    inline const cv::Mat& grad_x(void) const {
        return _gradient_x;
    }
    
    inline void write_grad_x(int x, int y, float val) {
        _gradient_x.at<float>(y, x) = val;
    }

    inline const cv::Mat& grad_y(void) const {
        return _gradient_y;
    }
    
    inline void write_grad_y(int x, int y, float val) {
        _gradient_y.at<float>(y, x) = val;
    }
    
    inline float grad_x(int x, int y) const {
        return _gradient_x.at<float>(y, x);
    }

    inline float grad_y(int x, int y) const {
        return _gradient_y.at<float>(y, x);
    }

    inline float grad_magnitude(int x, int y) const {
        size_t i = _gradient_x.cols * y + x;
        return SQR(*((float*)_gradient_x.data + i)) + SQR(*((float*)_gradient_y.data + i));
    }

    inline int width(void) const {
        return _width;
    }

    inline int height(void) const {
        return _height;
    }
    
    void release(void) {
        _gradient_x.release();
        _gradient_y.release();
    }
    
private:

    void _compute_gradients(const cv::Mat& smoothed_im);

protected:
    int _width;
    int _height;

    cv::Mat     _gradient_x;
    cv::Mat     _gradient_y;
};

#endif // GRADIENT_H


