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


#include "include/undistort_rectilinear.h"
#include "include/laguerre_roots.h"

Undistort_rectilinear::Undistort_rectilinear(const cv::Rect& r, const vector<double>& in_coeffs) : Undistort(r) {
    coeffs = vector<double>(2, 0);
    for (size_t i=0; i < 2; i++) coeffs[i] = in_coeffs[i];
    
    radius_norm = sqrt(centre.x*centre.x + centre.y*centre.y);
}

static double laguerre_smallest_positive_root(double ru, double k1, double k2) {
    vector<cplex> a = {ru, -1.0, ru*k1, 0.0, ru*k2};
    vector<cplex> roots(4);
    lroots(a, roots);
    
    double minroot = 1e30;
    for (size_t i=0; i < roots.size(); i++) {
        if (roots[i].imag() == 0 && roots[i].real() >= 0 && roots[i].real() < minroot) {
            minroot = roots[i].real();
        }
    }
    
    // try to refine the root
    cplex root(minroot, 0.0);
    int its;
    laguerre(a, root, its);
    
    return root.real();
}

static bool lagsolve(double ru, double k1, double k2, double& root) {
    // we are looking for the roots of
    // P(r) = ru*k2*root^4 + ru*k1*root^2 - root + ru = 0
    // we would prefer the smallest positive root
    
    if (fabs(k1) < 1e-8 && fabs(k2) < 1e-8) {
        root = ru;
    } else {
        if (fabs(k2) < 1e-8) { // we only have a quadratic
            // a == 1
            double b = -1.0/(k1*ru);
            double c = 1.0/k1;
            double q = -0.5*(b + std::copysign(sqrt(b*b - 4*c), b));
            // force negative roots to become very large
            double r1 = q;
            double r2 = c/q;
            if (r1 < 0) {
                root = r2;
            } else {
                if (r2 < 0) {
                    root = r1;
                } else {
                    root = std::min(r1, r2);
                }
            }
        } else {
            root = laguerre_smallest_positive_root(ru, k1, k2);    
        }
    }
    return true;
}

cv::Point2d Undistort_rectilinear::slow_transform_point(double col, double row) {
    double px = (col + offset.x - centre.x);
    double py = (row + offset.y - centre.y);
    
    double ru = sqrt((px)*(px) + (py)*(py)) / radius_norm; 
    
    if (ru == 0) {
        px = col;
        py = row;
    } else {
        double rd=0;
        bool s = lagsolve(ru, coeffs[0], coeffs[1], rd); // not the fastest method, but it should be robust
        
        if (s) {
            px = px*rd/ru + centre.x - offset.x;
            py = py*rd/ru + centre.y - offset.y;
        } else {
            px = col;
            py = row;
        }
    }
    
    return cv::Point2d(px, py);
}    

cv::Point2d Undistort_rectilinear::inverse_transform_point(double col, double row) { // explicit inverse model
    double px = (col + offset.x - centre.x);
    double py = (row + offset.y - centre.y);

    double rd = sqrt((px)*(px) + (py)*(py)) / radius_norm; 

    if (rd == 0) {
        return cv::Point2d(centre.x - offset.x, centre.y - offset.y);
    }

    double r2 = rd*rd;
    double ru = 1 + (coeffs[0] + coeffs[1]*r2)*r2;
    
    px = px/ru + centre.x - offset.x;
    py = py/ru + centre.y - offset.y;
    
    return cv::Point2d(px, py);
}

// note: second parameter is the raw Bayer image, which must also be padded out
cv::Mat Undistort_rectilinear::unmap(const cv::Mat& in_src, cv::Mat& rawimg) {
    
    const double buffer = 0.025;
    Point2d pi = inverse_transform_point(centre.x - (max_val.x+buffer*std::max(in_src.cols, in_src.rows)), centre.y);
    int pad_left = pi.x < 0 ? ceil(-pi.x) : 0;
    pi = inverse_transform_point(centre.x, centre.y - (max_val.y+buffer*std::max(in_src.cols, in_src.rows)));
    int pad_top = pi.y < 0 ? ceil(-pi.y) : 0;
    
    return unmap_base(in_src, rawimg, pad_left, pad_top);
}
    
    
