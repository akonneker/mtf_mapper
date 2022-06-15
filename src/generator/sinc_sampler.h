/*
Copyright 2021 Frans van den Bergh. All rights reserved.

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
#ifndef SINC_SAMPLER_H
#define SINC_SAMPLER_H

#include "include/common_types.h"
#include "sine_integral.h"

class Sinc_sampler {
  public:
    Sinc_sampler(double lambda, double pitch, double N, double N2) 
      : lambda(lambda), pitch(pitch), N(N), N2(N2) {
        
        diam = sine_integral.size() * sine_integral_step;
        diam_prob = 1; // TODO: can we improve on this
    }
    
    double rairy2d(double& sx, double& sy, normal_sampler& sampler) {
        double px, orig_px;
        double py, orig_py;
        
        do {
            sampler.runif2d(orig_px, orig_py, 1, 1);
            px = fabs(orig_px);
            py = fabs(orig_py);
        } while (px > diam_prob || py > diam_prob);
        
        // do the lookups
        size_t x_idx = lower_bound(sine_integral.begin(), sine_integral.end(), px) - sine_integral.begin();
        size_t y_idx = lower_bound(sine_integral.begin(), sine_integral.end(), py) - sine_integral.begin();
        
        // slope of sine_integral table will give us the estimate of sinc^2
        size_t next_x_idx = (x_idx + 1) % sine_integral.size();
        size_t next_y_idx = (y_idx + 1) % sine_integral.size();
        double sinc2_x = (sine_integral[next_x_idx] - sine_integral[x_idx]) / sine_integral_step;
        double sinc2_y = (sine_integral[next_y_idx] - sine_integral[y_idx]) / sine_integral_step;
        
        double pscale = (1.0 / M_PI) * diam / double(sine_integral.size());
        
        // note that this is the 1st quadrant value only
        sx = pscale * x_idx * (lambda/pitch)*N;
        sy = pscale * y_idx * (lambda/pitch)*N2;
        
        // move point to original quadrant
        sx *= orig_px < 0 ? -1 : 1;
        sy *= orig_py < 0 ? -1 : 1;
        
        return sine_integral_max * sine_integral_max * sinc2_x * sinc2_y;
    }
    
    inline static double sinc(double v) {
        if (fabs(v) == 0) {
            return 1.0;
        }
        return sin(v*M_PI)/(v*M_PI);
    }
    
    double weight(double px, double py) const {
        double sx = sinc(px*pitch/(lambda*N));
        double sy = sinc(py*pitch/(lambda*N2));
        return sx*sx*sy*sy;
    }
    
    double diam_prob;
    
    double lambda;
    double pitch;
    double N;
    double N2;
    double diam;
    
    size_t nsamples;
    vector<double> x;
    vector<double> y;
};

#endif

