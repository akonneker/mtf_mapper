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
#ifndef QMC_NORMAL_SAMPLING_H
#define QMC_NORMAL_SAMPLING_H

#include <cmath>
#include <stdint.h>

class normal_sampler {
  public:
  
    normal_sampler(uint64_t in_seed=20) : seed(in_seed) {
    }
  
    double rnorm(double sigma) {
        double r = Moro_norm_inv(corput_base_b(2, seed));
        seed++;
        return r * sigma;    
    }
    
    void rnorm2d(double &x, double& y, double sigma=1.0, double sigma_minor=1.0, double theta=0) {
        double lx = Moro_norm_inv(corput_base_b(2, seed)) * sigma;
        double ly = Moro_norm_inv(corput_base_b(3, seed)) * sigma_minor;
        seed++;
        x = lx*cos(theta) - ly*sin(theta);
        y = lx*sin(theta) + ly*cos(theta);
    }
    
    void runif2d(double &x, double& y, double sigma=1.0, double sigma_minor=1.0, double theta=0) {
        double lx = (corput_base_b(2, seed) - 0.5) * 2 * sigma;
        double ly = (corput_base_b(3, seed) - 0.5) * 2 * sigma_minor;
        seed++;
        x = lx*cos(theta) - ly*sin(theta);
        y = lx*sin(theta) + ly*cos(theta);
    }
    
    static double corput_base_b(uint64_t b, uint64_t N) {
        double c = 0.0;
        double ib = 1.0/double(b);
        uint64_t n1 = N;
        uint64_t n2;
        
        while (n1 > 0) {
            n2 = n1 / b;
            c = c + ib * (n1 % b);
            ib = ib / b;
            n1 = n2;
        }
        
        return c;
    }
  
    static inline double Moro_norm_inv(double u) { // assumes u in (0,1)
        
        double c[] = {
            0.337475482272615, 0.976169019091719, 0.160797971491821, 2.76438810333863E-02, 
            3.8405729373609E-03, 3.951896511919E-04, 3.21767881768E-05, 2.888167364E-07,
            3.960315187E-07
        };
                                         
        
        double x = u - 0.5;
        
        double r = 0;
        if (fabs(x) < 0.42) {
            double a[] = {2.50662823884, -18.61500062529, 41.39119773534, -25.44106049637};
            double b[] = {-8.4735109309, 23.08336743743, -21.06224101826, 3.13082909833};
          
            r = x*x;
            r = x * 
              (((a[3] * r + a[2]) * r + a[1]) * r + a[0]) / 
              ((((b[3] * r + b[2]) * r + b[1]) * r + b[0]) * r + 1);

        } else {
			if (x > 0) {
				r = log(-log(1-u));
			} else {
				r = log(-log(u));
			}
			if (!(r <= DBL_MAX && r >= -DBL_MAX)) { // catch infinities returned by log() on MSVC
				r = 0;
			}
            r = c[0] + r * (c[1] + r * (c[2] + r * (c[3] + r * (c[4] + r * (c[5] + r * (c[6] + r * (c[7] + r * c[8])))))));
            if (x <= 0) {
                r = -r;
            }
        }
        return r;
    }
    
    uint64_t seed;
};


#endif 
