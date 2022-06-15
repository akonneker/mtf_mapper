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
#ifndef AIRY_SAMPLER_H
#define AIRY_SAMPLER_H

#include "include/common_types.h"

class Airy_sampler {
  public:
    Airy_sampler(double lambda, double pitch, double N) 
      : lambda(lambda), pitch(pitch), N(N), diam(std::max(85.0, 45.0*pitch/(lambda*N))) {
        
        diam_prob = 1 - j0(diam*M_PI)*j0(diam*M_PI) - j1(diam*M_PI)*j1(diam*M_PI);
    }
    
    double rairy2d(double& sx, double& sy, normal_sampler& sampler) {
        double px;
        double py;
        
        do {
            sampler.runif2d(px, py, 0.5, 0.5);
            px += 0.5; // undo the scaling of runif2d
            py += 0.5;
            py *= 2*M_PI;
        } while (px > diam_prob);
        
        const double tol = 1e-11;
        double gr = (sqrt(5.0)-1.0)*0.5;
        double a = 0;
        double b = diam; // larger?
        double c = b - gr*(b - a);
        double d = a + gr*(b - a);
        while (fabs(c - d) > tol) {
            double fc = jinc_p(c, px);
            double fd = jinc_p(d, px);
            if (fc < fd) {
                b = d;
                d = c;
                c = b - gr*(b - a);
            } else {
                a = c;
                c = d;
                d = a + gr*(b - a);
            }
        }
        double rad = (0.5*(a+b)) * (lambda/pitch)*N;
        sx = rad*cos(py);
        sy = rad*sin(py);
        
        return 4*jinc(rad*pitch/(lambda*N))*jinc(rad*pitch/(lambda*N));
    }
    
    inline static double jinc(double v) {
        if (fabs(v) == 0) {
            return 0.5;
        }
        return j1(v*M_PI)/(v*M_PI);
    }
    
    inline double jinc_p(double v, double target) {
        return fabs((1 - j0(v*M_PI)*j0(v*M_PI) - j1(v*M_PI)*j1(v*M_PI)) - target);
    }
    
    double diam_prob;
    
    double lambda;
    double pitch;
    double N;
    double diam;
    
    size_t nsamples;
    vector<double> x;
    vector<double> y;
};

#endif

