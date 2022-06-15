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

#ifndef RATPOLY_FIT_H
#define RATPOLY_FIT_H

#include "include/logger.h"
#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;

#include "include/laguerre_roots.h"

class Sample {
  public:
    Sample(double x, double y, double weight, double yweight=1) 
    :x(x), y(y), weight(weight), yweight(yweight) {}
    
    bool operator< (const Sample& b) const {
        return x < b.x;
    }
    
    double x;
    double y;
    double weight;
    double yweight;
};

class Ratpoly_fit  {
  public:
    
    Ratpoly_fit(const vector<Sample>& data, int order_n, int order_m, bool silent=false)
    : data(data), order_n(order_n), order_m(order_m), base_value(1.0), 
      ysf(1), pscale(0.1), silent(silent), evaluation_count(0) {}
    
    
    virtual int dimension(void) {
        return (order_n+1 + order_m);
    }
    
    inline double rpeval(const VectorXd& v, double x) {
        double top_val = v[0];
        double p = x;
        for (int n=1; n <= order_n; n++) {
            top_val += cheb(n, p)*v[n];
        }
        double bot_val = base_value;
        p = x;
        for (int m=0; m < order_m; m++) {
            bot_val += cheb(m+1, p)*v[m+order_n+1];
        }
        
        return top_val / bot_val;
    }
    
    inline VectorXd rp_deriv(const VectorXd& v, double x, double& f) {
        // if x falls on a pole, we are in trouble
        // and should probably just return the zero vector?
        
        // TODO: we can probably combine this function with rp_deriv_eval ??
    
        VectorXd dp = VectorXd::Zero(v.rows());
        VectorXd dq = VectorXd::Zero(v.rows());
        
        double top_val = v[0];
        dp[0] = 1;
        for (int n=1; n <= order_n; n++) {
            dp[n] = cheb(n, x);
            top_val += cheb(n, x)*v[n];
        }
        double bot_val = base_value;
        for (int m=0; m < order_m; m++) {
            dq[m+order_n+1] = cheb(m+1, x);
            bot_val += cheb(m+1, x)*v[m+order_n+1];
        }
        
        double den = bot_val*bot_val;
        if (den < 1e-12) {
            return VectorXd::Zero(v.rows());
        }
        
        f = top_val / bot_val;
        den = 1.0/den;
        
        return (dp*bot_val - top_val*dq) * den;
    }
    
    const vector<Sample>& get_data(void) const {
        return data;
    }
    
    inline double cheb(int order, double x) {
        double y = 0;
        switch (order) {
            case 0: y = 1; break;
            case 1: y = x; break;
            case 2: y = 2*x*x - 1; break;
            case 3: y = 4*x*x*x - 3*x; break;
            case 4: y = 8*x*x*x*x - 8*x*x + 1; break;
            case 5: y = 16*x*x*x*x*x - 20*x*x*x + 5*x; break;
            case 6: y = 32*x*x*x*x*x*x - 48*x*x*x*x + 18*x*x - 1; break;
            default: fprintf(stderr, "Unsupported Cheb poly order requested!\n"); break;
        }
        return y;
    }
    
    double evaluate(VectorXd& v);
    VectorXd evaluate_derivative(VectorXd& v);
    VectorXd gauss_newton_direction(VectorXd& v, VectorXd& deriv, double& fsse);
    VectorXd gauss_newton_armijo(VectorXd& v);
    double peak(const VectorXd& v);
    bool has_poles(const VectorXd& v);
    
    const vector<Sample>& data;
    int order_n;
    int order_m;
    double base_value;
    
    double ysf;
    double pscale;
    bool silent;
    
    double xs_min;
    double xs_scale;
    
    inline double scale(double x) const {
        return (x - xs_min) * xs_scale;
    }
    
    inline double unscale(double x) const {
        return x / xs_scale + xs_min;
    }
    
    unsigned long long evaluations(void) {
        return evaluation_count;
    }
                                
  protected:
    unsigned long long evaluation_count;
};

#endif
