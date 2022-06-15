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

#include "include/ratpoly_fit.h"

double Ratpoly_fit::evaluate(VectorXd& v) {
    double err = 0;
    for (size_t i=0; i < data.size(); i++) {
        double w = data[i].weight * data[i].yweight;
        double z = rpeval(v, scale(data[i].x));
        double e = data[i].y*ysf - z;
        err += e*e*w;
    }
    evaluation_count++;
    return err*0.5;
}

VectorXd Ratpoly_fit::gauss_newton_direction(VectorXd& v, VectorXd& deriv, double& fsse) {
    MatrixXd J(data.size(), v.rows());
    J.setZero();
    fsse = 0; 
    
    VectorXd r(data.size());
    for (size_t m=0; m < data.size(); m++) {
        double w = data[m].weight * data[m].yweight;
        double fx = 0;
        
        J.row(m) = rp_deriv(v, scale(data[m].x), fx); 
        double e = fx - data[m].y*ysf;
        r[m] = e*w;
        fsse += e*e*w;
    }
    
    
    VectorXd direction = J.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(-r);
    
    fsse *= 0.5;
    deriv = J.transpose() * r;
    evaluation_count++;
    return direction;
}

VectorXd Ratpoly_fit::gauss_newton_armijo(VectorXd& v) {
    const double tau = 0.5;
    const double c = 1e-4;
    double fx = 0;
    
    VectorXd grad;
    VectorXd next;
    VectorXd pk;
    for (int k=0; k < 50; k++) {
        
        double alpha = 1.0;
        pk = gauss_newton_direction(v, grad, fx);
        
        double target = fx + c*alpha*pk.dot(grad);
        
        int max_steps = 30;
        next = v + alpha*pk;
        while (evaluate(next) > target && --max_steps > 0) { // iteratively step close until we have a sufficient decrease (Armijo condition)
            target = fx + c*alpha*pk.dot(grad);
            alpha *= tau;
            next = v + alpha*pk;
        }
        
        double stepsize = pk.array().abs().maxCoeff()*fabs(alpha);
        if (stepsize < 5e-8) {
            break;
        }
        v = next;
    }
    return v;
}

double Ratpoly_fit::peak(const VectorXd& v) {
    double xmin=1e50;
    double xmax=-1e50;
    for (size_t i=0; i < data.size(); i++) {
        xmin = std::min(data[i].x, xmin);
        xmax = std::max(data[i].x, xmax);
    }
    // bracket the maximum
    double peak_z = 0;
    double peak_x = (xmin + xmax)*0.5;
    double step = (xmax - xmin)/20.0;
    for (double x=xmin; x <= xmax; x += step) {
        double z = rpeval(v, scale(x));
        if (z > peak_z) {
            peak_x = x;
            peak_z = z;
        }
    }
    
    // golden section search
    const double phi = 0.61803398874989;
    double lower = peak_x - 2*step;
    double upper = peak_x + 2*step;
    double c = upper - phi*(upper - lower);
    double d = lower + phi*(upper - lower);
    const double tol = 1e-10;
    while ((upper - lower) > tol) {
        double fc = rpeval(v, scale(c));
        double fd = rpeval(v, scale(d));
        if (fc > fd) {
            upper = d;
            d = c;
            c = upper - phi*(upper - lower);
        } else {
            lower = c;
            c = d;
            d = lower + phi*(upper - lower);
        }
    }
    return 0.5*(upper + lower);
}

bool Ratpoly_fit::has_poles(const VectorXd& v) {
    double xmin=1e50;
    double xmax=-1e50;
    for (size_t i=0; i < data.size(); i++) {
        xmin = std::min(data[i].x, xmin);
        xmax = std::max(data[i].x, xmax);
    }
    
    // ensure the bounds are slightly wider than the actual data
    double span=xmax - xmin;
    xmin -= pscale*span;
    xmax += pscale*span;
    
    // now compute roots of bottom polynomial
    switch(order_m) {
    case 0:
        return false; // cannot have poles
    case 1:
        {
            double pole = unscale(-1.0 / v[order_n+1]);
            
            
            if (!silent && pole >= xmin && pole <= xmax) {
                logger.debug("pole at %lf on [%lf, %lf]\n", pole, xmin, xmax);
            }
            
            
            return pole >= xmin && pole <= xmax;
        }
    case 2:
        {
            double a = v[order_n+2];
            double b = v[order_n+1];
            double c = base_value - v[order_n+2]; // since this is a chebyshev polynomial a(x^2 - 1) -> c' = c - a
            double sb = b < 0 ? -1 : 1;
            double q = -0.5*(b + sb*sqrt(b*b - 4*a*c));
            double pole1 = unscale(q/a);
            double pole2 = unscale(c/q);
            
            if (!silent && ((pole1 >= xmin && pole1 <= xmax) || (pole2 >= xmin && pole2 <= xmax))) {
                logger.debug("pole at %lf or %lf on [%lf, %lf]\n", pole1, pole2, xmin, xmax);
            }
                          
            return (pole1 >= xmin && pole1 <= xmax) ||
                   (pole2 >= xmin && pole2 <= xmax);
        }
    default:
        // TODO: see NR chapter 5.6 for cubic roots
        logger.error("Warning: no implementation to compute roots of order-%d polynomial\n",
            order_m
        );
        return false;
    };
}

