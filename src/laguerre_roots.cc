/*
Copyright 2020 Frans van den Bergh. All rights reserved.

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

#include "include/laguerre_roots.h"
#include <limits>
#include <algorithm>

bool laguerre(const vector<cplex>& a, cplex& x, int& its) {
    const int MR=8;
    const int MT=10;
    const int MAXIT=MT*MR;
    
    const double EPS=std::numeric_limits<double>::epsilon();
    static const double frac[MR+1] = {0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0};
    
    int m = a.size() - 1;
    for (int iter=1; iter <= MAXIT; iter++) {
        its = iter;
        cplex b = a[m];
        double err = std::abs(b);
        cplex d(0,0);
        cplex f(0,0);
        double abx = std::abs(x);
        for (int j=m-1; j >= 0; j--) {
            f = x*f + d;
            d = x*d + b;
            b = x*b + a[j];
            err = std::abs(b) + abx*err;
        }
        err *= EPS;
        if (std::abs(b) <= err) return true; // on the root
        cplex g = d/b;
        cplex g2 = g*g;
        cplex h = g2 - 2.0*f/b;
        cplex sq = std::sqrt(double(m-1)*(double(m)*h-g2));
        cplex gp = g + sq;
        cplex gm = g - sq;
        double abp = std::abs(gp);
        double abm = std::abs(gm);
        gp = (abp < abm) ? gm : gp;
        cplex dx = std::max(abp, abm) > 0.0 ? double(m)/gp : std::polar(1+abx, double(iter));
        cplex x1 = x - dx;
        if (x == x1) return true;
        if (iter % MT != 0) {
            x = x1;
        } else {
            x -= frac[iter/MT]*dx;
        }
    }
    // no convergence ...
    return false;
}

void lroots(const vector<cplex>& a, vector<cplex>& roots, bool polish) {
    const double EPS = 1e-14;
    int its;
    
    int m = a.size() - 1;
    vector<cplex> ad(a);
    for (int j=m-1; j >= 0; j--) {
        cplex x(0,0);
        vector<cplex> ad_v(j+2);
        for (int jj=0; jj < j+2; jj++) {
            ad_v[jj] = ad[jj];
        }
        laguerre(ad_v, x, its);
        if (fabs(x.imag()) <= 2.0*EPS*fabs(x.real())) {
            x = cplex(x.real(), 0.0);
        }
        roots[j] = x;
        cplex b = ad[j+1];
        for (int jj=j; jj >= 0; jj--) { // deflate the poly
            cplex c = ad[jj];
            ad[jj] = b;
            b = x*b + c;
        }
    }
    
    if (polish) {
        for (int j=0; j < m; j++) {
            laguerre(a, roots[j], its);
        }
    }
}

