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
#ifndef AFFT_H
#define AFFT_H

#include <assert.h>

#include <vector>
using std::vector;
using std::pair;
using std::make_pair;

template< int n >
class AFFT {
  public:
    AFFT(void) {
    
        int j = 1;
        // digit (bit) reverse
        for (int i=1; i < n; i++) {
            if (i < j) {
                bitrev.push_back(make_pair(j,i));
            }
            int k = n / 2;
            while (k >= 1 && k < j) {
                j -= k;
                k /= 2;
            }
            j += k;
        }
        
        
        power = 0;
        unsigned int t = n/2;
        while (t > 0) {
            power++;
            t >>= 1;
        }
        
        cstab = vector< vector<double> >(power+1, vector<double>(n/2, 0)); // choose fixed upper limit to help compiler
        int n2 = 1;
        for (unsigned int k=2; k <= power; k++) {
            int n4 = n2;
            n2 = 2*n4;
            int n1 = 2*n2;
            double e = 2*M_PI / double(n1);
            int idx = 0;
            for (int i=1; i <= n; i += n1) {
                double a = e;
                for (int j=1; j <= (n4 - 1); j++) { 
                    cstab[k][idx++] = cos(a);
                    cstab[k][idx++] = sin(a);
                    a += e;
                }
            }
        }
    }
    
    // NB: Input is n real samples
    //     Output is [Re(0), Re(1), ..., Re(N/2), Im(N/2-1), ..., Im(1)]
    //     So DC is x[0], and complex frequency k is (x[k], x[N-k])
    void realfft(double* x) {
        x--; // simulate 1-based arrays
        
        // TODO: can we combine the first pass with the bit reversal?
        
        for (size_t i=0; i < bitrev.size(); i++) {
            std::swap(x[bitrev[i].first], x[bitrev[i].second]);
        }
        
        // length 2 butterflies have special twiddle factors, do them first
        double* xp = x + 1;
        double* xp_sent = x + n;
        for (; xp <= xp_sent; xp += 2) {
            double xt = *xp;
            *(xp)   = xt + *(xp+1);
            *(xp+1) = xt - *(xp+1);
        }
        
        // other stages
        int n2 = 1;
        for (unsigned int k=2; k <= power; k++) {
            int n4 = n2;
            n2 = n4 << 1;
            int n1 = n2 << 1;
            double* cs_ptr = cstab[k].data();
            
            double* xp = x + 1; // start at x[1]
            double* xp_sent = xp + n; // stage sentinel
            
            for (;xp < xp_sent;) { 
                double xt = *xp;
                *(xp)    = xt + *(xp+n2);
                *(xp+n2) = xt - *(xp+n2);
                *(xp+n4+n2) = -*(xp+n4+n2);
                
                for (int j=1; j <= (n4 - 1); j++) { 
                    double* i1 = xp + j;      
                    double* i2 = xp - j + n2;
                    double* i3 = i1 + n2;
                    double* i4 = i2 + n2;
                    
                    double t1 = *(i1+n2)*(*(cs_ptr))  +  *(i2+n2)*(*(cs_ptr+1));
                    double t2 = *(i1+n2)*(*(cs_ptr+1)) - *(i2+n2)*(*(cs_ptr));
                    *i4 =  *i2 - t2;
                    *i3 = -*i2 - t2;
                    
                    *i2 =  *i1 - t1;
                    *i1 += t1;
                    cs_ptr += 2;
                }
                xp += n1;
            }
        }
    }
    
    vector< pair<int,int> > bitrev;
    vector< vector<double> > cstab;
    unsigned int power;
};

#endif
