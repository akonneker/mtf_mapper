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
#ifndef ESF_MODEL_H
#define ESF_MODEL_H

#include "include/ordered_point.h"
#include "include/sampling_rate.h"
#include "include/snr.h"
#include <vector>
#include <string>
using std::vector;
using std::string;

class Esf_model {
  public:
    Esf_model(double alpha=13.7) 
    : alpha(alpha), w(NYQUIST_FREQ*4, 0.0)  {}
    
    virtual ~Esf_model(void) {}
    
    virtual int build_esf(vector< Ordered_point  >& ordered, double* sampled, 
        const int fft_size, double max_distance_from_edge, vector<double>& esf, 
        Snr& snr, bool allow_peak_shift=false) = 0;
        
    void moving_average_smoother(vector<double>& smoothed, double* sampled, int fft_size, 
        int fft_left, int fft_right, int left_trans, int right_trans, int width=16);
        
    int estimate_esf_clipping(vector< Ordered_point  >& ordered, double* sampled, 
        const int fft_size, bool allow_peak_shift, int effective_maxdot, vector<double>& mean,
        vector<double>& weights, int& fft_left, int& fft_right, int& twidth, Snr& snr);
        
    virtual void set_alpha(double a) {
        alpha = a;
        compute_mtf_corrections();
    }
    
    virtual void compute_mtf_corrections(void);
    
    const vector<double>& get_correction(void) const {
        return w;
    }
    
    void set_monotonic_filter(bool b) {
        apply_monotonic_filter = b;
    }
        
    const static std::array<string, 2> esf_model_names;
  protected:
    inline double get_alpha(void) const {
        return alpha;
    }
    
    inline double sinc(double x) {
        return x == 0 ? 1 : sin(x)/x;
    }
    
    double alpha = 13.5;
    vector<double> w; // MTF correction weight
    static constexpr double missing = -1e7;
    bool apply_monotonic_filter = false;
};

#endif
