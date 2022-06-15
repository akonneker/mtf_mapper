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
#ifndef SNR_H
#define SNR_H

class Snr {
  public:
    Snr(double dark_mean=0, double dark_sd=1, double bright_mean=1, double bright_sd=1) 
    : dark_mean(dark_mean), dark_sd(dark_sd), bright_mean(bright_mean), 
      bright_sd(bright_sd), contrast_val(bright_mean - dark_mean), oversampling_val(0) {
    }
    
    inline double mean_cnr(void) const {
        double denom = sqrt(0.5*dark_sd*dark_sd + 0.5*bright_sd*bright_sd);
        return std::min(10000.0, contrast_val / std::max(0.01, denom));
    }
    
    inline double dark_cnr(void) const {
        return std::min(10000.0, contrast_val / std::max(0.01, dark_sd));
    }
    
    inline double bright_cnr(void) const {
        return std::min(10000.0, contrast_val / std::max(0.01, bright_sd));
    }
    
    inline double dark_snr(void) const {
        return std::min(10000.0, dark_mean / std::max(0.01, dark_sd));
    }
    
    inline double bright_snr(void) const {
        return std::min(10000.0, bright_mean / std::max(0.01, bright_sd));
    }
    
    inline double contrast(void) const {
        return contrast_val;
    }
    
    inline double oversampling(void) const {
        return oversampling_val;
    }
    
    void set_oversampling(double val) {
        oversampling_val = val;
    }
    
  private:
    double dark_mean;
    double dark_sd;
    double bright_mean;
    double bright_sd;
    double contrast_val;
    double oversampling_val;
};

#endif
