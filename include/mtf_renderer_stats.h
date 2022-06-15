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
#ifndef MTF_RENDERER_STATS_H
#define MTF_RENDERER_STATS_H

#include "include/logger.h"
#include "mtf_renderer.h"
#include "common_types.h"

#include <vector>
using std::vector;

#include <cmath>

class Mtf_renderer_stats : public Mtf_renderer {
  public:
    Mtf_renderer_stats(bool lpmm_mode=false, double pixel_size=1.0)
    : pixel_size(lpmm_mode ? pixel_size : 1)  {
      
    }
    
    void render(const vector<Block>& blocks) {
        
        if (blocks.size() == 0) {
            return;
        }
    
        vector<double> unfiltered;
        vector<double> filtered;
        for (size_t i=0; i < blocks.size(); i++) {
            if (!blocks[i].valid) continue;
            for (size_t k=0; k < 4; k++) {
            
                double val = blocks[i].get_mtf50_value(k);
                
                if (val == 1.0) continue;
                
                unfiltered.push_back(val);
                
                if (blocks[i].get_quality(k) >= 0.5) {
                    filtered.push_back(val);
                }
            }
        }
        
        if (unfiltered.size() < 1) {
            return;
        }
        
        print_stats(unfiltered, filtered);
    }
    
    void render(const vector<Mtf_profile_sample>& samples) {
        
        if (samples.size() == 0) {
            return;
        }
    
        vector<double> unfiltered;
        vector<double> filtered;
        for (size_t i=0; i < samples.size(); i++) {
            
            double val = samples[i].mtf;
            
            if (val == 1.0) continue;
            
            unfiltered.push_back(val);
            
            if (samples[i].quality >= 0.5) {
                filtered.push_back(val);
            }
        }
        
        if (unfiltered.size() < 2) {
            return;
        }
        
        print_stats(unfiltered, filtered);
    }
    
    void print_stats(vector<double>& unfiltered, vector<double>& filtered) {
        sort(filtered.begin(), filtered.end());
        sort(unfiltered.begin(), unfiltered.end());
        
        logger.info("    Quantiles ->                   %5d%% %5d%% %5d%% %5d%% %5d%%\n", 5, 25, 50, 75, 95);
        logger.info("Statistics on all edges:           %5.4lf %5.4lf %5.4lf %5.4lf %5.4lf  (total=%u)\n", 
            quantile(unfiltered, 0.05)*pixel_size,
            quantile(unfiltered, 0.25)*pixel_size,
            quantile(unfiltered, 0.5)*pixel_size,
            quantile(unfiltered, 0.75)*pixel_size,
            quantile(unfiltered, 0.95)*pixel_size,
            (unsigned int)unfiltered.size()
        );
        
        if (filtered.size() < 2) {
            return;
        }
        
        logger.info("Statistics on all filtered edges : %5.4lf %5.4lf %5.4lf %5.4lf %5.4lf  (total=%u)\n", 
            quantile(filtered, 0.05)*pixel_size,
            quantile(filtered, 0.25)*pixel_size,
            quantile(filtered, 0.5)*pixel_size,
            quantile(filtered, 0.75)*pixel_size,
            quantile(filtered, 0.95)*pixel_size,
            (unsigned int)filtered.size()
        );
    }
    
  private:
    double quantile(const vector<double>& d, double q) {
        size_t idx = (int)floor(d.size()*q);
        return d[idx];
    }

    double  pixel_size;
};

#endif
