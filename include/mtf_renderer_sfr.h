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
#ifndef MTF_RENDERER_SFR_H
#define MTF_RENDERER_SFR_H

#include "mtf_renderer.h"
#include "common_types.h"

class Mtf_renderer_sfr : public Mtf_renderer {
  public:
    Mtf_renderer_sfr(const std::string& fname, 
      bool lpmm_mode=false, double pixel_size=1.0) 
      :  ofname(fname), lpmm_mode(lpmm_mode), pixel_size(pixel_size) {
      
    }
    
    void render(const vector<Block>& blocks) {
        FILE* fout = fopen(ofname.c_str(), "wt");
        for (size_t i=0; i < blocks.size(); i++) {
            for (size_t k=0; k < 4; k++) {
                const vector<double>& sfr = blocks[i].get_sfr(k);
                fprintf(fout, "%lf ", angle_reduce(atan2(-blocks[i].get_normal(k).x, blocks[i].get_normal(k).y)));
                for (size_t j=0; j < sfr.size(); j++) {
                    fprintf(fout, "%lf ", sfr[j]);
                }
                fprintf(fout, "\n");
            }
        }    
        fclose(fout);
    }
    
    string ofname;
    bool    lpmm_mode;
    double  pixel_size;

  private:
    double angle_reduce(double x) {
        double quad1 = fabs(fmod(x, M_PI/2.0));
        if (quad1 > M_PI/4.0) {
            quad1 = M_PI/2.0 - quad1;
        }
        quad1 = quad1 / M_PI * 180;
        return quad1;
    }
};

#endif
