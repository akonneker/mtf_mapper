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
#ifndef MTF_RENDERER_LENSPROFILE_H
#define MTF_RENDERER_LENSPROFILE_H

#include "include/logger.h"
#include "mtf_renderer.h"
#include "common_types.h"
#include "include/loess_fit.h"
#include "include/savitzky_golay_tables.h"

class Mtf_renderer_lensprofile : public Mtf_renderer {
  public:
    Mtf_renderer_lensprofile(
        const std::string& img_filename,
        const std::string& wdir, 
        const std::string& prof_fname, 
        const std::string& gnuplot_binary,
        const cv::Mat& img,
        const vector<double>& in_resolution,
        int gnuplot_width,
        bool lpmm_mode=false, double pixel_size=1.0) 
      :  Mtf_renderer(img_filename),
         wdir(wdir), prname(prof_fname), 
         gnuplot_binary(gnuplot_binary), img(img), 
         lpmm_mode(lpmm_mode), pixel_size(pixel_size),
         gnuplot_failure(false),
         in_resolution(in_resolution), gnuplot_width(gnuplot_width) {
    }
    
    void set_sparse_chart(bool s) {
        sparse_chart = s;
    }
    
    void set_fixed_size(bool f) {
        fixed_size = f;
    }
    
    void render(const vector<Block>& blocks);
    
    void lsfit(const vector<Ordered_point>& in_data, vector<Ordered_point>& recon, 
        vector<Ordered_point>& spread, int recon_samples=64); 
        
  private:
    double angle_reduce(double x) {
        double quad1 = fabs(fmod(x, M_PI/2.0));
        if (quad1 > M_PI/4.0) {
            quad1 = M_PI/2.0 - quad1;
        }
        quad1 = quad1 / M_PI * 180;
        return quad1;
    }
    
    inline double polyeval(double x, const vector<double>& a) const {
        double px = x;
        double s = a[0];
        for (size_t i=1; i < a.size(); i++) {
            s += a[i]*px;
            px *= x;
        }
        return s;
    }
    
    string wdir;
    string prname; 
    string pfname;
    string gnuplot_binary;
    const cv::Mat& img;            
    bool    lpmm_mode;
    double  pixel_size;
    bool gnuplot_failure;
    
    vector<double> in_resolution;
    int gnuplot_width;
    bool sparse_chart = false;
    bool fixed_size = false;
};

#endif
