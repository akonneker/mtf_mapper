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
#ifndef MTF_RENDERER_ANNOTATE_H
#define MTF_RENDERER_ANNOTATE_H

#include "mtf_renderer.h"
#include "common_types.h"
#include "include/srgb_render.h"

class Mtf_renderer_annotate : public Mtf_renderer {
  public:
    Mtf_renderer_annotate(const cv::Mat& in_img, const std::string& fname,
      bool lpmm_mode, double pixel_size, bool jpeg_output) 
      : img(in_img), ofname(fname), lpmm_mode(lpmm_mode), 
        pixel_size(pixel_size), jpeg_output(jpeg_output) {
      
          out_img = Srgb_render::linear_to_sRGB(in_img);
    }
    
    void render(const vector<Block>& blocks);
    void write_number(cv::Mat& img, int px, int py, double val, double quality, double font_scale);
    
    const cv::Mat& img;
    cv::Mat out_img;
    string ofname;
    bool    lpmm_mode;
    double  pixel_size;
    bool jpeg_output = false;
};

#endif
