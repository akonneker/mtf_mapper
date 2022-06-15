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
#ifndef DISPLAY_PROFILE_H
#define DISPLAY_PROFILE_H

#include <stdio.h>
#include <string>
using std::string;

#include <vector>
using std::vector;
using std::pair;

#include <array>

#include <opencv2/core/core.hpp>

class Display_profile {
  public:
    Display_profile(void);
    Display_profile(const vector<double>& gparm, const vector<double>& luminance_weights);
    Display_profile(const vector< pair<uint16_t, uint16_t> >& gtable, const vector<double>& luminance_weights);
    
    void force_linear(void) { is_linear = true; }
    void force_sRGB(void);
    cv::Mat to_luminance(const cv::Mat& img);
    vector<cv::Mat> to_linear_rgb(const cv::Mat& img);
    
  private:
    void render_parametric(const vector<double>& gparm);
  
    std::array<uint16_t, 65536> lut;
    vector<double> luminance_weights;
    bool is_linear = false;
    
};

#endif
