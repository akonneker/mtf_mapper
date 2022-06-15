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
#ifndef CA_CORE_H
#define CA_CORE_H

#include "include/mtf_core.h"

class Ca_core {
  public:
    Ca_core(Mtf_core& mtf_core)
    : mtf_core(mtf_core) {
    }
    
    void set_rgb_channels(vector<cv::Mat> in_channels) {
        channels = in_channels;
    }
    
    void calculate_ca(Block& block);
    void set_allow_all_edges(void);
    
    Mtf_core& mtf_core;
    vector<cv::Mat> channels;
    
  private:
    void extract_rgb_lsf_bayer(Block& block, const cv::Mat& img, const cv::Mat& bayer_img,
        vector<vector<double>>& red_lsf, vector<vector<double>>& green_lsf, vector<vector<double>>& blue_lsf);
        
    void extract_rgb_lsf(Block& block, const cv::Mat& img, const vector<cv::Mat>& channels,
        vector<vector<double>>& red_lsf, vector<vector<double>>& green_lsf, vector<vector<double>>& blue_lsf);
        
    bool allow_all_edges = false;
};

#endif
