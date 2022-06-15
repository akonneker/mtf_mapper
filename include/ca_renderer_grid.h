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
#ifndef CA_RENDERER_GRID_H
#define CA_RENDERER_GRID_H

#include "mtf_renderer.h"
#include "common_types.h"

class Ca_renderer_grid : public Mtf_renderer {
  public:
    Ca_renderer_grid(const std::string& img_filename, const std::string& wdir, 
        const std::string& fname, const std::string& gnuplot_binary, const cv::Mat& img, 
        int gnuplot_width, bool lpmm_mode, double pixel_size, bool fraction_mode, 
        bool allow_all_edges);
    
    void render(const vector<Block>& blocks);
    
    
    void set_sparse_chart(bool mode) {
        sparse_mode = mode;
    }
    
    string img_filename;
    string wdir;
    string fname;
    string gnuplot_binary;
    cv::Mat img; // Maybe we can drop this one?
    int gnuplot_width;
    bool lpmm_mode;
    double pixel_size;
    bool fraction_mode;
    
    Point2d img_centre;
    cv::Size img_dims;
    
    bool sparse_mode;
    bool allow_all_edges;
    
    vector<cv::Mat> grid_coarse;
    vector<cv::Mat> grid_fine;
};

#endif
