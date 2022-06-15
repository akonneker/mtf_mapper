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
#ifndef MTF_RENDERER_BLOCKS_H
#define MTF_RENDERER_BLOCKS_H

#include "mtf_renderer.h"
#include "common_types.h"

class Mtf_renderer_blocks : public Mtf_renderer {
  public:
    Mtf_renderer_blocks(const cv::Mat& in_img, const std::string& fname) 
      : img(in_img), ofname(fname) {
      
          out_img = img.clone();
    }
    
    void render(const vector<Block>& blocks) {
        cv::Point* pts = new cv::Point[4];;
        for (size_t i=0; i < blocks.size(); i++) {
            int npts=4;
            const Mrectangle& rect = blocks[i].rect;
            for (size_t k=0; k < 4; k++) { 
                pts[k] = cv::Point(rect.corners[k].x, rect.corners[k].y);
            }
            cv::polylines(out_img, (const cv::Point**)&pts, &npts, 1, true, CV_RGB(65535, 65535, 65535 ), 1, CV_AA);
        }    
        delete [] pts;
        
        imwrite(ofname, out_img);
    }
    
    const cv::Mat& img;
    cv::Mat out_img;
    string ofname;
};

#endif
