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

#include "include/esf_sampler_line.h"

void Esf_sampler_line::sample(Edge_model& edge_model, vector<Ordered_point>& local_ordered, 
    const map<int, scanline>& scanset, double& edge_length,
    const cv::Mat& geom_img, const cv::Mat& sampling_img,
    Bayer::cfa_mask_t cfa_mask) {
    
    cfa_mask = cfa_mask == Bayer::DEFAULT ? default_cfa_mask : cfa_mask;
    
    double max_along_edge = -1e50;
    double min_along_edge = 1e50;
    
    for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); ++it) {
        int y = it->first;
        if (y < border_width || y > geom_img.rows-1-border_width) continue;
        int rowcode = (y & 1) << 1;
        
        for (int x=it->second.start; x <= it->second.end; ++x) {
            
            if (x < border_width || x > geom_img.cols-1-border_width) continue;
            
            int code = 1 << ( (rowcode | (x & 1)) ^ 3 );
            if ((code & cfa_mask) == 0) continue;
            
            Point2d d = Point2d(x, y) - edge_model.get_centroid();
            double perp = d.ddot(edge_model.get_normal()); 
            double par = d.ddot(edge_model.get_direction());
            if (fabs(perp) < max_dot && fabs(par) < max_edge_length) {
                local_ordered.push_back(Ordered_point(perp, sampling_img.at<uint16_t>(y,x) ));
                max_along_edge = max(max_along_edge, par);
                min_along_edge = min(min_along_edge, par);
            }
        }
    }
        
    edge_length = max_along_edge - min_along_edge;
}

