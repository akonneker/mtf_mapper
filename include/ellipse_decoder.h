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
#ifndef ELLIPSE_DECODER_H
#define ELLIPSE_DECODER_H

#include "include/logger.h"
#include "include/common_types.h"
#include "include/ellipse.h"
#include "include/component_labelling.h"


class Ellipse_decoder {
  public:
    
    Ellipse_decoder(const Ellipse_detector& e, const cv::Mat& img) 
      : code(-1), valid(false), ratio(e.minor_axis/e.major_axis) {
        
        extract_code(e, img);
    }
    
    void extract_code(const Ellipse_detector&e, const cv::Mat& img) {
        
        // collect histogram stats inside ellipse
        map<int, int> histo;
        double total = 0;
        for (map<int, scanline>::const_iterator it=e.scanset.begin(); it != e.scanset.end(); it++) {
            int y = it->first;
            for (int x=it->second.start; x <= it->second.end; x++) {
                int val = img.at<uint16_t>(y, x);
            
                map<int, int>::iterator hi=histo.find(val);
                if (hi == histo.end()) {
                    histo[val] = 0;
                }
                histo[val]++;
                total++;
            }
        }
        double sum = 0;
        // compute Otsu threshold
        for (map<int, int>::iterator it=histo.begin(); it != histo.end(); it++) {
            sum += it->first * it->second;
        }
        double sum_b = 0;
        double w_b = 0;
        double max = 0;
        double thresh1 = 0;
        double thresh2 = 0;
        for (map<int, int>::iterator it=histo.begin(); it != histo.end(); it++) {
            w_b += it->second;
            double w_f = total - w_b;
            if (w_f == 0) break;
            
            sum_b += it->first * it->second;
            double m_b = sum_b / w_b;
            double m_f = (sum - sum_b) / w_f;
            double between = w_b * w_f * (m_b - m_f) * (m_b - m_f);
            if (between >= max) {
                thresh1 = it->first;
                if (between > max) {
                    thresh2 = it->first;
                }
                max = between;
            }
        }
        int otsu = lrint((thresh1 + thresh2) * 0.5);
        
        double cs = cos(e.angle);
        double ss = sin(e.angle);
        
        // take 30 samples along inner track
        int steps = 30;
        int ones = 0;
        for (int i=0; i < steps; i++) {
            double theta = i*2.0*M_PI/double(steps);
            
            double synth_x = 0.3*e.major_axis * cos(theta);
            double synth_y = 0.3*e.minor_axis * sin(theta);
            double px = cs*synth_x - ss*synth_y + e.centroid_x;
            double py = ss*synth_x + cs*synth_y + e.centroid_y;
            
            int bit = img.at<uint16_t>(lrint(py), lrint(px)) > otsu ? 1 : 0;
            ones += bit;
        }
        double inner_ratio = ones/double(steps);
        
        steps = 50;
        ones = 0;
        for (int i=0; i < steps; i++) {
            double theta = i*2.0*M_PI/double(steps);
            double synth_x = 0.5*e.major_axis * cos(theta);
            double synth_y = 0.5*e.minor_axis * sin(theta);
            double px = cs*synth_x - ss*synth_y + e.centroid_x;
            double py = ss*synth_x + cs*synth_y + e.centroid_y;
            
            int bit = img.at<uint16_t>(lrint(py), lrint(px)) > otsu ? 1 : 0;
            ones += bit;
        }
        double outer_ratio = ones/double(steps);
        
        steps = 10;
        ones = 0;
        for (int i=0; i < steps && !ones; i++) {
            double theta = i*2.0*M_PI/double(steps);
            
            double synth_x = 0.07*e.major_axis * cos(theta);
            double synth_y = 0.07*e.minor_axis * sin(theta);
            double px = cs*synth_x - ss*synth_y + e.centroid_x;
            double py = ss*synth_x + cs*synth_y + e.centroid_y;
            
            int bit = img.at<uint16_t>(lrint(py), lrint(px)) > otsu ? 1 : 0;
            ones += bit;
        }
        
        int ix = lrint(e.centroid_x);
        int iy = lrint(e.centroid_y);
        if (!ones && ix > 0 && iy > 0 && ix < img.cols - 1 && iy < img.rows - 1 &&
            img.at<uint16_t>(iy, ix) < otsu) {
            
            code = -1;
            valid = false;
            return;
        }
        
        int inner_code = lrint(inner_ratio*3);
        int outer_code = lrint(outer_ratio*5);
        int final_code = outer_code*4 + inner_code;
        code = final_code;
        valid = true;
        
        
        if (e.fg_fraction > 0.9999) {
            logger.debug("%s\n", "ellipse too solid, cannot be a valid code");
            valid = false;
        }
    }
    
    int code;
    bool valid;
    
    double ratio;
};

#endif
