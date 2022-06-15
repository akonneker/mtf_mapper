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

#include "include/mtf_renderer_annotate.h"
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>

static vector<double> calculate_edge_lengths(const Block& block) {
    vector<Point2d> base(4);
    vector<Point2d> n(4);
    for (size_t k=0; k < 4; k++) {
        size_t kp = (k + 1) % 4;
        
        base[k] = block.get_corner(k);
        Point2d dir(block.get_corner(kp) - base[k]);
        n[k] = Point2d(-dir.y, dir.x);
    }
    
    vector<double> lengths;
    int num_valid = 0;
    
    for (size_t k=0; k < 4; k++) {
        if (block.get_edge_valid(k)) {
            num_valid++;
            lengths.push_back(block.get_edge_length(k)*1.1);
        }
    }
    
    if (num_valid == 4) {
        lengths.resize(4, 0);
        for (size_t k=0; k < 4; k++) {
            size_t min_j = 0;
            double min_dist = 1e50;
            for (size_t j=0; j < 4; j++) {
                Point2d delta = block.get_edge_centroid(k) - base[j];
                double dist = delta.dot(n[j])/(n[j].dot(n[j]));
                if (fabs(dist) < min_dist) {
                    min_dist = fabs(dist);
                    min_j = j;
                }
            }
            lengths[min_j] = cv::norm(n[k]);
        }
    }
    
    return lengths;
}

static double estimate_font_scale(const cv::Mat& img, const char* buffer, double edge_length, double mean_edge_length) {
    int baseline = 0;
    double font_scale = 0.5;
    double font_thickness_scale = 1.0;
    
    // initial estimate of font width
    cv::Size ts = cv::getTextSize(buffer, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness_scale, &baseline);
    
    if (edge_length > 0) {
    
        double wr = ts.width / double(img.cols);
        constexpr double target_wr = 0.016; // magic number for a pleasing font-size-to-image-width ratio
        
        font_scale = std::max(0.5, 0.5 * target_wr / wr);
        // if this is a vertical image, increase the font scale proportionally to the aspect ratio
        if (double(img.rows)/img.cols > 1.55) {
            font_scale = std::max(0.5, 0.5 * target_wr / wr * double(img.rows)/img.cols);
        }
        
        font_thickness_scale = font_scale / 0.5;
        ts = cv::getTextSize(buffer, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness_scale, &baseline);
        
        // first try to adjust size based on mean edge length (actually, mean of the three shortest edges)
        while (ts.width > mean_edge_length*0.8 && font_scale > 0.5) {
            font_scale = std::max(0.5, font_scale*0.95);
            font_thickness_scale = std::max(1.0, font_scale / 0.5);
            ts = cv::getTextSize(buffer, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness_scale, &baseline);
        }
        
        // now reduce the font scale further if this specific edge is short
        while (ts.width > edge_length*0.65 && font_scale > 0.5) {
            font_scale = std::max(0.5, font_scale*0.95);
            font_thickness_scale = std::max(1.0, font_scale / 0.5);
            ts = cv::getTextSize(buffer, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness_scale, &baseline);
        }
        
        // lastly we deal with really short edges, where we allow font scales below 0.5
        while (edge_length < 50 && ts.width > edge_length*0.65 && font_scale > 0.35) {
            font_scale = std::max(0.35, font_scale*0.95);
            font_thickness_scale = std::max(1.0, font_scale / 0.5);
            ts = cv::getTextSize(buffer, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness_scale, &baseline);
        }
        
    }
    
    return font_scale;
}

void Mtf_renderer_annotate::render(const vector<Block>& blocks) {
    vector<vector<double>> all_edge_lengths(blocks.size());
    vector<double> edge_length_stats;
    for (size_t i=0; i < blocks.size(); i++) {
        all_edge_lengths[i] = calculate_edge_lengths(blocks[i]);
        for (double v: all_edge_lengths[i]) {
            if (v > 0) {
                edge_length_stats.push_back(v);
            }
        }
    }
    sort(edge_length_stats.begin(), edge_length_stats.end());
    double quantile_edge_length = 0;
    if (edge_length_stats.size() > 0) {
        quantile_edge_length = edge_length_stats[(size_t)floor(edge_length_stats.size()/5)];
    }
    char buffer[] = "8.88";
    double initial_font_scale = estimate_font_scale(out_img, buffer, quantile_edge_length, quantile_edge_length);
    
    for (size_t i=0; i < blocks.size(); i++) {
        vector<double>& edge_lengths = all_edge_lengths[i];
        vector<double> ordered_edge_lengths(edge_lengths);
        sort(ordered_edge_lengths.begin(), ordered_edge_lengths.end());
        double mean_edge_length = (
            ordered_edge_lengths[0] + 
            ordered_edge_lengths[1 % ordered_edge_lengths.size()] + 
            ordered_edge_lengths[2 % ordered_edge_lengths.size()]
        )/3.0;
        for (size_t k=0; k < 4; k++) {
            double val = blocks[i].get_mtf50_value(k);
            if (val > 0) {
                Point2d cent = blocks[i].get_edge_centroid(k);
                if (cent.x > 1 && cent.y > 1) {
                
                    // if this is a short edge, adjust font_scale to match
                    double font_scale = initial_font_scale;
                    if (edge_lengths[k] < quantile_edge_length * 0.9 || edge_lengths[k] < 50) {
                        char buffer[30];
                        double freq_scale = lpmm_mode ? pixel_size : 1.0;
                        if (val < 1.0) {
                            if (lpmm_mode) {
                                sprintf(buffer, "%.1lf", val*freq_scale);
                            } else {
                                sprintf(buffer, "%.2lf", val*freq_scale);
                            }
                        } else {
                            sprintf(buffer, "N/A");
                        }
                        font_scale = std::min(font_scale, estimate_font_scale(out_img, buffer, edge_lengths[k], mean_edge_length));
                    }
                    write_number(out_img, lrint(cent.x), lrint(cent.y), val, blocks[i].get_quality(k), font_scale);
                }
            }
        }
    }
    if (jpeg_output) {
        vector<int> parms = {cv::IMWRITE_JPEG_QUALITY, 98};
        imwrite(ofname + ".jpg", out_img, parms);
    } else {
        imwrite(ofname + ".png", out_img);
    }
    
}

void Mtf_renderer_annotate::write_number(cv::Mat& img, int px, int py, double val, double quality, double font_scale) {
    char buffer[10];
    
    double freq_scale = lpmm_mode ? pixel_size : 1.0;
    
    if (val < 1.0) {
        if (lpmm_mode) {
            sprintf(buffer, "%.1lf", val*freq_scale);
        } else {
            sprintf(buffer, "%.2lf", val*freq_scale);
        }
    } else {
        sprintf(buffer, "N/A");
    }
    
    int baseline = 0;
    double font_thickness_scale = std::max(1.0, font_scale / 0.5);
    
    cv::Size ts = cv::getTextSize(buffer, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness_scale, &baseline);
    
    cv::Point to(-ts.width/2,  ts.height/2);
    to.x += px;
    to.y += py;
    
    cv::putText(img, buffer, to, 
        cv::FONT_HERSHEY_SIMPLEX, font_scale, 
        CV_RGB(20, 20, 20), 2*font_thickness_scale + 0.5, cv::LINE_AA
    );
    
    cv::Scalar col = CV_RGB(0, 255, 255);
    if (quality < 0.8) {
        col = CV_RGB(255, 255, 0);
    }
    if (quality <= 0.2 || val == 1.0) {
        col = CV_RGB(255, 0, 0);
    }
    
    cv::putText(img, buffer, to, 
        cv::FONT_HERSHEY_SIMPLEX, font_scale, 
        col, font_thickness_scale > 1 ? font_thickness_scale + 0.5 : font_thickness_scale, cv::LINE_AA
    );
    
}