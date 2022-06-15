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

#include "include/logger.h"
#include "include/demosaic.h"
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include "include/threadpool.h"
#include <string>
using std::string;
#include <vector>
using std::vector;

void simple_demosaic_green(cv::Mat& cvimg, cv::Mat& rawimg, bool unbalanced_scene, bool swap_diag=false);
void simple_demosaic_redblue(cv::Mat& cvimg, cv::Mat& rawimg, Bayer::bayer_t bayer, Bayer::cfa_pattern_t cfa_pattern);

void simple_demosaic(cv::Mat& cvimg, cv::Mat& rawimg, Bayer::cfa_pattern_t cfa_pattern, Bayer::bayer_t bayer, bool unbalanced_scene) {
    // switch on Bayer subset
    if (bayer == Bayer::GREEN) {
        simple_demosaic_green(cvimg, rawimg, unbalanced_scene, cfa_pattern == Bayer::GRBG || cfa_pattern == Bayer::GBRG);
    } else {
        if (bayer == Bayer::RED || bayer == Bayer::BLUE) {
            simple_demosaic_redblue(cvimg, rawimg, bayer, cfa_pattern);
        } else {
            logger.error("%s\n", "Fatal error: Unknown bayer subset requested. Aborting");
            exit(-1);
        }
    }
}

static void match(const vector<int>& source, const vector<int>& target, vector<int>& m) {
    int cur_ix = 0;
    for (int i=0; i < 65536; i++) {
        while(source[cur_ix] < target[i]) cur_ix++;
        if (cur_ix > 0) {
            if ( (source[cur_ix] - target[i]) < (target[i] - source[cur_ix - 1]) ) {
                m[i] = cur_ix;
            } else {
                m[i] = cur_ix - 1;
            }
        } else {
            m[i] = cur_ix;
        }
    }
}

void simple_demosaic_green(cv::Mat& cvimg, cv::Mat& rawimg, bool unbalanced_scene, bool swap_diag) {
    rawimg = cvimg.clone();
    
    // target subsets
    int tss1 = 0;
    int tss2 = 3;
    if (swap_diag) {
        tss1 = 1;
        tss2 = 2;
    }
    
    if (!unbalanced_scene) {
        logger.debug("%s\n", "Green Bayer subset specified, performing quick-and-dirty balancing of green channels");
        vector < vector<int> > hist(4, vector<int>(65536, 0));
        for (size_t row=0; row < (size_t)cvimg.rows; row++) {
            for (int col=0; col < cvimg.cols; col++) {
                int val = cvimg.at<uint16_t>(row, col);
                int subset = ((row & 1) << 1) | (col & 1);
                hist[subset][val]++;
            }
        }
        // convert histograms to cumulative histograms
        for (int subset=0; subset < 4; subset++) {
            int acc = 0;
            for (size_t i=0; i < hist[subset].size(); i++) {
                acc += hist[subset][i];
                hist[subset][i] = acc;
            }
        }
        
        int from_ss = 1;
        int to_ss = 2;
        if (swap_diag) {
            from_ss = 0;
            to_ss = 3;
        }
        
        // NB: cumulative histogram totals must match
        // this can be distorted on cropped images
        if (hist[from_ss].back() > hist[to_ss].back()) {
            for (auto& v: hist[from_ss]) {
                v = int64_t(v)*int64_t(hist[to_ss].back()) / int64_t(hist[from_ss].back());
            }
        }
        
        vector<int> m(65536, 0);
        match(hist[to_ss], hist[from_ss], m);
        for (size_t row=0; row < (size_t)cvimg.rows; row++) {
            for (int col=0; col < cvimg.cols; col++) {
                int subset = ((row & 1) << 1) | (col & 1);
                if (subset == from_ss) { // TODO: optimize access pattern
                    int val = cvimg.at<uint16_t>(row, col);
                    cvimg.at<uint16_t>(row, col) = m[val];
                }
            }
        }
    } else {
        logger.debug("%s\n", "ROI mode, not performing G1/G2 Bayer subset matching");
    }

    ThreadPool& tp = ThreadPool::instance();
    vector<std::future<void>> futures;
    size_t n_blocks = std::max(1, int(cvimg.rows < 50 ? 1 : tp.size()));
    size_t block_size = cvimg.rows / n_blocks;
    for (size_t block = 0; block < n_blocks; block++) {
        futures.emplace_back(
            tp.enqueue([&, block] {
                size_t start_row = std::min(4 + block*block_size, size_t(cvimg.rows - 1 - 4));
                size_t end_row = std::min(4 + (block + 1) * block_size, size_t(cvimg.rows - 4));
                uint16_t* cv_base = (uint16_t*)cvimg.data;
                for (size_t row = start_row; row < end_row; row++) {
                    for (int col = 4; col < cvimg.cols - 4; col++) {
                        int subset = ((row & 1) << 1) | (col & 1);
                        uint64_t cidx = row * cvimg.cols + col;
                        uint16_t* cvp = cv_base + cidx;
                        if (subset == tss1 || subset == tss2) {

                            double hgrad = fabs(double(int32_t(cvp[-3]) + 3 * int32_t(cvp[-1]) - 3 * int32_t(cvp[1]) - int32_t(cvp[3])));
                            double vgrad = fabs(double(int32_t(cvp[-3 * cvimg.cols]) + 3 * int32_t(cvp[-cvimg.cols]) - 3 * int32_t(cvp[cvimg.cols]) - int32_t(cvp[3 * cvimg.cols])));

                            if (max(hgrad, vgrad) < 1 || fabs(hgrad - vgrad) / max(hgrad, vgrad) < 0.1) {
                                *cvp = (int32_t(cvp[-cvimg.cols]) + int32_t(cvp[cvimg.cols]) + int32_t(cvp[-1]) + int32_t(cvp[1])) / 4;
                            }
                            else {
                                double l = (hgrad * hgrad + vgrad * vgrad);
                                if (hgrad > vgrad) {
                                    l = hgrad * hgrad / l;
                                    if (l > 0.92388 * 0.92388) { // more horizontal than not
                                        *cvp = (int32_t(cvp[-cvimg.cols]) + int32_t(cvp[cvimg.cols])) / 2;
                                    }
                                    else { // in between, blend it
                                        *cvp = (2 * (int32_t(cvp[-cvimg.cols]) + int32_t(cvp[cvimg.cols])) + (int32_t(cvp[-1]) + int32_t(cvp[1]))) / 6.0;
                                    }
                                }
                                else {
                                    l = vgrad * vgrad / l;
                                    if (l > 0.92388 * 0.92388) {
                                        *cvp = (int32_t(cvp[-1]) + int32_t(cvp[1])) / 2;
                                    }
                                    else {
                                        *cvp = (2 * (int32_t(cvp[-1]) + int32_t(cvp[1])) + (int32_t(cvp[-cvimg.cols]) + int32_t(cvp[cvimg.cols]))) / 6.0;
                                    }
                                }
                            }
                        }
                    }
                }
            })
        );
    }
    for (size_t i = 0; i < futures.size(); i++) {
        futures[i].wait();
    }
    
    
    // since the relative row/column offsets are hardcoded, just swap the subsets
    if (swap_diag) {
        std::swap(tss1, tss2);
    }
    
    // pad out the border using three-sample averaging; this leaves 
    // some zippering, but this should not matter too much downstream
    for (size_t row=0; row < (size_t)cvimg.rows; row++) {
        for (size_t col=0; col < 4; col++) {
            int subset = ((row & 1) << 1) | (col & 1);
            
            if (row == 0 || row == (size_t)(cvimg.rows-1)) {
                if (subset == tss1) { 
                    cvimg.at<uint16_t>(row, col) = cvimg.at<uint16_t>(row, col+1);
                } else if (subset == tss2) {
                    cvimg.at<uint16_t>(row, col) = cvimg.at<uint16_t>(row, col-1);
                }
            } else {
                if (subset == tss1 || subset == tss2) { 
                    cvimg.at<uint16_t>(row, col) = (cvimg.at<uint16_t>(row-1, col) + cvimg.at<uint16_t>(row+1, col) + cvimg.at<uint16_t>(row, col+1))/3;
                }
            }
            
        }
        for (size_t col=size_t(cvimg.cols-4); col < size_t(cvimg.cols); col++) {
            int subset = ((row & 1) << 1) | (col & 1);
            
            if (row == 0 || row == (size_t)(cvimg.rows-1)) {
                if (subset == tss1) { 
                    cvimg.at<uint16_t>(row, col) = cvimg.at<uint16_t>(row, col < size_t(cvimg.cols-1) ? col+1 : col-1);
                } else if (subset == tss2) {
                    cvimg.at<uint16_t>(row, col) = cvimg.at<uint16_t>(row, col-1);
                }
            } else {
                if (subset == tss1 || subset == tss2) { 
                    cvimg.at<uint16_t>(row, col) = (cvimg.at<uint16_t>(row-1, col) + cvimg.at<uint16_t>(row+1, col) + cvimg.at<uint16_t>(row, col-1))/3;
                }
            }
        }
    }
    
    for (size_t col=4; col < (size_t)(cvimg.cols-4); col++) {
        for (size_t row=0; row < 4; row++) {
            int subset = ((row & 1) << 1) | (col & 1);
            if (subset == tss1 || subset == tss2) { 
                cvimg.at<uint16_t>(row, col) = (cvimg.at<uint16_t>(row, col+1) + cvimg.at<uint16_t>(row, col-1) + cvimg.at<uint16_t>(row+1, col))/3;
            }
        }
        for (size_t row=size_t(cvimg.rows-4); row < size_t(cvimg.rows); row++) {
            int subset = ((row & 1) << 1) | (col & 1);
            if (subset == tss1 || subset == tss2) { 
                cvimg.at<uint16_t>(row, col) = (cvimg.at<uint16_t>(row, col+1) + cvimg.at<uint16_t>(row, col-1) + cvimg.at<uint16_t>(row-1, col))/3;
            }
        }
    }
}

void simple_demosaic_redblue(cv::Mat& cvimg, cv::Mat& rawimg, Bayer::bayer_t bayer, Bayer::cfa_pattern_t cfa_pattern) {
    // no need to white balance, only one channel?
    rawimg = cvimg.clone();
    
    int h_ss = 0;
    int v_ss = 3;
    
    // first subset is the complement, i.e., the one to replace
    int first_subset = 3;
    if (bayer == Bayer::RED) {
        switch(cfa_pattern) {
        case Bayer::RGGB: first_subset = 3; h_ss = 1; v_ss = 2; break;
        case Bayer::BGGR: first_subset = 0; h_ss = 2; v_ss = 1; break;
        case Bayer::GRBG: first_subset = 2; h_ss = 0; v_ss = 3; break;
        case Bayer::GBRG: first_subset = 1; h_ss = 3; v_ss = 0; break;
        }
    } else { // blue, of course
        switch(cfa_pattern) {
        case Bayer::RGGB: first_subset = 0; h_ss = 2; v_ss = 1; break;
        case Bayer::BGGR: first_subset = 3; h_ss = 1; v_ss = 2; break;
        case Bayer::GRBG: first_subset = 1; h_ss = 3; v_ss = 0; break;
        case Bayer::GBRG: first_subset = 2; h_ss = 0; v_ss = 3; break;
        }
    }
    
    // interpolate the other channel (Red if we want blue, or Blue if we want red)
    ThreadPool& tp = ThreadPool::instance();
    vector<std::future<void>> futures;
    size_t n_blocks = std::max(1, int(cvimg.rows < 50 ? 1 : tp.size()));
    size_t block_size = cvimg.rows / n_blocks;
    for (size_t block = 0; block < n_blocks; block++) {
        futures.emplace_back(
            tp.enqueue([&, block] {
            size_t start_row = std::min(4 + block * block_size, size_t(cvimg.rows - 1 - 4));
            size_t end_row = std::min(4 + (block + 1) * block_size, size_t(cvimg.rows - 4));
            for (size_t row = start_row; row < end_row; row++) {
                for (int col = 4; col < cvimg.cols - 4; col++) {
                    int subset = ((row & 1) << 1) | (col & 1);
                    if (subset == first_subset) { // TODO: really should optimize access pattern
                        double d1grad = fabs(double((int32_t)cvimg.at<uint16_t>(row - 1, col - 1) - (int32_t)cvimg.at<uint16_t>(row + 1, col + 1)));
                        double d2grad = fabs(double((int32_t)cvimg.at<uint16_t>(row - 1, col + 1) - (int32_t)cvimg.at<uint16_t>(row + 1, col - 1)));

                        if (max(d1grad, d2grad) < 1 || fabs(d1grad - d2grad) / max(d1grad, d2grad) < 0.001) {
                            cvimg.at<uint16_t>(row, col) =
                                ((int32_t)cvimg.at<uint16_t>(row - 1, col - 1) +
                                (int32_t)cvimg.at<uint16_t>(row + 1, col + 1) +
                                    (int32_t)cvimg.at<uint16_t>(row - 1, col + 1) +
                                    (int32_t)cvimg.at<uint16_t>(row + 1, col - 1)) / 4;
                        }
                        else {
                            if (d1grad > d2grad) {
                                cvimg.at<uint16_t>(row, col) = ((int32_t)cvimg.at<uint16_t>(row - 1, col + 1) + (int32_t)cvimg.at<uint16_t>(row + 1, col - 1)) / 2;
                            }
                            else {
                                cvimg.at<uint16_t>(row, col) = ((int32_t)cvimg.at<uint16_t>(row - 1, col - 1) + (int32_t)cvimg.at<uint16_t>(row + 1, col + 1)) / 2;
                            }
                        }
                    }
                }
            }
        })
        );
    }
    for (size_t i = 0; i < futures.size(); i++) {
        futures[i].wait();
    }
    
    // interpolate the two green channels
    futures.clear();
    for (size_t block = 0; block < n_blocks; block++) {
        futures.emplace_back(
            tp.enqueue([&, block] {
                size_t start_row = std::min(4 + block * block_size, size_t(cvimg.rows - 1 - 4));
                size_t end_row = std::min(4 + (block + 1) * block_size, size_t(cvimg.rows - 4));
                for (size_t row = start_row; row < end_row; row++) {
                    for (int col = 4; col < cvimg.cols - 4; col++) {
                        int subset = ((row & 1) << 1) | (col & 1);

                        if (subset == h_ss || subset == v_ss) {

                            double hgrad = fabs(double((int32_t)cvimg.at<uint16_t>(row, col - 1) - (int32_t)cvimg.at<uint16_t>(row, col + 1)));
                            double vgrad = fabs(double((int32_t)cvimg.at<uint16_t>(row - 1, col) - (int32_t)cvimg.at<uint16_t>(row + 1, col)));

                            if (max(hgrad, vgrad) < 1 || fabs(hgrad - vgrad) / max(hgrad, vgrad) < 0.001) {
                                cvimg.at<uint16_t>(row, col) =
                                    ((int32_t)cvimg.at<uint16_t>(row - 1, col) +
                                    (int32_t)cvimg.at<uint16_t>(row + 1, col) +
                                        (int32_t)cvimg.at<uint16_t>(row, col - 1) +
                                        (int32_t)cvimg.at<uint16_t>(row, col + 1)) / 4;
                            }
                            else {
                                if (hgrad > vgrad) {
                                    cvimg.at<uint16_t>(row, col) = ((int32_t)cvimg.at<uint16_t>(row - 1, col) + (int32_t)cvimg.at<uint16_t>(row + 1, col)) / 2;
                                }
                                else {
                                    cvimg.at<uint16_t>(row, col) = ((int32_t)cvimg.at<uint16_t>(row, col - 1) + (int32_t)cvimg.at<uint16_t>(row, col + 1)) / 2;
                                }
                            }
                        }
                    }
                }
            })
        );
    }
    for (size_t i = 0; i < futures.size(); i++) {
        futures[i].wait();
    }
    
    // on the first pass, interpolate the green channels
    for (size_t row=1; row < (size_t)cvimg.rows-1; row++) {
        for (size_t col=0; col < 4; col++) {
            int subset = ((row & 1) << 1) | (col & 1);
            
            if (subset == h_ss || subset == v_ss) {
                if (subset == h_ss) {
                    cvimg.at<uint16_t>(row, col) = ((int32_t)cvimg.at<uint16_t>(row, col + (col > 0 ? -1 : 1) ) +
                        (int32_t)cvimg.at<uint16_t>(row, col+1) )/2;
                } else {
                    cvimg.at<uint16_t>(row, col) = ((int32_t)cvimg.at<uint16_t>(row-1, col) + (int32_t)cvimg.at<uint16_t>(row+1, col))/2;
                }
            }
        }
        for (size_t col=size_t(cvimg.cols-5); col < size_t(cvimg.cols); col++) {
            int subset = ((row & 1) << 1) | (col & 1);
            
            if (subset == h_ss || subset == v_ss) {
                if (subset == h_ss) {
                    cvimg.at<uint16_t>(row, col) = ((int32_t)cvimg.at<uint16_t>(row, col-1) +
                        (int32_t)cvimg.at<uint16_t>(row, col + (col < size_t(cvimg.cols - 1) ? 1 : -1)) )/2;
                } else {
                    cvimg.at<uint16_t>(row, col) = ((int32_t)cvimg.at<uint16_t>(row-1, col) + (int32_t)cvimg.at<uint16_t>(row+1, col))/2;
                }
            }
        }
    }
    
    for (size_t col=1; col < size_t(cvimg.cols - 1); col++) {
        for (size_t row=0; row <= 4; row++) {
            int subset = ((row & 1) << 1) | (col & 1);
            
            if (subset == h_ss || subset == v_ss) {
                if (subset == h_ss) {
                    cvimg.at<uint16_t>(row, col) = ((int32_t)cvimg.at<uint16_t>(row, col-1) + (int32_t)cvimg.at<uint16_t>(row, col+1))/2;
                } else {
                    cvimg.at<uint16_t>(row, col) = ((int32_t)cvimg.at<uint16_t>(row + (row > 0 ? -1 : 1), col) +
                        (int32_t)cvimg.at<uint16_t>(row+1, col) )/2;
                }
            }
        }
        for (size_t row=size_t(cvimg.rows-5); row < size_t(cvimg.rows); row++) {
            int subset = ((row & 1) << 1) | (col & 1);
            
            if (subset == h_ss || subset == v_ss) {
                if (subset == h_ss) {
                    cvimg.at<uint16_t>(row, col) = ((int32_t)cvimg.at<uint16_t>(row, col-1) + (int32_t)cvimg.at<uint16_t>(row, col + 1))/2;
                } else {
                    cvimg.at<uint16_t>(row, col) = ((int32_t)cvimg.at<uint16_t>(row-1, col) +
                        (int32_t)cvimg.at<uint16_t>(row + (row < size_t(cvimg.rows-1) ? 1 : -1), col) )/2;
                }
            }
        }
    }
    
    // second pass, interpolate the remaining red or blue channel
    for (size_t row=1; row < (size_t)cvimg.rows-1; row++) {
        for (size_t col=0; col < 4; col++) {
            int subset = ((row & 1) << 1) | (col & 1);
            if (subset == first_subset) {
                cvimg.at<uint16_t>(row, col) = 
                        ((int32_t)cvimg.at<uint16_t>(row-1, col) +
                         (int32_t)cvimg.at<uint16_t>(row+1, col) +
                         (int32_t)cvimg.at<uint16_t>(row, col+1)) / 3;
                
            }
        }
        for (size_t col=size_t(cvimg.cols-4); col < size_t(cvimg.cols); col++) {
            int subset = ((row & 1) << 1) | (col & 1);
            if (subset == first_subset) {
                cvimg.at<uint16_t>(row, col) = 
                        ((int32_t)cvimg.at<uint16_t>(row-1, col) +
                         (int32_t)cvimg.at<uint16_t>(row+1, col) +
                         (int32_t)cvimg.at<uint16_t>(row, col-1)) / 3;
                
            }
        }
    }
    for (size_t col=1; col < (size_t)cvimg.cols-1; col++) {
        for (size_t row=0; row < 4; row++) {
            int subset = ((row & 1) << 1) | (col & 1);
            if (subset == first_subset) {
                cvimg.at<uint16_t>(row, col) = 
                        ((int32_t)cvimg.at<uint16_t>(row, col+1) +
                         (int32_t)cvimg.at<uint16_t>(row+1, col) +
                         (int32_t)cvimg.at<uint16_t>(row, col-1)) / 3;
                
            }
        }
        
        for (size_t row=size_t(cvimg.rows-5); row < size_t(cvimg.rows); row++) {
            int subset = ((row & 1) << 1) | (col & 1);
            if (subset == first_subset) {
                cvimg.at<uint16_t>(row, col) = 
                        ((int32_t)cvimg.at<uint16_t>(row-1, col) +
                         (int32_t)cvimg.at<uint16_t>(row, col+1) +
                         (int32_t)cvimg.at<uint16_t>(row, col-1)) / 3;
                
            }
        }
    }
    
    // nearest-neighbour interpolate the corners
    int self = 3 - first_subset;
    int corner[4][2] = { {0,0}, {cvimg.cols-2, 0}, {0, cvimg.rows-2}, {cvimg.cols-2, cvimg.rows-2}};
    uint16_t self_v[4];
    
    for (int cnr=0; cnr < 4; cnr++) {
        self_v[cnr] = cvimg.at<uint16_t>(
            corner[cnr][1] + ((self >> 1) & 1), 
            corner[cnr][0] + (self & 1)
        );
    }
    
    for (int dr=0; dr <= 1; dr++) {
        for (int dc=0; dc <= 1; dc++) {
            for (int cnr=0; cnr < 4; cnr++) {
                cvimg.at<uint16_t>(corner[cnr][1] + dr, corner[cnr][0] + dc) = self_v[cnr];
            }
        }
    }
}

inline void hv_cross(cv::Mat& cvimg, int row0, int col0, int centre_val, double& topsum, double& botsum) {
    double I1 = cvimg.at<uint16_t>(row0, col0-2);
    double I4 = cvimg.at<uint16_t>(row0, col0+2);
    double I2 = cvimg.at<uint16_t>(row0-2, col0);
    double I3 = cvimg.at<uint16_t>(row0+2, col0);
    double Ii = centre_val;
    
    topsum += (I1 + I4 - I2 - I3)*(Ii - 0.5*I2 - 0.5*I3);
    botsum += (I1 + I4 - I2 - I3)*(I1 + I4 - I2 - I3);
}

void geometric_demosaic(cv::Mat& cvimg, cv::Mat& rawimg, [[maybe_unused]] int target_subset) {
    logger.debug("%s\n", "Bayer subset specified, performing geometric demosaic");
    rawimg = cvimg.clone();
    
    // TODO: WB does not seem to work very well on IQ180 images ...
    
    vector < vector<int> > hist(4, vector<int>(65536, 0));
    for (size_t row=0; row < (size_t)cvimg.rows; row++) {
        for (int col=0; col < cvimg.cols; col++) {
            int val = cvimg.at<uint16_t>(row, col);
            int subset = ((row & 1) << 1) | (col & 1);
            if (subset > 3 || subset < 0) {
                logger.debug("subset = %d\n", subset);
            }
            if (val < 0 || val > 65535) {
                logger.debug("val = %d\n", val);
            }
            hist[subset][val]++;
        }
    }
    // convert histograms to cumulative histograms
    for (int subset=0; subset < 4; subset++) {
        int acc = 0;
        for (size_t i=0; i < hist[subset].size(); i++) {
            acc += hist[subset][i];
            hist[subset][i] = acc;
        }
    }
    
    // find 90% centre
    vector<double> l5(4, 0);
    vector<double> u5(4, 0);
    vector<double> mean(4, 0);
    for (int subset=0; subset < 4; subset++) {
        int lower = 0;
        int upper = hist[subset].size() - 1;
        while (hist[subset][lower] < 0.05*hist[subset].back() && lower < upper) lower++;
        while (hist[subset][upper] > 0.95*hist[subset].back() && upper > lower) upper--;
        l5[subset] = lower;
        u5[subset] = upper;
        
        mean[subset] = upper - lower;
        
        logger.debug("subset %d: %d %d, mean=%lf\n", subset, lower, upper, mean[subset]);
    }
    const int targ = 1;
    for (int i=0; i < 4; i++) {
        if (i != targ) {
            mean[i] = mean[targ] / mean[i];
        }
    }
    
    double maderr = 0;
    int interp_count = 0;
    vector<int> ncount(13, 0);
    int outlier_count = 0;
    // bilnear interpolation to get rid op R and B channels?
    for (size_t row=8; row < (size_t)cvimg.rows-8; row++) {
        for (int col=8; col < cvimg.cols-8; col++) {
            int subset = ((row & 1) << 1) | (col & 1);
            
            if (subset == 0 || subset == 3) {
            
                // obtain initial estimate of interpolated value using modfied gradient interpolation
                int gi_estimate = 0;
            
                double hgrad = fabs(double(cvimg.at<uint16_t>(row, col-3) + 3*cvimg.at<uint16_t>(row, col-1) - 
                                    3*cvimg.at<uint16_t>(row, col+1) - cvimg.at<uint16_t>(row, col+3)));
                double vgrad = fabs(double(cvimg.at<uint16_t>(row-3, col) + 3*cvimg.at<uint16_t>(row-1, col) - 
                                    3*cvimg.at<uint16_t>(row+1, col) - cvimg.at<uint16_t>(row+3, col)));
                                    
                if (max(hgrad, vgrad) < 1 || fabs(hgrad - vgrad)/max(hgrad, vgrad) < 0.001) {
                    gi_estimate = 
                        (cvimg.at<uint16_t>(row-1, col) + 
                        cvimg.at<uint16_t>(row+1, col) + 
                        cvimg.at<uint16_t>(row, col-1) + 
                        cvimg.at<uint16_t>(row, col+1)) / 4;
                } else {
                    double l = sqrt(hgrad*hgrad + vgrad*vgrad);
                    if (hgrad > vgrad) {
                        l = hgrad / l;
                        if (l > 0.92388) { // more horizontal than not
                            gi_estimate = (cvimg.at<uint16_t>(row-1, col) + cvimg.at<uint16_t>(row+1, col)) / 2;
                        } else { // in between, blend it
                            gi_estimate = (2*(cvimg.at<uint16_t>(row-1, col) + cvimg.at<uint16_t>(row+1, col)) +
                                (cvimg.at<uint16_t>(row, col-1) + cvimg.at<uint16_t>(row, col+1)) ) / 6.0;
                        }
                    } else {
                        l = vgrad / l;
                        if (l > 0.92388) {
                            gi_estimate = (cvimg.at<uint16_t>(row, col-1) + cvimg.at<uint16_t>(row, col+1)) / 2;
                        } else {
                            gi_estimate = (2*(cvimg.at<uint16_t>(row, col-1) + cvimg.at<uint16_t>(row, col+1)) +
                                (cvimg.at<uint16_t>(row-1, col) + cvimg.at<uint16_t>(row+1, col)) ) / 6.0;
                        }    
                    }
                }
                
                double topsum = 0;
                double botsum = 0;
                
                double topsums[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
                double botsums[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
                
                // cheat a little
                //gi_estimate = cvimg.at<uint16_t>(row, col);
                
                hv_cross(cvimg, row, col, gi_estimate, topsums[0], botsums[0]);
                
                hv_cross(cvimg, row - 1, col, cvimg.at<uint16_t>(row - 1, col), topsums[1], botsums[1]);
                hv_cross(cvimg, row + 1, col, cvimg.at<uint16_t>(row + 1, col), topsums[2], botsums[2]);
                hv_cross(cvimg, row, col - 1, cvimg.at<uint16_t>(row, col - 1), topsums[3], botsums[3]);
                hv_cross(cvimg, row, col + 1, cvimg.at<uint16_t>(row, col + 1), topsums[4], botsums[4]);
                
                hv_cross(cvimg, row + 1, col - 2, cvimg.at<uint16_t>(row + 1, col - 2), topsums[5], botsums[5]);
                hv_cross(cvimg, row + 1, col + 2, cvimg.at<uint16_t>(row + 1, col + 2), topsums[6], botsums[6]);
                hv_cross(cvimg, row - 1, col - 2, cvimg.at<uint16_t>(row - 1, col - 2), topsums[7], botsums[7]);
                hv_cross(cvimg, row - 1, col + 2, cvimg.at<uint16_t>(row - 1, col + 2), topsums[8], botsums[8]);
                
                
                /*
                hv_cross(cvimg, row - 3, col, cvimg.at<uint16_t>(row - 3, col), topsums[1], botsums[1]);
                hv_cross(cvimg, row + 3, col, cvimg.at<uint16_t>(row + 3, col), topsums[2], botsums[2]);
                hv_cross(cvimg, row, col - 3, cvimg.at<uint16_t>(row, col - 3), topsums[3], botsums[3]);
                hv_cross(cvimg, row, col + 3, cvimg.at<uint16_t>(row, col + 3), topsums[4], botsums[4]);
                */
                
                // we have some idea what the correct structure is, meaning topsum[0], botsum[0] will
                // give us a pretty good idea of what the preferred value of a1 is
                // select from the other four sites the ones that are compatible
                
                double a1_0 = 0.25;
                if (botsums[0] > 2) {
                    a1_0 = topsums[0]/botsums[0];
                    
                    // TODO: we can reduce the weight of this one?
                    topsum += topsums[0];
                    botsum += botsums[0];
                }
                
                
                // zero neighbours appears to be random (i.e., caused by noise)
                // which estimate do we use in these cases?
                // perhaps GI is safer?
                
                int nneighbours = 0;
                for (int n=1; n <= 8; n++) {
                    double nweight = (n > 4) ? ((n > 8) ? 0.125 : 0.25) : 1;
                    if (botsums[n] > 0) {
                        double a = topsums[n]/botsums[n];
                        if (fabs(a - a1_0) < 0.5 && a*a1_0 > 0) { 
                            nweight = 1;
                            nneighbours++;
                        }
                    }
                    //nneighbours++;
                    topsum += nweight*topsums[n];
                    botsum += nweight*botsums[n];
                }
                ncount[nneighbours]++;
                
                double ndiff = 0;
                
                ndiff = max(ndiff, fabs(double(cvimg.at<uint16_t>(row, col-1)) - double(cvimg.at<uint16_t>(row, col+1))));
                ndiff = max(ndiff, fabs(double(cvimg.at<uint16_t>(row, col-1)) - double(cvimg.at<uint16_t>(row-1, col))));
                ndiff = max(ndiff, fabs(double(cvimg.at<uint16_t>(row, col-1)) - double(cvimg.at<uint16_t>(row+1, col))));
                ndiff = max(ndiff, fabs(double(cvimg.at<uint16_t>(row, col+1)) - double(cvimg.at<uint16_t>(row-1, col))));
                ndiff = max(ndiff, fabs(double(cvimg.at<uint16_t>(row, col+1)) - double(cvimg.at<uint16_t>(row+1, col))));
                ndiff = max(ndiff, fabs(double(cvimg.at<uint16_t>(row-1, col)) - double(cvimg.at<uint16_t>(row+1, col))));
                
                double a1 = 0.25;
                double a2 = 0.25;
                double a3 = 0.25;
                double a4 = 0.25;
                //if (ndiff > 0 && norm_botsum / ndiff > 0.2 && botsum > 0 /*botsum > 1314329*/) { // magical empirical values??
                if (botsum > 2) {
                    a1 = topsum / botsum;
                    a1 = max(-1.0, a1);
                    a1 = min(a1, 1.0);
                    a4 = a1;
                    a2 = a3 = 0.5 - a1;
                }
                
                
                
                // we can re-normalize the weights, thereby allowing a larger range?
                double wsum = a1 + a2 + a3 + a4;
                a1 /= wsum;
                a2 /= wsum;
                a3 /= wsum;
                a4 /= wsum;
                
                double interp = 
                    cvimg.at<uint16_t>(row, col-1)*a1 +
                    cvimg.at<uint16_t>(row, col+1)*a4 +
                    cvimg.at<uint16_t>(row-1, col)*a2 +
                    cvimg.at<uint16_t>(row+1, col)*a3;
                    
                interp = max(0.0, min(interp, 65535.0));
                
                // difference between two interpolants should be small, or we fall back on the simpler gradient interpolant
                if (fabs(interp - gi_estimate) > 2*ndiff) {
                    outlier_count++;
                    interp = gi_estimate;
                }
                
                //interp = gi_estimate; // override geo interpolant
                
                int proposed = lrint(interp); 
                
                maderr += fabs(double(cvimg.at<uint16_t>(row, col) - proposed));
                interp_count++;
                    
                cvimg.at<uint16_t>(row, col) = proposed;
            }
        }
    }
    
    logger.debug("MAD interpolation error: %lg\n", maderr/double(interp_count));
    logger.debug("Fraction of outliers: %lf\n", outlier_count/double(interp_count));
    for (size_t i=0; i < ncount.size(); i++) {
        logger.debug("%lu neighours = %lf\n", i, ncount[i]/double(interp_count));
    }
    
    imwrite(string("prewhite.png"), rawimg);
    imwrite(string("white.png"), cvimg);
    exit(1);
}
