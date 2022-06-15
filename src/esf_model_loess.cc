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

#include "include/esf_model_loess.h"
#include "include/fastexp.h"
#include <Eigen/Dense>

#include "include/pava.h"

// Note: The kernel function is used during ESF construction, but remember that this kernel
// is not used (at all) when calculating MTF corrections
static inline double local_kernel(double d, double alpha, double w=0.125) {
    if (fabs(d) < w) return 1.0;
    return fastexp(-fabs(fabs(d)-w)*alpha);
}

static double interpolate(double cnr, const vector< std::pair<double, double> >& lut) {
    if (cnr < lut.back().first) {
        if (cnr > lut.front().first) {
            size_t i = lut.size() - 1;
            while (i > 0 && cnr < lut[i].first) {
                i--;
            }
            if (i < lut.size() - 1) {
                return lut[i].second + (lut[i+1].second - lut[i].second) * (cnr - lut[i].first) / (lut[i+1].first - lut[i].first);
            } else {
                return lut.back().second;
            }
        } else {
            return lut.front().second;
        }
    } else {
        return lut.back().second;
    }
}

static void local_moving_average_smoother(vector<double>& smoothed, double* sampled, int fft_size, 
    int fft_left, int fft_right, int left_trans, int right_trans, int width) {
    
    if (width < 1) return;
    width = std::min(width, 32);
    if (left_trans - fft_left <= 1) return;
    if (fft_right - right_trans <= 1) return;
    
    smoothed[0] = sampled[0];
    for (int idx=1; idx < fft_size; idx++) {
        smoothed[idx] = smoothed[idx-1] + sampled[idx];
    }
    
    const int bhw = width;
    left_trans = std::max(left_trans, fft_left + bhw);
    for (int idx=std::max(fft_left, bhw + 1); idx < left_trans; idx++) {
        sampled[idx] = (smoothed[idx+bhw] - smoothed[idx-bhw-1])/double(2*bhw+1);
    }
    double delta = sampled[left_trans] - (smoothed[left_trans+bhw] - smoothed[left_trans-bhw-1])/double(2*bhw+1);
    for (int idx=left_trans; idx < fft_size; idx++) {
        sampled[idx] -= delta;
        smoothed[idx] -= delta;
    }
    right_trans = std::min(right_trans, fft_right - bhw - 1);
    delta = sampled[right_trans] - (smoothed[right_trans+bhw] - smoothed[right_trans-bhw-1])/double(2*bhw+1);
    for (int idx=right_trans; idx < std::min(fft_right, fft_size-bhw-1); idx++) {
        sampled[idx] = (smoothed[idx+bhw] - smoothed[idx-bhw-1])/double(2*bhw+1);
    }
    for (int idx=right_trans; idx < fft_size; idx++) {
        sampled[idx] += delta;
    }
}

static void gauss_smooth(vector<double>& smoothed, double* sampled, size_t start_idx, size_t end_idx, int strength=1, double edge_value=0.5) {
    
    const int sgh = std::min(strength, 8);
    double w[2*8+1] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    const double c = -log(edge_value);
    for (int x=-sgh; x <= sgh; x++) {
        w[x+sgh] = exp(-c * (double(x)/double(sgh)) * (double(x)/double(sgh)));
    }
    for (int idx=start_idx; idx <= int(end_idx); idx++) {
        smoothed[idx] = 0;
        double wsum = 0;
        for (int x=-sgh; x <= sgh; x++) {
            smoothed[idx] += sampled[idx+x] * w[x+sgh];
            wsum += w[x+sgh];
        }
        smoothed[idx] /= wsum;
    }
    
    for (int idx=start_idx; idx <= int(end_idx); idx++) {
        sampled[idx] = smoothed[idx];
    }
}

int Esf_model_loess::build_esf(vector< Ordered_point  >& ordered, double* sampled, 
    const int fft_size, double max_distance_from_edge, vector<double>& esf, Snr& snr,
    bool allow_peak_shift) {
    
    thread_local vector<double> weights(fft_size, 0);
    thread_local vector<double> mean(fft_size, 0);

    for (int i=0; i < fft_size; i++) {
        sampled[i] = missing;
    }
    
    int rval = 0;
    
    int fft_left = 0;
    int fft_right = fft_size-1;
    int twidth = 32;
    
    rval = estimate_esf_clipping(ordered, sampled, fft_size, allow_peak_shift, 
        max_distance_from_edge, mean, weights, fft_left, fft_right, twidth, snr
    );
    
    double cnr = snr.mean_cnr();
    double contrast = snr.contrast();
    
    if (rval < 0) return rval;
    
    fill(weights.begin(), weights.end(), 0.0);
    fill(mean.begin(), mean.end(), 0.0);
    
    auto left_it = ordered.begin();
    auto right_it = ordered.end();
    double alpha = get_alpha();
    constexpr double model_switch_bias = 3*0.125;
    int min_left_bin = ordered.front().first*8 + fft_size/2 + 1;
    int max_right_bin = ordered.back().first*8 + fft_size/2 - 1;
    
    vector< std::pair<double, double> > ridge_lut = {
        {4, 95}, {5, 92}, {6, 89}, {7, 86}, {8, 83}, {9, 79}, {10, 76}, {12, 70}, {14, 63}, {16, 58}, 
        {18, 53}, {20, 49}, {22, 46}, {24, 43}, {26, 42}, {28, 40}, {30, 39}, {35, 38}, {40, 38}, 
        {45, 36}, {50, 34}, {60, 30}, {70, 27}, {80, 22}, {90, 14}, {100, 7}, {124, 2}, {149, 2},
        {794, 1}, {894, 1}, {1000, 5e-8}
    };
    const double ridge_parm = interpolate(cnr, ridge_lut);
    
    constexpr int order = 6;
    for (int b=min_left_bin; b < max_right_bin; b++) {
        double mid = (b - fft_size/2)*0.125;
        
        constexpr double loess_span = 4.5;
        left_it = lower_bound(ordered.begin(), ordered.end(), mid - 0.5*loess_span);
        right_it = lower_bound(left_it, ordered.end(), mid + 0.5*loess_span);
        
        size_t npts = right_it - left_it;
        if (npts < (order+1)) {
            printf("empty interval in bin %d\n", b);
            weights[b] = 0;
        } else {
            if (fabs(mid) < 0.125*twidth + model_switch_bias) {
                Eigen::MatrixXd design(npts, order + 1);
                Eigen::VectorXd v(npts);
                
                double fw = 0.125;
                if (fabs(mid) >= 0.125*0.5*twidth) {
                    fw = 0.250;
                }
                
                size_t row = 0;
                for (auto it=left_it; it != right_it; it++, row++) {
                    double d = fabs(it->first - mid);
                    double w = local_kernel(d, alpha, fw); 
                    
                    double x = (it->first - mid)/(0.5*loess_span);
                    v[row] = w*(it->second - mean[b-1])/contrast;
                    double x2 = x*x;
                    double x3 = x2*x;
                    double x4 = x2*x2;
                    double x5 = x4*x;
                    double x6 = x4*x2;
                    
                    design(row, 0) = w*1;
                    design(row, 1) = w*(x);
                    design(row, 2) = w*(2*x2 - 1);
                    design(row, 3) = w*(4*x3 - 3*x);
                    design(row, 4) = w*(8*x4 - 8*x2 + 1);
                    design(row, 5) = w*(16*x5 - 20*x3 + 5*x);
                    design(row, 6) = w*(32*x6 - 48*x4 + 18*x2 - 1);
                }
                const double phi = (fabs(mid) > 0.75*twidth*0.125) ? ridge_parm : 5e-8;
                Eigen::VectorXd sol = (design.transpose() * design + phi*Eigen::MatrixXd::Identity(order+1, order+1)).llt().solve(design.transpose() * v);
                
                mean[b] = (sol[0] + sol[4])*contrast - (sol[2] + sol[6])*contrast + mean[b-1];
                weights[b] = 1.0;
            } else {
                // blend in the box filter, otherwise f/16 LSFs lose some detail right at the foot
                const double upper_thresh = 2.25*0.125*twidth;
                const double mid_thresh = 2.0*0.125*twidth;
                const double lower_thresh = 0.125*twidth + model_switch_bias;
                constexpr double min_span = 0.24; // such a narrow interval can produce empty bins, handled below
                constexpr double mid_span = 1.0;
                constexpr double max_span = 1.5;
                double span = min_span;
                if (fabs(mid) >= upper_thresh) {
                    span = max_span;
                } else {
                    if (fabs(mid) >= mid_thresh) {
                        double t = (fabs(mid) - mid_thresh)/(upper_thresh - mid_thresh);
                        span = mid_span + t*(max_span - mid_span);
                    } else {
                        double t = (fabs(mid) - lower_thresh)/(mid_thresh - lower_thresh);
                        span = min_span + t*(mid_span - min_span);
                    }
                }
                left_it = lower_bound(left_it, right_it, mid - span);
                right_it = lower_bound(left_it, right_it, mid + span);
                
                double sum = 0;
                for (auto it=left_it; it != right_it; it++) {
                    sum += it->second;
                }
                if (right_it - left_it > 0) {
                    mean[b] = sum / double(right_it - left_it);
                    weights[b] = 1.0;
                } else {
                    mean[b] = mean[b-1];
                    weights[b] = weights[b-1];
                }
            } 
        }
    }

    int left_non_missing = 0;
    int right_non_missing = 0;
    
    // some housekeeping to take care of missing values
    for (int i=0; i < fft_size; i++) {
        sampled[i] = 0;
    }
    for (int idx=fft_left-1; idx <= fft_right+1; idx++) {
        if (weights[idx] > 0) {
            sampled[idx] = mean[idx];
            if (!left_non_missing) {
                left_non_missing = idx; // first non-missing value from left
            }
            right_non_missing = idx; // last non-missing value
        } else {
            sampled[idx] = missing;
        }
    }
    
    // estimate a reasonable value for the non-missing samples
    int nm_target = cnr > 500 ? 8 : 8*3;
    int nm_count = 1;
    double l_nm_mean = sampled[left_non_missing];
    for (int idx=left_non_missing+1; idx < fft_size/2 && nm_count < nm_target; idx++) {
        if (sampled[idx] != missing) {
            l_nm_mean += sampled[idx];
            nm_count++;
        }                                                                                                                                                                   
    }                                                                                                                                                                       
    l_nm_mean /= nm_count;
    
    nm_count = 1;
    double r_nm_mean = sampled[right_non_missing];
    for (int idx=right_non_missing-1; idx > fft_size/2 && nm_count < nm_target; idx--) {
        if (sampled[idx] != missing) {
            r_nm_mean += sampled[idx];
            nm_count++;
        }
    }
    r_nm_mean /= nm_count;
    
    // now just pad out the ends of the sequences with the last non-missing values
    for (int idx=left_non_missing+4; idx >= left_non_missing; idx--) {
        double w = (idx - left_non_missing + 1)/5.0;
        sampled[idx] = w*sampled[idx] + (1.0 - w)*l_nm_mean;
    }
    for (int idx=left_non_missing-1; idx >= 0; idx--) {
        sampled[idx] = l_nm_mean;
    }
    for (int idx=right_non_missing-4; idx <= right_non_missing; idx++) {
        double w = -(idx - right_non_missing - 1)/5.0;
        sampled[idx] = w*sampled[idx] + (1.0 - w)*r_nm_mean;
    }
    for (int idx=right_non_missing+1; idx < fft_size; idx++) {
        sampled[idx] = r_nm_mean;
    }
    
    thread_local vector<double> smoothed(fft_size, 0);
    
    if (apply_monotonic_filter) {
        vector<double> pre_pava(fft_size);
        for (size_t i=0; i < (size_t)fft_size; i++) {
            pre_pava[i] = sampled[i];
        }
        
        // first, reverse the vector if necessary
        bool reversed = false;
        if (l_nm_mean > r_nm_mean) {
            reversed = true;
            std::reverse(sampled, sampled + fft_size);
        }
        jbk_pava(sampled, fft_size);
        if (reversed) {
            std::reverse(sampled, sampled + fft_size);
        }
    }
    
    // smooth out the ESF before the calculate the LSF
    // the degree of smoothing adapts to the estimated CNR
    // the ESF is split into segments, with different types and levels of smoothing applied to different segments
    const vector< std::pair<double, double> > extreme_tails_lut = {
        {4.012, 27}, {5.006, 26}, {6.001, 25}, {6.996, 25}, {7.990, 24}, {8.985, 23}, {9.979, 22}, {11.968, 22}, 
        {13.958, 21}, {15.947, 20}, {17.937, 20}, {19.926, 19}, {21.915, 19}, {23.904, 19}, {25.894, 18}, 
        {27.884, 18}, {29.874, 17}, {34.846, 16}, {39.822, 16}, {44.795, 15}, {49.770, 14}, {59.715, 12}, 
        {69.671, 12}, {79.614, 11}, {89.571, 10}, {99.506, 9}, {124.374, 8}, {149.222, 8}, {174.085, 7}, 
        {198.973, 7}, {298.480, 7}, {397.781, 6}, {497.074, 6}, {595.965, 5}, {695.846, 4}, {794.242, 4}
    };
    int sw_width = lrint(interpolate(cnr, extreme_tails_lut));
    int lt = std::max(fft_left + 4, int(fft_size/2 - 2*twidth));
    int rt = std::min(fft_right - 5, int(fft_size/2 + 2*twidth));
    int prev_lt = lt - 4;
    int prev_rt = rt + 4;
    for (int rep=0; rep < 1; rep++) {
        local_moving_average_smoother(smoothed, sampled, fft_size, 0, fft_size, lt, rt, sw_width);
    }
    
    const vector< std::pair<double, double> > transition_lut = {
        {99.5, 8},{124.4, 8},{149.2, 7},{174.1, 7},{199.0, 6},{298.5, 5},{397.8, 4},{497.1, 3},{596.0, 3},
        {695.8, 3},{794.2, 2},{893.9, 2},{1000,0}, {1100,0}
    };
    sw_width = lrint(interpolate(cnr, transition_lut));
    lt = std::max(fft_left + 2, int(fft_size/2 - 1.25*twidth));
    rt = std::min(fft_right - 3, int(fft_size/2 + 1.25*twidth));
    if (sw_width > 0) {
        gauss_smooth(smoothed, sampled, prev_lt, lt, sw_width, 0.5);
        gauss_smooth(smoothed, sampled, rt, prev_rt, sw_width, 0.5);
    }
    prev_lt = lt - 2;
    prev_rt = rt + 2;
    
    const vector< std::pair<double, double> > core_lut = {
        {4.0, 7},{5.0, 6},{6.0, 6},{7.0, 5},{8.0, 5},{9.0, 5},{10.0, 5},{12.0, 4},{14.0, 4},{15.9, 4},
        {17.9, 3},{19.9, 3},{21.9, 3},{23.9, 3},{25.9, 3},{27.9, 3},{29.9, 3},{34.8, 2},{39.8, 2},{44.8, 2},
        {49.8, 2},{59.7, 2},{69.7, 1},{79.6, 1},{89.6, 1},{99.5, 1},{124.4, 1},{149.2, 1},{174.1, 1},
        {199.0, 1},{298.5, 0},{397.8, 0}
    };
    sw_width = lrint(interpolate(cnr, core_lut));
    lt = std::max(fft_left, int(fft_size/2 - 0.75*twidth));
    rt = std::min(fft_right - 1, int(fft_size/2 + 0.75*twidth));
    if (sw_width > 0) {
        gauss_smooth(smoothed, sampled, prev_lt, lt, sw_width, 0.5);
        gauss_smooth(smoothed, sampled, rt, prev_rt, sw_width, 0.5);
    }
    
    int lidx = 0;
    for (int idx=0; idx < fft_size; idx++) {
        esf[lidx++] = sampled[idx];
    }
    
    double old = sampled[fft_left];
    for (int idx=fft_left; idx <= fft_right; idx++) {
        double temp = sampled[idx];
        sampled[idx] = (sampled[idx+1] - old);
        old = temp;
    }
    for (int idx=0; idx < fft_left+2; idx++) {
        sampled[idx] = 0;
    }
    for (int idx=fft_right-3; idx < fft_size; idx++) {
        sampled[idx] = 0;
    }
    
    return rval;
}

