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
#include "include/cubic_spline_surface.h"
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>

Cubic_spline_surface::Cubic_spline_surface(int ncells_x, int ncells_y)
: ncells_x(ncells_x), ncells_y(ncells_y), 
  kmax_x(ceil(ncells_x/4.0)), kmax_y(ceil(ncells_y/4.0)) {

    for (double y=-4; y < kmax_y + 2; y += 1.0) {
        knots_y.push_back(y);
    }
    
    for (double x=-4; x < kmax_x + 2; x += 1.0) {
        knots_x.push_back(x);
    }

}

size_t find_knot(const vector<double>& knots, double x) {
    int idx = std::upper_bound(knots.begin(), knots.end(), x) - knots.begin();
    return std::max(0, idx - 1);
}

Eigen::MatrixXd Cubic_spline_surface::spline_fit(const vector<Eigen::Vector3d>& in_data, double lambda, 
    const Eigen::Vector2d& img_dims, int pruning_threshold) {
    
    vector<Eigen::Vector4d> data(in_data.size());
    
    Eigen::Vector2d dmin(0, 0);
    Eigen::Vector2d dmax(img_dims[0], img_dims[1]);
    Eigen::MatrixXd occupied(knots_y.size(), knots_x.size());
    occupied.setZero();
    
    cv::Mat occ_count(knots_y.size(), knots_x.size(), CV_32FC1, cv::Scalar::all(0.0));
    cv::Mat occ_val(knots_y.size(), knots_x.size(), CV_32FC1, cv::Scalar::all(0.0));
    
    double min_z = 1e50;
    double max_z = -1e50;
    // rescale input data into knot space
    for (size_t i=0; i < in_data.size(); i++) {
        if ((in_data[i][0] - dmin[0])/(dmax[0] - dmin[0]) >= 1 || 
            (in_data[i][1] - dmin[1])/(dmax[1] - dmin[1]) >= 1) {
            
            printf("bounds exceeded with (%lf, %lf) vs (%lf, %lf), normalized = (%lf, %lf)\n",
                in_data[i][0], in_data[i][1],
                dmax[0], dmax[1],
                (in_data[i][0] - dmin[0])/(dmax[0] - dmin[0]), (in_data[i][1] - dmin[1])/(dmax[1] - dmin[1])
            );
            exit(-1);
        }
        data[i][0] = kmax_x * (in_data[i][0] - dmin[0])/(dmax[0] - dmin[0]);
        data[i][1] = kmax_y * (in_data[i][1] - dmin[1])/(dmax[1] - dmin[1]);
        data[i][2] = in_data[i][2];
        data[i][3] = 1.0;
        min_z = std::min(min_z, data[i][2]);
        max_z = std::max(max_z, data[i][2]);
        
        int col = find_knot(knots_x, data[i][0]); 
        int row = find_knot(knots_y, data[i][1]);
        occupied(row, col) = 1;
        occ_count.at<float>(row, col) += 1;
        occ_val.at<float>(row, col) += data[i][2];
    }
    
    size_t occ_fr_count = 0;
    for (int r=0; r < occ_val.rows; r++) {
        for (int c=0; c < occ_val.cols; c++) {
            if (occ_count.at<float>(r, c) > 0) {
                occ_val.at<float>(r, c) /= occ_count.at<float>(r, c);
                occ_fr_count++;
            }
        }
    }
    double occ_fraction = occ_fr_count / (kmax_x * kmax_y);
    
    cv::Mat first_val = occ_val.clone();
    cv::Mat occ_dist(knots_y.size(), knots_x.size(), CV_32FC1, cv::Scalar::all(0.0));
    
    for (int i=0; i < std::max(occ_val.rows, occ_val.cols); i++) {
        // hard-coded dilation with 3-pixel cross, but skipping pixels with zero counts
        cv::Mat dilated = cv::Mat(occ_val.rows, occ_val.cols, CV_32FC1, cv::Scalar::all(0.0));
        for (int r=0; r < occ_val.rows; r++) {
            for (int c=0; c < occ_val.cols; c++) {
                
                float local_max = -1e10;
                int delta[5][2] = { {0, 0}, {-1, 0}, {1, 0}, {0, 1}, {0, -1} };
                for (int d=0; d < 5; d++) {
                    int dr = r + delta[d][0];
                    int dc = c + delta[d][1];
                    
                    if (dr >= 0 && dc >= 0 && dr < occ_val.rows && dc < occ_val.cols && occ_count.at<float>(dr, dc) > 0) {
                        local_max = std::max(occ_val.at<float>(dr, dc), local_max);
                    }
                }
                if (local_max > -1e10) {
                    dilated.at<float>(r, c) = local_max;
                }
            }
        }
        
        for (int r=0; r < occ_val.rows; r++) {
            for (int c=0; c < occ_val.cols; c++) {
                if (occ_count.at<float>(r, c) == 0 && occ_val.at<float>(r, c) != dilated.at<float>(r, c)) {
                    occ_count.at<float>(r, c) = 1;
                    occ_dist.at<float>(r, c) = i + 1;
                    occ_val.at<float>(r, c) = dilated.at<float>(r, c);
                    first_val.at<float>(r, c) = dilated.at<float>(r, c);
                }
            }
        }
        occ_val = dilated.clone();
    }
    
    constexpr double small_wt = 1e-4;
    for (int r=0; r < occupied.rows(); r++) {
        for (int c=0; c < occupied.cols(); c++) {
            double fill_value = first_val.at<float>(r, c);
            if (occupied(r, c) == 0) {
                data.push_back(Eigen::Vector4d(knots_x[c] + 0.5, knots_y[r] + 0.5, fill_value, small_wt));
            }
        }
    }
    
    size_t L = knots_x.size();
    size_t K = knots_y.size();
    size_t m = data.size();
    Eigen::MatrixXd Bx(m, L);
    Eigen::MatrixXd By(m, K);
    
    // TODO: Actually, we do not need Bx and By at all, and we can probably do an in-place construction directly into C ...
    Eigen::MatrixXd C;
    {
        for (size_t col=0; col < L; col++) {
            for (size_t row=0; row < m; row++) {
                Bx(row, col) = B3(data[row][0] - knots_x[col]);
            }
        }
        for (size_t col=0; col < K; col++) {
            for (size_t row=0; row < m; row++) {
                By(row, col) = B3(data[row][1] - knots_y[col]);
            }
        }
        
        C = Eigen::kroneckerProduct(Bx, Eigen::MatrixXd::Ones(1, K)).cwiseProduct(
            Eigen::kroneckerProduct(Eigen::MatrixXd::Ones(1, L), By)
        );
    }
    
    Eigen::MatrixXd P;
    {
        Eigen::MatrixXd Px;
        Eigen::MatrixXd Py;
        {   
            Eigen::MatrixXd Dx = Eigen::MatrixXd::Zero(L-1, L);
            for (size_t col=0; col < L-1; col++) {
                double scale = knots_x[col] >= -1 ? 1.0 : 0.5;
                Dx(col, col) = -scale;
                Dx(col, col + 1) = scale;
            }
            Px = Dx.transpose()*Dx;
        }
        
        {
            Eigen::MatrixXd Dy = Eigen::MatrixXd::Zero(K-1, K);
            for (size_t col=0; col < K-1; col++) {
                double scale = knots_y[col] >= -1 ? 1.0 : 0.5;
                Dy(col, col) = -scale;
                Dy(col, col + 1) = scale;
            }
            Py = Dy.transpose()*Dy;
        }
        
        P = Eigen::kroneckerProduct(Eigen::MatrixXd::Identity(K, K), Px) + 
            Eigen::kroneckerProduct(Eigen::MatrixXd::Identity(L, L), Py);
    }
    
    Eigen::VectorXd rhs(m);
    Eigen::VectorXd W(m);
    for (size_t i=0; i < data.size(); i++) {
        rhs[i] = data[i][2];
        W[i] = data[i][3];
    }
    
    // if the chart appears a bit sparse, we increase the amount of smoothing a bit
    if (occ_fraction < 0.6) {
        lambda *= 10;
    }
    
    Eigen::VectorXd sol = (C.transpose()*W.asDiagonal()*C + lambda*P).colPivHouseholderQr().solve(C.transpose()*W.asDiagonal()*rhs);
    
    Eigen::MatrixXd grid(ncells_y, ncells_x);
    grid.setZero();
    for (size_t gc=0; gc < (size_t)grid.cols(); gc++) {
        double px = double(gc) * kmax_x / (grid.cols());
        for (size_t gr=0; gr < (size_t)grid.rows(); gr++) {
            double py = double(gr) * kmax_y / (grid.rows());
            double sum = 0;
            
            size_t sol_idx = 0;
            for (size_t row=0; row < L; row++) {
                for (size_t col=0; col < K; col++) {
                    sum += sol[sol_idx++] * B3(px - knots_x[row]) * B3(py - knots_y[col]);
                }
            }
            
            int kx = find_knot(knots_x, px);
            int ky = find_knot(knots_y, py);
            if (occ_dist.at<float>(ky, kx) > pruning_threshold) {
                sum = 0;
            }
            grid(gr, gc) = sum;
        }
    }
    
    return grid;
}