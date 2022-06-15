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
#include "include/grid_interpolator.h"

#include <Eigen/Dense>
#include "include/cubic_spline_surface.h"


void interpolate_grid(const Grid_functor& ftor, Edge_type target_edge_type, cv::Mat& grid_coarse, cv::Mat& grid_fine,
    cv::Size img_dims, const vector<Block>& blocks, double upper, double smoothing_factor, int pruning_threshold) {
    
    grid_fine = ftor.nodata();
    
    Point2d centr(0,0);
    for (size_t i=0; i < blocks.size(); i++) {
        centr += blocks[i].get_centroid();
    }
    centr = Point2d(img_dims.width/2, img_dims.height/2);
    
    for (size_t i=0; i < blocks.size(); i++) {
        for (size_t k=0; k < 4; k++) {
            double val = ftor.value(blocks[i], k);
            if (ftor.in_range(val)) {
                Point2d cent = blocks[i].get_edge_centroid(k);

                Point2d dir = cent - centr;
                dir = dir * (1.0/norm(dir));

                Point2d norm = blocks[i].get_normal(k);
                double delta = dir.x*norm.x + dir.y*norm.y;

                Edge_type edge_type = NEITHER;
                if (target_edge_type == MERIDIONAL) {
                    if (fabs(delta) > cos(45.0/180*M_PI)) {
                        edge_type = MERIDIONAL;
                    } else {
                        edge_type = SAGITTAL;
                    }
                } else if (target_edge_type == SAGITTAL) {
                    if (fabs(delta) > cos(45.0/180*M_PI)) {
                        edge_type = MERIDIONAL;
                    } else {
                        edge_type = SAGITTAL;
                    }
                }
                
                int ly = lrint(cent.y*(grid_fine.rows-1)/img_dims.height); 
                int lx = lrint(cent.x*(grid_fine.cols-1)/img_dims.width);

                if ((fabs(grid_fine.at<float>(ly, ly) - ftor.nodata()) < 1e-6 || fabs(val) > fabs(grid_fine.at<float>(ly,lx))) 
                    && edge_type == target_edge_type) { // max composite
                    grid_fine.at<float>(ly,lx) = val;
                }
            }
            
        }
    }
    
    vector< Eigen::VectorXd > solutions(grid_coarse.rows * grid_coarse.cols, Eigen::VectorXd::Zero(1));
    
    double global_max = -1e50;
    double global_min = 1e50;
    vector<Eigen::Vector3d> samples;
    vector<double> all_vals;
    for (int lr=0; lr < grid_fine.rows; lr++) {
        for (int lc=0; lc < grid_fine.cols; lc++) {
            double val = grid_fine.at<float>(lr, lc);
            if (ftor.in_range(val)) {
                samples.push_back(Eigen::Vector3d(lc, lr, val));
                all_vals.push_back(val);
                global_max = std::max(global_max, val);
                global_min = std::min(global_min, val);
            }
        }
    }
    
    sort(all_vals.begin(), all_vals.end());
    double m_upper = all_vals[97*all_vals.size()/100];
    double m_lower = all_vals[5*all_vals.size()/100];
    m_upper = upper < -100 ? global_max : upper;
    m_lower = global_min;
    
    Cubic_spline_surface csurf(grid_coarse.cols, grid_coarse.rows);
    Eigen::MatrixXd cg = csurf.spline_fit(samples, 
        smoothing_factor,
        Eigen::Vector2d(grid_fine.cols, grid_fine.rows),
        pruning_threshold
    ); 
    
    for (int r=0; r < grid_coarse.rows; r++) {
        for (int c=0; c < grid_coarse.cols; c++) {
            grid_coarse.at<float>(r,c) = ftor.clamp(cg(r, c), m_upper, m_lower);
        }
    }
    
    // use bilinear interpolation to produce a smoother output image
    for (int row=0; row < grid_fine.rows; row++) {
        for (int col=0; col < grid_fine.cols; col++) {
            double y = row*grid_coarse.rows / double(grid_fine.rows);
            double x = col*grid_coarse.cols / double(grid_fine.cols);
            
            int fx = (int)floor(x);
            int fy = (int)floor(y);

            int cx = fx + 1;
            int cy = fy + 1;

            cx = min(cx, grid_coarse.cols - 1);
            cy = min(cy, grid_coarse.rows - 1);

            double xfac = x - fx;
            double yfac = y - fy;

            double w1 = (1 - xfac) * (1 - yfac);
            double w2 =      xfac  * (1 - yfac);
            double w3 = (1 - xfac) *      yfac;
            double w4 =      xfac  *      yfac;
            
            
            double val =
                w1 * grid_coarse.at<float>(fy, fx) + w2 * grid_coarse.at<float>(fy, cx) +
                w3 * grid_coarse.at<float>(cy, fx) + w4 * grid_coarse.at<float>(cy, cx);
            
            grid_fine.at<float>(row, col) = ftor.clamp(val, m_upper, m_lower);
        }
    }
}