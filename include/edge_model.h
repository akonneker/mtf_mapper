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
#ifndef EDGE_MODEL_H
#define EDGE_MODEL_H

#include "include/logger.h"
#include <opencv2/imgproc/imgproc.hpp>

#include <map>
using std::multimap;
using std::pair;
using std::make_pair;

#include <cmath>
#include <array>

// Concept: an edge is modelled as a line (point on line and unit direction vector)
// The edge can be refined by populating the coefficients of a quadratic polynomial
// that allows for some curvature in the edge
// The ultimate definition of the edge is a dense set of points fitted locally to the
// gradient of the input image, called the "ridge"

// We can usually obtain estimates of the centroid, direction and quadratic coefficients
// during the quadrangle detection phase
// This information is not available in --single-roi mode, but we could extend that mode
// to produce a rough estimate


class Edge_model {
  public:
    Edge_model(void) {}
    
    Edge_model(const Point2d& centroid, const Point2d& direction) 
    : centroid{centroid}, direction{direction}, normal{-direction.y, direction.x}  {
        coeff[0] = coeff[1] = coeff[2] = 0; 
    }
    
    Edge_model(const Point2d& centroid, const Point2d& direction, const double coeffs[3]) 
    : centroid{centroid}, direction{direction}, normal{-direction.y, direction.x} {
        coeff[0] = coeffs[0];
        coeff[1] = coeffs[1];
        coeff[2] = coeffs[2];
    }

    void hint_point_set_size(int par_dist_bias, int max_edge_len, int est_roi_width) {
        par_bias = par_dist_bias;
        int est_par_length = std::max(max_edge_len, 2 * std::abs(par_bias) + 2);

        points.resize(est_par_length);
        for (auto& v : points) {
            v.reserve(est_roi_width);
        }
        points_hinted = true;        
    }
    
    void add_point(double x, double y, double weight, double distance_threshold=100) {
        Point2d delta = Point2d(x, y) - centroid;
        double perp = delta.ddot(normal);
        double par = delta.ddot(direction);

        double pred_perp = coeff[0]*par*par + coeff[1]*par + coeff[2];
        if (fabs(perp - pred_perp) < distance_threshold) {
            int ipar = rint(par) + par_bias;
            if (ipar >= 0 && ipar < (int)points.size()) {
                points[ipar].push_back(cv::Point3d(x, y, weight));
            }
        }
    }
    
    void estimate_ridge(void) {
        if (!points_hinted) {
            logger.error("%s\n", "did you forget to call hint_point_set_size() before adding points?");
            return;
        }
        
        if (points.size() < 9) {
            logger.debug("too few points in estimate_ridge (only %ld points!)\n", points.size());
            release_points();
            return;
        }
    
        vector<cv::Point3d> samples;
        samples.reserve(points.size());
        for (int ipar=0; ipar < (int)points.size(); ipar++) {
            if (points[ipar].size() > 12) {
                cv::Point3d sample(0, 0, 0);
                for (const auto& p : points[ipar]) {
                    double w = p.z*p.z*p.z*p.z;
                    sample.x += p.x * w;
                    sample.y += p.y * w;
                    sample.z += w;
                }
                samples.push_back(sample);
            }
        }
        
        if (samples.size() < 9) {
            logger.debug("%s\n", "too few points in estimate_ridge");
            release_points();
            return;
        }
        
        // At this point we have an estimate of the ridge curve, but we are almost guaranteed to have some contamination
        // at the edges of the ridge, caused by large-gradient points near the corners
        // Because the connected component contour can differ substantially from the true edge location, especially in the corners
        // we have little choice but to clean up the edges of the ridge after the fact
        
        vector<Point2d> est_ridge;
        est_ridge.reserve(samples.size());
        for (int i=0; i < (int)samples.size(); i++) {
            cv::Point3d sum(0, 0, 0);
            for (int d=-2; d <= 2; d++) {
                if ((i+d) >= 0 && (i+d) < (int)samples.size()) {
                    sum += samples[i+d];
                }
            }
            est_ridge.push_back(Point2d(sum.x/sum.z, sum.y/sum.z));
        }
        
        double mean_delta = pnorm(est_ridge[est_ridge.size()/4] - est_ridge[est_ridge.size()/4-1]);
        double prev_mean_delta = mean_delta;
        double count_delta = 1;
        double var_delta = 0;
        double prev_var_delta = 0;
        for (size_t i=est_ridge.size()/4+1; i < 3*est_ridge.size()/4; i++) {
            count_delta += 1;
            double delta = pnorm(est_ridge[i] - est_ridge[i-1]);
            mean_delta = prev_mean_delta + (delta - prev_mean_delta)/count_delta;
            var_delta = prev_var_delta + (delta - prev_mean_delta) * (delta - mean_delta);
            
            prev_mean_delta = mean_delta;
            prev_var_delta = var_delta;
        }
        var_delta /= count_delta - 1;
        double sdev_delta = sqrt(var_delta);
        
        // trim the edges if necessary
        size_t lower_idx = 1;
        while (lower_idx < est_ridge.size()/4 && 
            fabs(pnorm(est_ridge[lower_idx] - est_ridge[lower_idx-1]) - mean_delta) > 2*sdev_delta) {
            
            lower_idx++;
        }
        
        size_t upper_idx = est_ridge.size() - 1;
        while (upper_idx > 3*est_ridge.size()/4 && 
            fabs(pnorm(est_ridge[upper_idx] - est_ridge[upper_idx-1]) - mean_delta) > 2*sdev_delta) {
            
            upper_idx--;
        }
        
        ridge.reserve(upper_idx - lower_idx + 2);
        for (size_t i=lower_idx; i <= upper_idx; i++) {
            ridge.push_back(est_ridge[i]);
        }

        release_points();
        return;
    }
    
    void update_location(const Point2d& new_centroid, const Point2d& new_direction) {
        coeffs_invalidated = true;
        centroid = new_centroid;
        direction = new_direction;
        normal = Point2d(-direction.y, direction.x);
    }
    
    std::array<double, 3> quad_coeffs(void) { // rather unfortunate, but this method cannot be const
        if (coeffs_invalidated) {
            fit_quad_to_ridge();
            if (coeffs_invalidated) { // if the fit failed, fall back to a pure line model
                coeff = {0,0,0};
            }
        }
        return coeff;
    }
    
    cv::Point3d line_deviation(void) {
        quad_coeffs(); // force update, if necessary
        
        Point2d min_recon(0, 1e20);
        Point2d max_recon(0, -1e20);
        for (size_t i=0; i < ridge.size(); i++) {
            Point2d delta = ridge[i] - centroid;
            double par = delta.ddot(direction);
            
            double pred = coeff[0]*par*par + coeff[1]*par + coeff[2];
            if (pred < min_recon.y) {
                min_recon = Point2d(par, pred);
            }
            if (pred > max_recon.y) {
                max_recon = Point2d(par, pred);
            }
        }
        if (fabs(max_recon.x - min_recon.x) != 0) {
            cv::Point3d deviation_rise_run(0, fabs(max_recon.y - min_recon.y), fabs(max_recon.x - min_recon.x));
            deviation_rise_run.x =  deviation_rise_run.y / deviation_rise_run.z;
            return deviation_rise_run;
        }
        return cv::Point3d(0, 0, 1);
    }
    
    const Point2d& get_centroid(void) const {
        return centroid;
    }
    
    const Point2d& get_direction(void) const {
        return direction;
    }
    
    const Point2d& get_normal(void) const {
        return normal;
    }
    
    bool quad_fit_valid(void) const {
        return !coeffs_invalidated;
    }
   
    vector<Point2d> ridge;

  private:

    void release_points(void) {
        points.clear();
    }
  
    inline double pnorm(const Point2d& dir) const {
        return fabs(normal.ddot(dir));
    }
  
    void fit_quad_to_ridge(void) {
        if (ridge.size() <= 3) {
            coeffs_invalidated = true;
            return; // simply refuse to re-estimate quad
        }
    
        double D[3][3];
        double B[3];
        for (int r=0; r < 3; r++) {
            B[r] = 0;
            for (int c=0; c < 3; c++) {
                D[r][c] = 0;
            }
        }
        
        // scale the x values to improve the condition number
        double span = (ridge.back() - centroid).ddot(direction) - (ridge.front() - centroid).ddot(direction);
        double x_scale = fabs(span)*0.5/sqrt(0.5);
        
        for (size_t i=0; i < ridge.size(); i++) {
            Point2d delta = ridge[i] - centroid;
            double par = delta.ddot(direction) / x_scale;
            double perp = delta.ddot(normal);
            
            D[0][0] += 1;
            D[0][1] += par;
            D[0][2] += par*par;
            D[1][2] += par*par*par;
            D[2][2] += par*par*par*par;
            
            B[0] += perp;
            B[1] += par*perp;
            B[2] += par*par*perp;
        }
        
        // complete the design matrix D, and calculate quadratic coefficients
        D[1][0] = D[0][1];
        D[1][1] = D[0][2];
        D[2][0] = D[0][2];
        D[2][1] = D[1][2];
        
        // what do we do if the determinant is zero?
        
        double det_D =
             D[0][0] * (D[1][1]*D[2][2] - D[1][2]*D[2][1]) +
            -D[0][1] * (D[1][0]*D[2][2] - D[1][2]*D[2][0]) +
             D[0][2] * (D[1][0]*D[2][1] - D[1][1]*D[2][0]);
             
        if (fabs(det_D) < 1e-12) {
            logger.debug("%s\n", "Weird. Quadratic fit failed on ridge fitting");
            coeffs_invalidated = true;
            return;
        }
        
        double det_0 = // replace column 0 with B
                B[0] * (D[1][1]*D[2][2] - D[1][2]*D[2][1]) +
            -D[0][1] * (   B[1]*D[2][2] - D[1][2]*B[2]   ) +
             D[0][2] * (   B[1]*D[2][1] - D[1][1]*B[2]   );
             
        double det_1 = // replace column 1 with B
             D[0][0] * (   B[1]*D[2][2] - D[1][2]*B[2]   ) +
            -   B[0] * (D[1][0]*D[2][2] - D[1][2]*D[2][0]) +
             D[0][2] * (D[1][0]*B[2]    -    B[1]*D[2][0]);
             
        double det_2 = // replace column 2 with B
             D[0][0] * (D[1][1]*B[2]    - B[1]   *D[2][1]) +
            -D[0][1] * (D[1][0]*B[2]    - B[1]   *D[2][0]) +
                B[0] * (D[1][0]*D[2][1] - D[1][1]*D[2][0]);
        
        coeff[0] = det_2/det_D / (x_scale * x_scale);
        coeff[1] = det_1/det_D / x_scale;
        coeff[2] = det_0/det_D;
        
        coeffs_invalidated = false;
    }
  
    Point2d centroid = {0, 0};
    Point2d direction = {0, 0};
    Point2d normal = {0, 0};
    std::array<double, 3> coeff; // least squares quadratic fit
    bool coeffs_invalidated = false;
    
    int par_bias = 0;
    bool points_hinted = false;
    vector<vector<cv::Point3d>> points;
};

#endif
