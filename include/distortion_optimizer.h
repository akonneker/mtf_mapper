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

#ifndef DISTORTION_OPTIMIZER_H
#define DISTORTION_OPTIMIZER_H

#include "Eigen/Dense"
#include <random>
#include <limits>
#include <cmath>

#include <opencv2/imgproc/imgproc.hpp>

#include "include/block.h"

class Ridge {
  public:
    Ridge(const vector<Point2d>& ridge, const Point2d& centroid, const Point2d& normal, double weight=0) 
    : ridge(ridge), centroid(centroid), normal(normal), weight(weight) {}
    
    vector<Point2d> ridge;
    Point2d centroid;
    Point2d normal;
    double weight;
    double residual = 0;
};

class Distortion_optimizer {
  public:
    Distortion_optimizer(const vector<Block>& in_blocks, const Point2d& prin);
    
    Point2d get_max_val(void) const {
        return max_val;
    }
    
    void solve(void);
    
    Point2d inv_warp(const Point2d& p, const Eigen::VectorXd& v);
    double model_not_invertible(const Eigen::VectorXd& v);
    
    double medcouple(vector<float>& x);
    double evaluate(const Eigen::VectorXd& v, double penalty=1.0);
    
    void seed_simplex(Eigen::VectorXd& v, const Eigen::VectorXd& lambda);
    void simplex_sum(Eigen::VectorXd& psum);
    
    
    void nelder_mead(const double ftol, int& num_evals);
    double try_solution(Eigen::VectorXd& psum, const int ihi, const double fac);
    Eigen::VectorXd iterate(double tol);
    
    bool optimization_failure(void) const {
        return nelder_mead_failed;
    }
    
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    vector<Ridge> ridges;
    Point2d prin;
    double radius_norm;
    Eigen::VectorXd best_sol;
    Point2d max_val;
    double maxrad;
    
    // variables used by nelder-mead
    vector<Eigen::VectorXd> np;
    Eigen::VectorXd ny;
    bool nelder_mead_failed;
    Eigen::VectorXd initial;
    double focal_lower;
    double focal_upper;
    double focal_mode_constraint;
};

#endif
