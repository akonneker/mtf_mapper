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

#include <cmath>
#include <vector>
using std::vector;
#include <algorithm>

#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

class Cubic_spline_surface {
  public:
    Cubic_spline_surface(int ncells_x, int ncells_y);
    
    Eigen::MatrixXd spline_fit(const vector<Eigen::Vector3d>& data, double lambda, 
        const Eigen::Vector2d& img_dims, int pruning_threshold=1);

  private:
    inline double B3(double t) {
        if (t < 0) return 0;
        if (t >= 4) return 0;
        if (t < 1) return t*t*t/6.0;
        if (t < 2) return (-3*t*t*t + 12*t*t - 12*t + 4)/6.0;
        if (t < 3) return ( 3*t*t*t - 24*t*t + 60*t - 44)/6.0;
        return (4 - t)*(4 - t)*(4 - t)/6.0;
    }
    
    vector<double> knots_x;
    vector<double> knots_y;
    
    int ncells_x;
    int ncells_y;
    
    double kmax_x;
    double kmax_y;
};




