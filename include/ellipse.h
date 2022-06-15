#ifndef ELLIPSE_H
#define ELLIPSE_H

#include "include/common_types.h"
#include "include/gradient.h"

#include <Eigen/Dense>

#include <cmath>
#include <string>
using std::string;

#include <set>
using std::set;

const double ellipse_fit_threshold = 0.2;
const double ellipse_fit_deviation_threshold = 0.04;

const double min_major_axis = 26/2.0;
const double min_minor_axis = min_major_axis/3.0; // we expect a 50% reduction on a 45 degree slope
const double max_ellipse_gradient_error = 15;

using std::pair;
using std::make_pair;
typedef std::pair<int, int> iPoint;

class Component_labeller;

class Ellipse_detector {
  public:
    Ellipse_detector(void)
      : centroid_x(0), centroid_y(0), major_axis(1), minor_axis(1), angle(0),
        quality(0), fg_fraction(0), valid(false), solid(false), code(-1), _C(3,3) {
                          
        _C.setZero();
    }

    ~Ellipse_detector(void) {}
    
    static Point2d calculate_centroid(Pointlist& p);
    
    int fit(const Component_labeller& cl, const Gradient& gradient,
        const Pointlist& raw_points, int tl_x, int tl_y, int dilate=1);
    
    bool gradient_check(const Component_labeller& cl, const Gradient& gradient, const Pointlist& raw_points);
    
    void set_code(int incode) {
        code = incode;
    }
    
    int _matrix_to_ellipse(Eigen::Matrix3d& C);
    void _dilate(set<iPoint>& s, int width, int height, int iters);
    void _dilate_outer_only(set<iPoint>& s, int width, int height);
    void _correct_eccentricity(double major_scale, double bp_x, double bp_y);

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    double centroid_x;
    double centroid_y;
    double major_axis;
    double minor_axis;
    double angle;
    double quality;
    double fg_fraction;
    bool valid; // flagged by detector at a later stage
    bool solid;
    int code;
    
    map<int, scanline> scanset;
    
    Eigen::Vector3d pos1;
    Eigen::Vector3d pos2;
    
    Eigen::Vector3d n1;
    Eigen::Vector3d n2;

  private:
    Eigen::Matrix3d _C;
};

#endif
