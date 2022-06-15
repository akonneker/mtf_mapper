#include "include/logger.h"
#include "include/ellipse.h"
#include "include/component_labelling.h"
#include "include/point_helpers.h"

#include "opencv2/imgproc/imgproc.hpp" // required for fitEllipse

#include <cmath>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
using namespace Eigen;

#include <limits>

typedef Matrix<double, 5, 5> Matrix5d;
typedef Matrix<double, 5, 1> Vector5d;


int Ellipse_detector::fit(const Component_labeller& cl, const Gradient& gradient,
    const Pointlist& raw_points, int tl_x, int tl_y, int dilate) {
    
    if (raw_points.size() < 2) { // don't even bother
        return 0;
    }
    
    int width  = gradient.width();
    int height = gradient.height();
    
    const cv::Mat& grad_x = gradient.grad_x();
    const cv::Mat& grad_y = gradient.grad_y();
    //const cv::Mat& grad_m = gradient.grad_magnitude();

    set<iPoint> boundary;

    const int border = 1;
    bool edge_touched = false;
    double mx = 0;
    double my = 0;
    double wsum = 0;
    double cov_xx = 0;
    double cov_yy = 0;
    double cov_xy = 0;
    double circ = 0; // circumference of curve
    cv::Point2d prev_point = raw_points.back();
    for (size_t i=0; i < raw_points.size(); i++) {
    
        double temp = wsum + 1.0;
        double delta_x = raw_points[i].x - mx;
        double delta_y = raw_points[i].y - my;
        double rx = delta_x / temp;
        double ry = delta_y / temp;
        mx += rx;
        my += ry;
        
        cov_xx += wsum * delta_x * rx;
        cov_yy += wsum * delta_y * ry;
        cov_xy += wsum * delta_x * ry;
        
        wsum = temp;
    
        boundary.insert(iPoint(lrint(raw_points[i].x), lrint(raw_points[i].y)) );
        
        if (raw_points[i].x <= border || raw_points[i].x >= width - 1 - border ||
            raw_points[i].y <= border || raw_points[i].y >= height - 1 - border) {
            
            edge_touched = true;
        }
        
        
        circ += sqrt(SQR(raw_points[i].x - prev_point.x) + SQR(raw_points[i].y - prev_point.y));
        prev_point = raw_points[i];
    }
    cov_xx /= wsum - 1;
    cov_xy /= wsum - 1;
    cov_yy /= wsum - 1;
    
    if (edge_touched) return 0; // skip objects that touch the edge of the image (ellipses not really allowed there)

    if (std::isnan(mx) || std::isnan(my)) {
        return 0; // 0 -> not a circle
    }
    
    // use eigenvalues of covariance matrix to estimate ellipse
    double tr = cov_xx + cov_yy;
    double det = cov_xx*cov_yy - cov_xy*cov_xy;
    
    double pa=1.0;
    double pb=-tr;
    double pc=det;
    
    double q = -0.5 * (pb + sgn(pb)*sqrt(pb*pb - 4*pa*pc) );
    double l1 = q/pa;
    double l2 = pc / q;
    
    double i_maj = std::max(sqrt(2*l1), sqrt(2*l2));
    double i_min = std::min(sqrt(2*l1), sqrt(2*l2));
    
    double h = SQR(i_maj - i_min)/SQR(i_maj + i_min);
    double ell_circ = M_PI*(i_maj + i_min)*(1 + 3*h/(10 + sqrt(4 - 3*h))); // Ramanujan's approximation
    
    // early exit when the curve circumference is too different from the expected
    // circumference of a best-fit ellipse estimated using PCA
    if (fabs(circ / ell_circ - 1) > 0.25) { // this might be a fairly tight bound ... maybe it should depend on the size of the ellipse?
        return 0;
    }

    _dilate(boundary, width, height, dilate);
    _dilate_outer_only(boundary, width, height);

    Matrix<double, 5, 5> wK;
    Matrix<double, 5, 1> K; 
    K.setZero();
    Vector3d l;
    l.setZero();
    Matrix<double, 5, 1> rhs;
    rhs.setZero();

    double sum_c4 = 0;

    wK.setZero();
    int count = (int)boundary.size();
    Matrix<double, Eigen::Dynamic, 3> L(count, 3);
    L.setZero();

    double mean_dist = 0;
    int counter = 0;
    for (set<iPoint>::const_iterator it=boundary.begin(); it != boundary.end(); it++) {
        const int& x_pos = it->first;
        const int& y_pos = it->second;

        if (gradient.grad_magnitude(x_pos, y_pos) > 1e-7) {
            Vector3d v;
            v << x_pos,y_pos,1;
            double dist = sqrt( SQR(x_pos-mx) + SQR(y_pos-my) );
            mean_dist += dist;
            counter++;
        }
    }
    mean_dist /= counter;
    double iso_scale = sqrt(2.0) / mean_dist;

    size_t idx = 0;
    for (set<iPoint>::const_iterator it=boundary.begin(); it != boundary.end(); it++) {
        const int& x_pos = it->first;
        const int& y_pos = it->second;
        
        if (gradient.grad_magnitude(x_pos, y_pos) > 1e-7) {


            l[0] = grad_x.at<float>(y_pos, x_pos);
            l[1] = grad_y.at<float>(y_pos, x_pos);
            l[2] = -(grad_x.at<float>(y_pos, x_pos) * (x_pos-mx) * iso_scale + grad_y.at<float>(y_pos, x_pos) * (y_pos-my) * iso_scale);

            L.row(idx++) = l;

        }
    }

    for (size_t r=0; r < (size_t)L.rows(); r++) {

        Vector3d l = L.row(r);

        sum_c4 += SQR(SQR(l[2]));

        K[0] = l[0]*l[0];
        K[1] = l[0]*l[1];
        K[2] = l[1]*l[1];
        K[3] = l[0]*l[2];
        K[4] = l[1]*l[2];

        double weight = (l[0]*l[0] + l[1]*l[1]);

        wK += weight * K * K.transpose();

        rhs -= weight * K*(l[2]*l[2]);
    }
    
    if (rhs.rows() != 5 || rhs.cols() != 1) {
        logger.debug("%s\n", "rhs undefined");
        exit(1);
    }
    
    for (size_t ri=0; ri < size_t(rhs.rows()); ri++) {
        if (std::isnan(rhs(ri,0))) {
            return 0;
        }
    }

    Vector5d sol;
    sol = wK.fullPivHouseholderQr().solve(rhs);

    Matrix3d Cstar;
    Cstar(0, 0) = sol[0];
    Cstar(0, 1) = sol[1]*0.5;
    Cstar(0, 2) = sol[3]*0.5;
    Cstar(1, 0) = sol[1]*0.5;
    Cstar(1, 1) = sol[2];
    Cstar(1, 2) = sol[4]*0.5;
    Cstar(2, 0) = sol[3]*0.5;
    Cstar(2, 1) = sol[4]*0.5;
    Cstar(2, 2) = 1;      // F* == 1


    Matrix3d C = Cstar.inverse();
    C *= 1.0/C(2,2);
    _C = C;

    Matrix<double, 1, 5> s;
    s.row(0) = sol;
    double sAAs = (s * wK * (s.transpose()))(0,0);
    double R = (sAAs - 2*(sol.dot(rhs)) + sum_c4) / double(count - 5);
    Matrix<double, 5, 5> cov = wK.inverse();
    cov = cov * R;
    Matrix2d cov2;
    cov2(0, 0) = cov(3, 3);
    cov2(0, 1) = cov(3, 4);
    cov2(1, 0) = cov(4, 3);
    cov2(1, 1) = cov(4, 4);
    JacobiSVD<Matrix2d> svd(cov2, ComputeFullU | ComputeFullV);
    Matrix2d Vs = svd.matrixV();
    Vector2d ws = svd.singularValues();

    Matrix2d V;
    V.row(0) = Vs.row(1);
    V.row(1) = Vs.row(0);

    Vector2d w;
    w[0] = ws[1];
    w[1] = ws[0];

    Vector2d centre_uncertainty;
    centre_uncertainty[0] = sqrt(w[0]) * 0.25;
    centre_uncertainty[1] = sqrt(w[1]) * 0.25;

    // shift the ellipse back to the original pixel coordinate system
    Matrix3d S;
    S.setIdentity();
    S(0, 0) = iso_scale;
    S(1, 1) = iso_scale;
    C = _C = S.transpose()*C*S;
    S.setIdentity();
    S(0, 2) = -(tl_x + mx);
    S(1, 2) = -(tl_y + my);
    C = _C = S.transpose()*C*S;

    int result = _matrix_to_ellipse(C);

    if (minor_axis > major_axis) {
        _C = -_C;
        _matrix_to_ellipse(C);
    }

    bool gradient_ok = gradient_check(cl, gradient, raw_points);

    int is_circle = 0;

    if (result == 0 &&
        major_axis >= min_major_axis &&
        minor_axis >= min_minor_axis &&
        gradient_ok) {

        is_circle = 1;
    }


    if (std::isnan(centroid_x) || std::isnan(centroid_y) || std::isnan(major_axis) || std::isnan(minor_axis)) {
        is_circle = 0;
    }
    
    logger.debug("centre (%.2lf, %.2lf), major = %lf, minor = %lf, angle = %lf, is_circle = %d\n",
        centroid_x, centroid_y, major_axis, minor_axis, angle/M_PI*180.0, is_circle);
        
    if (is_circle) {
        scanset.clear();
        for (size_t i=0; i < raw_points.size(); i++) {
            int ix = lrint(raw_points[i].x);
            int iy = lrint(raw_points[i].y);
            
            map<int, scanline>::iterator it = scanset.find(iy);
            if (it == scanset.end()) {
                scanline sl(ix,ix);
                scanset.insert(make_pair(iy, sl));
            }
            if (ix < scanset[iy].start) {
                scanset[iy].start = ix;
            }
            if (ix > scanset[iy].end) {
                scanset[iy].end = ix;
            }
        }
        int clabel = cl(lrint(raw_points[0].x), lrint(raw_points[0].y));
        logger.debug("label used for scanset: %d\n", clabel);
        int total = 0;
        int foreground = 0;
        for (map<int, scanline>::iterator it=scanset.begin(); it != scanset.end(); it++) {
            int y=it->first;
            for (int x=it->second.start; x <= it->second.end; x++) {
                if (cl(x,y) == clabel) {
                    foreground++;
                }
                total++;
            }
        }
        fg_fraction = double(foreground)/double(total);
    }        
    
    return is_circle;
}


void Ellipse_detector::_dilate(set<iPoint>& s, int width, int height, int iters) {
    const int border = 1;

    if (iters > 0) {

        set<iPoint> gen_s;

        for (int k=0; k < iters; k++) {
            for (set<iPoint>::const_iterator it=s.begin(); it != s.end(); it++) {

                gen_s.insert(*it);

                int left = max(border, it->first-1);
                int right = min(width-1-border, it->first+1);
                int top = max(border, it->second-1);
                int bottom = min(height-1-border, it->second+1);

                gen_s.insert(iPoint(right, it->second));
                gen_s.insert(iPoint(right, bottom));
                gen_s.insert(iPoint(it->first, bottom));
                gen_s.insert(iPoint(left, bottom));
                gen_s.insert(iPoint(left, it->second));
                gen_s.insert(iPoint(left, top));
                gen_s.insert(iPoint(it->first, top));
                gen_s.insert(iPoint(right, top));

            }
            s = gen_s; // copy it back
        }
    }
}

void Ellipse_detector::_dilate_outer_only(set<iPoint>& s, int width, int height) {
    const int border = 1;
    
    double cx = 0;
    double cy = 0;
    for (set<iPoint>::const_iterator it=s.begin(); it != s.end(); it++) {
        cx += it->first;
        cy += it->second;
    }
    cx /= s.size();
    cy /= s.size();

    set<iPoint> gen_s;

    for (set<iPoint>::const_iterator it=s.begin(); it != s.end(); it++) {

        gen_s.insert(*it);
        
        Point2d dir(it->first - cx, it->second - cy); // current radial direction

        int left = max(border, it->first-1);
        int right = min(width-1-border, it->first+1);
        int top = max(border, it->second-1);
        int bottom = min(height-1-border, it->second+1);

        if ((right - it->first)*dir.x >= 0) gen_s.insert(iPoint(right, it->second));
        if ((right - it->first)*dir.x >= 0 && (bottom - it->second)*dir.y >= 0) gen_s.insert(iPoint(right, bottom));
        if ((bottom - it->second)*dir.y >= 0) gen_s.insert(iPoint(it->first, bottom));
        if ((left - it->first)*dir.x >= 0 && (bottom - it->second)*dir.y >= 0) gen_s.insert(iPoint(left, bottom));
        if ((left - it->first)*dir.x >= 0) gen_s.insert(iPoint(left, it->second));
        if ((left - it->first)*dir.x >= 0 && (top - it->second)*dir.y >= 0) gen_s.insert(iPoint(left, top));
        if ((top - it->second)*dir.y >= 0) gen_s.insert(iPoint(it->first, top));
        if ((right - it->first)*dir.x >= 0 && (top - it->second)*dir.y >= 0) gen_s.insert(iPoint(right, top));

    }
    s = gen_s; // copy it back
}

int Ellipse_detector::_matrix_to_ellipse(Matrix3d& C) {

    double a = C(0,0);
    double b = C(0,1)*2;
    double c = C(1,1);  
    double d = C(0,2)*2;
    double e = C(1,2)*2;
    double f = C(2,2);  


    double thetarad = 0.5*atan2(b, a - c);
    double cost = cos(thetarad);
    double sint = sin(thetarad);
    double sin_squared = sint*sint;
    double cos_squared = cost*cost;
    double cos_sin = sint*cost;

    double Ao = f;
    double Au =   d * cost + e * sint;
    double Av = - d * sint + e * cost;
    double Auu = a * cos_squared + c * sin_squared + b * cos_sin;
    double Avv = a * sin_squared + c * cos_squared - b * cos_sin;

    if(Auu==0 || Avv==0) {
        // invalid ellipse
        return -1;
    }

    double tuCentre = - Au/(2*Auu);
    double tvCentre = - Av/(2*Avv);
    double wCentre = Ao - Auu*tuCentre*tuCentre - Avv*tvCentre*tvCentre;

    double uCentre = tuCentre * cost - tvCentre * sint;
    double vCentre = tuCentre * sint + tvCentre * cost;

    double Ru = -wCentre/Auu;
    double Rv = -wCentre/Avv;

    Ru = sqrt(fabs(Ru))*(Ru < 0 ? -1 : 1);
    Rv = sqrt(fabs(Rv))*(Rv < 0 ? -1 : 1);

    centroid_x = uCentre;
    centroid_y = vCentre;
    major_axis = Ru;
    minor_axis = Rv;
    angle = thetarad;

    return 0;
}

bool Ellipse_detector::gradient_check([[maybe_unused]] const Component_labeller& cl, const Gradient& gradient, const Pointlist& raw_points) {
    
    // gradient just on ellipse perimeter must be perpendicular to ellipse tangent
    double cosa = cos(-angle);
    double sina = sin(-angle);
    
    vector<double> phi_diff;
    for (size_t i=0; i < raw_points.size(); i++) {
        // now generate a point just outside the ellipse ....
        int px = lrint(raw_points[i].x);
        int py = lrint(raw_points[i].y);
        
        Point2d d(raw_points[i].x - centroid_x, raw_points[i].y - centroid_y);
        double rx = cosa*d.x - sina*d.y;
        double ry = sina*d.x + cosa*d.y;
        double theta = atan2(ry, rx); // ellipse curve parameter theta
        Point2d tangent = normalize(Point2d(-major_axis*sin(theta), minor_axis*cos(theta))); 
        // rotate the tangent vector back to image coordinates
        rx = cosa*tangent.x + sina*tangent.y;
        ry = (-sina)*tangent.x + cosa*tangent.y;
        tangent.x = rx; 
        tangent.y = ry; 
        
        Point2d grad = normalize(Point2d(gradient.grad_x().at<float>(py, px), gradient.grad_y().at<float>(py, px)));
        
        double dot = tangent.x*grad.x + tangent.y*grad.y;
        double phi = acos(dot);
        
        phi_diff.push_back(phi/M_PI*180 - 90);
    }
    
    sort(phi_diff.begin(), phi_diff.end());
    const double phi_percentile = 0.9;
    
    double phi_delta = phi_diff[phi_percentile*phi_diff.size()];
    
    return (phi_delta < max_ellipse_gradient_error) && phi_delta >= 0;
}


