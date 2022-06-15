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

#include "include/esf_sampler_quad.h"

void Esf_sampler_quad::quad_tangency(const Point2d& p, const std::array<double, 3>& qp, vector<double>& roots) {
    const double& a = qp[0];
    const double& b = qp[1];
    const double& c = qp[2];
    
    double ca = 2*a*a;
    double cb = 3*a*b;
    double cc = 1 + 2*a*c - 2*a*p.y + b*b;
    double cd = b*c - p.y*b - p.x;
    
    // convert to normalized cubic with cubic coefficient = 1.0
    double pa = cb / ca;
    double pb = cc / ca;
    double pc = cd / ca;
    
    double Q = (pa*pa - 3*pb)/9.0;
    double R = (2*pa*pa*pa - 9*pa*pb + 27*pc)/54.0;
    
    if (R*R < Q*Q*Q) { // three real roots
        double theta = acos(R/sqrt(Q*Q*Q));
        roots.push_back(-2*sqrt(Q)*cos(theta/3.0) - pa/3.0);
        roots.push_back(-2*sqrt(Q)*cos((theta+2*M_PI)/3.0) - pa/3.0);
        roots.push_back(-2*sqrt(Q)*cos((theta-2*M_PI)/3.0) - pa/3.0);
    } else { // one real root, two complex roots
        double A = -std::copysign(std::cbrt(fabs(R) + sqrt(R*R - Q*Q*Q)), R);
        double B = A == 0 ? 0 : Q/A;
        roots.push_back(A+B - pa/3.0);
    }
}

void Esf_sampler_quad::sample(Edge_model& edge_model, vector<Ordered_point>& local_ordered, 
    const map<int, scanline>& scanset, double& edge_length,
    const cv::Mat& geom_img, const cv::Mat& sampling_img,
    Bayer::cfa_mask_t cfa_mask) {
    
    cfa_mask = cfa_mask == Bayer::DEFAULT ? default_cfa_mask : cfa_mask;
    
    double max_along_edge = -1e50;
    double min_along_edge = 1e50;
    
    const std::array<double, 3>& qp = edge_model.quad_coeffs();
    vector<double> roots;
    roots.reserve(3);
    for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); ++it) {
        int y = it->first;
        if (y < border_width || y > geom_img.rows-1-border_width) continue;
        int rowcode = (y & 1) << 1;
        
        for (int x=it->second.start; x <= it->second.end; ++x) {
            
            if (x < border_width || x > geom_img.cols-1-border_width) continue;
            
            int code = 1 << ( (rowcode | (x & 1)) ^ 3 );
            if ((code & cfa_mask) == 0) continue;
            
            Point2d d = Point2d(x, y) - edge_model.get_centroid();
            double perp = d.ddot(edge_model.get_normal()); 
            double par = d.ddot(edge_model.get_direction());
            
            // (par, perp) is in the local coordinate frame of the parabola qp
            // so find the closest point on qp
            roots.clear();
            quad_tangency(Point2d(par, perp), qp, roots);
            perp = 1e20;
            for (auto r: roots) {
                Point2d on_qp(r, r*r*qp[0] + r*qp[1] + qp[2]);
                Point2d recon = on_qp.x*edge_model.get_direction() + on_qp.y*edge_model.get_normal() + edge_model.get_centroid();
                Point2d delta = Point2d(x, y) - recon;
                double dist = std::copysign(norm(delta), delta.ddot(edge_model.get_normal()));
                if (fabs(dist) < fabs(perp)) {
                    perp = dist;
                }
            }
            
            if (fabs(perp) < max_dot) {
                local_ordered.push_back(Ordered_point(perp, sampling_img.at<uint16_t>(y,x) ));
                max_along_edge = max(max_along_edge, par);
                min_along_edge = min(min_along_edge, par);
                
                /*
                if ((x&1) || (y&1)) {
                    cv::Vec3b& color = od_img.at<cv::Vec3b>(y, x);
                    color[0] = 255;
                    color[1] = 0;
                    color[2] = 0;
                }
                */
            }
        }
    }
    
    edge_length = max_along_edge - min_along_edge;
}

