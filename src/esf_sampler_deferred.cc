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

#include "include/esf_sampler_deferred.h"

Point2d Esf_sampler_deferred::bracket_minimum(double t0, const Point2d& l, const Point2d& p, const Point2d& pt) {
    double h = 0.01;
    Point2d p_t0 = undistort->transform_point(t0*l + p);
    double d_t0 = norm(p_t0 - pt);
    Point2d p_tph = undistort->transform_point((t0+h)*l + p);
    double d_tph = norm(p_tph - pt);
    
    double a;
    double b;
    double eta;
    
    if (d_t0 > d_tph) {
        // forward step
        a = t0;
        eta = t0 + h;
        while (true) { // TODO: add a way to break out
            h = h*2;
            b = a + h;
            
            Point2d p_b = undistort->transform_point(b*l + p);
            double d_b = norm(p_b - pt);
            Point2d p_eta = undistort->transform_point(eta*l + p);
            double d_eta = norm(p_eta - pt);
            
            if (d_b >= d_eta) {
                return Point2d(a, b);
            }
            a = eta;
            eta = b;
        }
    } else {
        // backward step
        eta = t0;
        b = t0 + h;
        while (true) {
            h = h*2;
            a = b - h;
            
            Point2d p_a = undistort->transform_point(a*l + p);
            double d_a = norm(p_a - pt);
            Point2d p_eta = undistort->transform_point(eta*l + p);
            double d_eta = norm(p_eta - pt);
            
            if (d_a >= d_eta) {
                return Point2d(a, b);
            }
            b = eta;
            eta = a;
        }
    }
}

Point2d Esf_sampler_deferred::derivative(double t0, const Point2d& l, const Point2d& p) {
    const double epsilon = 1e-8;
    Point2d p_t = undistort->transform_point(t0*l + p);
    Point2d p_th = undistort->transform_point((t0+epsilon)*l + p);
    return (1.0/epsilon) * (p_th - p_t);
}

double Esf_sampler_deferred::quadmin(const Point2d& a, const Point2d& b, const Point2d& c) {
    double denom = (b.x - a.x)*(b.y - c.y) - (b.x - c.x)*(b.y - a.y);
    double num = (b.x - a.x)*(b.x - a.x)*(b.y - c.y) - (b.x - c.x)*(b.x - c.x)*(b.y - a.y);
    return b.x - 0.5*num/denom;
}

void Esf_sampler_deferred::sample(Edge_model& edge_model, vector<Ordered_point>& local_ordered, 
    const map<int, scanline>& scanset, double& edge_length,
    const cv::Mat& geom_img, const cv::Mat& sampling_img,
    Bayer::cfa_mask_t cfa_mask) {
    
    cfa_mask = cfa_mask == Bayer::DEFAULT ? default_cfa_mask : cfa_mask;
    
    double max_along_edge = -1e50;
    double min_along_edge = 1e50;
    
    map<int, scanline> m_scanset;
    for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); ++it) {
        int y = it->first;
        for (int x=it->second.start; x <= it->second.end; ++x) {
            cv::Point2i tp = undistort->transform_pixel(x, y);
            auto fit = m_scanset.find(tp.y);
            if (fit == m_scanset.end()) {
                m_scanset.insert(make_pair(tp.y, scanline(tp.x, tp.x)));
            } else {
                fit->second.update(tp.x);
            }
        }
    }
    
    for (map<int, scanline>::const_iterator it=m_scanset.begin(); it != m_scanset.end(); ++it) {
        int y = it->first;
        if (y < border_width || y > sampling_img.rows-1-border_width) continue;
        int rowcode = (y & 1) << 1;
        
        for (int x=it->second.start; x <= it->second.end; ++x) {
            
            if (x < border_width || x > sampling_img.cols-1-border_width) continue;
            
            int code = 1 << ( (rowcode | (x & 1)) ^ 3 );
            if ((code & cfa_mask) == 0) continue;
            
            // transform warped pixel to idealized rectilinear
            cv::Point2d tp = undistort->inverse_transform_point(x, y);
            
            bool inside_image = lrint(tp.x) >= 0 && lrint(tp.x) < geom_img.cols && lrint(tp.y) >= 0 && lrint(tp.y) < geom_img.rows;
            if (!inside_image) continue;
            
            Point2d d = tp - edge_model.get_centroid();
            double perp = d.ddot(edge_model.get_normal()); 
            double par = d.ddot(edge_model.get_direction());
            
            if (!undistort->rectilinear_equivalent()) {
                // 'par' is gamma from the paper
                // apply bracketing
                Point2d pd(x, y);
                Point2d bracketed = bracket_minimum(par, edge_model.get_direction(), edge_model.get_centroid(), pd);
                
                // apply quadratic interpolation
                Point2d p1(bracketed.x, norm(undistort->transform_point(bracketed.x*edge_model.get_direction() + edge_model.get_centroid()) - pd));
                Point2d p2(0.5*(bracketed.x+bracketed.y), norm(undistort->transform_point((0.5*(bracketed.x+bracketed.y))*edge_model.get_direction() + edge_model.get_centroid()) - pd));
                Point2d p3(bracketed.y, norm(undistort->transform_point(bracketed.y*edge_model.get_direction() + edge_model.get_centroid()) - pd));
                double tau_star = quadmin(p1, p2, p3);
                
                // find tangent, then project onto normal
                Point2d tangent = derivative(tau_star, edge_model.get_direction(), edge_model.get_centroid());
                tangent *= 1.0/norm(tangent);
                Point2d lnorm(-tangent.y, tangent.x);
                Point2d ppd = undistort->transform_point(tau_star*edge_model.get_direction() + edge_model.get_centroid());
                
                perp = (pd - ppd).ddot(lnorm);
            }
            if (fabs(perp) < max_dot) {
                local_ordered.push_back(Ordered_point(perp, sampling_img.at<uint16_t>(y,x) ));
                max_along_edge = max(max_along_edge, par);
                min_along_edge = min(min_along_edge, par);
            }
        }
    }
        
    edge_length = max_along_edge - min_along_edge;
}

