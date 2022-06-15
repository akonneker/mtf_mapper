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
#ifndef EDGE_RECORD_H
#define EDGE_RECORD_H

#include "include/logger.h"
#include "common_types.h"


class Edge_record {
  public:
    Edge_record(void) : pooled(false) {
    }

    static void pool_edges(Edge_record& a, Edge_record& b) {
        Edge_record m;
        m.points.resize(a.points.size() + b.points.size());
        m.weights.resize(a.weights.size() + b.weights.size());
        
        for (size_t i=0; i < a.points.size(); i++) {
            m.points[i].first = a.points[i].first - a.centroid.x;
            m.points[i].second = a.points[i].second - a.centroid.y;
            m.weights[i] = a.weights[i];
        }
        size_t as = a.points.size();
        for (size_t i=0; i < b.points.size(); i++) {
            m.points[i+as].first = b.points[i].first - b.centroid.x;
            m.points[i+as].second = b.points[i].second - b.centroid.y;
            m.weights[i+as] = b.weights[i];
        }
        
        m.compute_eigenvector_angle();
        a.angle = m.angle;
        b.angle = m.angle;
        a.pooled = true;
        b.pooled = true;
    }

    bool compatible(const Edge_record& b) {
        double Z;
        
        if (fabs(slope) > 2) {
            Z = (1.0/slope - 1.0/b.slope)/sqrt(sB*sB + b.sB*b.sB);
        } else {
            Z = (slope - b.slope)/sqrt(sB*sB + b.sB*b.sB);
        }

        return fabs(Z) < 1.66; // ~90% confidence interval on t-distribution with ~80 dof
    }
    
    inline double relative_orientation(const Edge_record& b) {
        Point2d d1(cos(angle), sin(angle));
        Point2d d2(cos(b.angle), sin(b.angle));
        
        return fabs(d1.x*d2.x + d1.y*d2.y);
    }
    
    void add_point(double x, double y, double gx, double gy) {
        points.push_back(make_pair(x, y));
        double mag = (gx*gx + gy*gy);
        weights.push_back(mag);
    }
    
    pair<double, double> compute_eigenvector_angle(void) {
        double covxx = 0;
        double covxy = 0;
        double covyy = 0;
        
        centroid.x = 0;
        centroid.y = 0;
        wsum = 0;
        for (size_t i=0; i < points.size(); i++) {
            double w = weights[i];
            
            if (w > 0) {
                double temp = w + wsum;
                double delta_x = points[i].first - centroid.x;
                double delta_y = points[i].second - centroid.y;
                double rx = delta_x * w / temp;
                double ry = delta_y * w / temp;
                centroid.x += rx;
                centroid.y += ry;
                
                covxx += wsum * delta_x * rx;
                covyy += wsum * delta_y * ry;
                covxy += wsum * delta_x * ry;
                
                wsum = temp;
            }
        }
        
        covxx /= wsum;
        covxy /= wsum;
        covyy /= wsum;
        
        
        // char. poly: l^2 - (a+d)l + (ad-bc) = 0
        // thus l^2 -(tr)l + det = 0
        double tr = covxx + covyy;
        double det = covxx*covyy - covxy*covxy;
        
        double pa=1.0;
        double pb=-tr;
        double pc=det;
        
        double q = -0.5 * (pb + sgn(pb)*sqrt(pb*pb - 4*pa*pc) );
        double l1 = q/pa;
        double l2 = pc / q;
        
        double l = max(l1,l2);
        assert(l >= 0);
        
        double ev[2];
        if (fabs(covxy) > 1e-10) {
            ev[0] = l - covyy;
            ev[1] = covxy;
            slope = ev[0] / ev[1]; // TODO: check this?
        } else {
            logger.info("%s\n", "Warning: edge is perfectly horizontal / vertical");
            if (covxx > covyy) {
                ev[0] = 1;
                ev[1] = 0;
            } else {
                ev[0] = 0;
                ev[1] = 1;
            }
            slope = 0;
        }
        
        angle = atan2(-ev[0], ev[1]);
        
        return make_pair(
            sqrt(max(fabs(l1), fabs(l2))),
            sqrt(min(fabs(l1), fabs(l2)))
        );
    }

    bool reduce(void) { // compute orientation, and remove weak points
        if (weights.size() < 20) { // not enough points to really extract an edge here
            angle = 0;
            rsq = 1.0;
            return false;
        }
        
        // TODO: could look for peaks etc. here?
        
        renormalize_weights();
        vector<double> inweights(weights);
        pair<double,double> dims = compute_eigenvector_angle();
        Point2d dir(cos(angle), sin(angle));
        
        vector<double> histo(2*16*8, 0);
        double total_weight = 0;
        
        for (size_t i=0; i < points.size(); i++) {
            double dx = points[i].first - centroid.x;
            double dy = points[i].second - centroid.y;
            
            double dot = dx * dir.x + dy * dir.y;
            
            // dot is "across edge" direction
            // histogram of dot vs weight?
            double idot = lrint(dot * 8 + 16*8);
            if (idot >= 3 && idot < (int)histo.size() - 3) {
                histo[idot] += weights[i];
                histo[idot-1] += weights[i];
                histo[idot+1] += weights[i];
                histo[idot-2] += weights[i];
                histo[idot+2] += weights[i];
                histo[idot-3] += weights[i];
                histo[idot+3] += weights[i];
                
                total_weight += weights[i]*7;
            }
        }
        // find peak value in central 5 pixel band
        size_t central_idx = 16*8;
        for (size_t j=(-5*8+16*8); j <= (5*8+16*8); j++) {
            if (histo[j] > histo[central_idx]) {
                central_idx = j;
            }
        }
        // now find 5% cut-off on both sides of peak
        size_t lower5p = central_idx-1;
        while (lower5p > 1 && histo[lower5p] > 0.05*histo[central_idx]) lower5p--;
        size_t upper5p = central_idx+1;
        while (upper5p < histo.size()-1 && histo[upper5p] > 0.05*histo[central_idx]) upper5p++;
        
        // move outwards a little, then look for another value greater than the 5% threshold
        size_t rise_lower5p = (size_t)max(int(1), int(lower5p) - 8);
        while (rise_lower5p > 1 && histo[rise_lower5p] <= histo[lower5p]) rise_lower5p--;
        size_t rise_upper5p = min(histo.size()-1, upper5p+8);
        while (rise_upper5p < histo.size()-1 && histo[rise_upper5p] <= histo[upper5p]) rise_upper5p++;
        
        double csum = 0;
        for (size_t i = 0; i < histo.size(); i++) {
            csum += histo[i];
            histo[i] = csum;
        }
        total_weight = csum;
        int p10idx = 0;
        for (size_t i = 1; i < histo.size(); i++) {
            if (fabs(histo[i] - 0.1*csum) < fabs(histo[p10idx] - 0.1*csum)) {
                p10idx = i;
            }
        }
        int p90idx = histo.size()-1;
        for (size_t i = histo.size() - 2; i > 0; i--) {
            if (fabs(histo[i] - 0.9*csum) < fabs(histo[p90idx] - 0.9*csum)) {
                p90idx = i;
            }
        }
        
        
        
        double lower = p10idx / 8.0 - 16.0;
        double upper = p90idx / 8.0 - 16.0;
        double span = upper - lower;
        
        lower -= span * 0.7;
        upper += span * 0.7;
        
        lower = max(double(rise_lower5p + lower5p)/16.0 - 16.0, lower);
        upper = min(double(rise_upper5p + upper5p)/16.0 - 16.0, upper);
        
        // trim weights to 10%-90% region?
        for (size_t i=0; i < points.size(); i++) {
            double dx = points[i].first - centroid.x;
            double dy = points[i].second - centroid.y;
            
            double dot = dx * dir.x + dy * dir.y;
            weights[i] = 0;
            if (dot >= lower && dot <= upper) {
                weights[i] =  SQR(SQR(inweights[i])) * (1.0 / (10.0 + fabs(dot)));
            }
        }
        dims = compute_eigenvector_angle();
        
        radii = dims;
        
        sB = radii.second / (radii.first*sqrt(wsum));
        
        rsq = 0.0;  // TODO: what is the appropriate measure of uncertainty in the angle estimate?
        
        return true;
    }

    inline bool is_pooled(void) {
        return pooled;
    }
    
    void clear(void) {
        weights.clear();
        points.clear();
    }

    double slope;
    double angle;
    double rsq;
    Point2d  centroid;
    pair<double, double>  radii;

  private:
  
    void renormalize_weights(void) {
        double maxw = 0;
        for (size_t i=0; i < weights.size(); i++) {
            maxw = max(weights[i], maxw);
        }
        if (maxw > 0) {
            for (size_t i=0; i < weights.size(); i++) {
                weights[i] /= maxw;
            }
        }
    }

    vector< pair<double, double> > points;
    vector< double > weights;
    
    double wsum;

    double sB; // standard error in slope estimate

    bool pooled;
};


#endif 


