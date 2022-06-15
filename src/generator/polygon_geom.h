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
#ifndef POLYGON_GEOM_H
#define POLYGON_GEOM_H

#ifndef M_PI
    #define M_PI 3.14159265358979
#endif

using namespace cv;

#include <vector>
using std::vector;

#include "geom.h"
#include "gh_clipping.h"

class Multipolygon_geom;

typedef enum {
    INSIDE = 1,
    OUTSIDE = 0,
    ON = 2
} Vertex_type;

//==============================================================================
class Polygon_geom : public Geometry {
  public:
    static const int max_verts_per_poly = 400; // max number of point in intersected polygon

    friend class Multipolygon_geom;
    
    
  
    Polygon_geom(double cx=0, double cy=0, double width=1, double height=1, double angle=0, int nverts=4) 
    : Geometry(cx, cy, 0), nvertices(nverts), convex(false) {

        construct_regular_polygon(width, height, angle);
        precompute_point_test();
        own_area = compute_area();
        //printf("built polygon with %d vertices, area=%lf\n\n", nverts, own_area);
    }

    Polygon_geom(const vector<cv::Vec2d>& verts) 
    : convex(false) {
    
        // although it is possible for an arbitrary polygon specified through
        // "verts" to be convex, we choose to not take chances here
        // so the slower (but concave-capable) point-in-poly test will be used

        nvertices = verts.size();

        bases = verts;
        normals = vector<cv::Vec2d>(nvertices);

        // compute normals
        int prev = nvertices - 1;
        for (int i=0; i < nvertices; i++) {
            cv::Vec2d d = bases[i] - bases[prev];
            normals[i][0] = -d[1] / norm(d);
            normals[i][1] = d[0] / norm(d);
            prev = (prev + 1) % nvertices;
        }
        
        convex = is_convex();
        if (!convex) {
            precompute_point_test();
        }

        compute_bounding_box();
        
        own_area = compute_area();
    }
    
    Polygon_geom(const Aa_bb& b) // construct polygon geometry from bounding box
    : Geometry(b), nvertices(4), convex(true) { 
    
        bases = vector<cv::Vec2d>(4);
    
        bases[0] = cv::Vec2d(bounds.max_x, bounds.max_y);
        bases[1] = cv::Vec2d(bounds.max_x, bounds.min_y);
        bases[2] = cv::Vec2d(bounds.min_x, bounds.min_y);
        bases[3] = cv::Vec2d(bounds.min_x, bounds.max_y);
        
        normals = vector<cv::Vec2d>(4);
        
        // compute normals
        int prev = 3;
        for (int i=0; i < 4; i++) {
            cv::Vec2d d = bases[i] - bases[prev];
            normals[i][0] = -d[1] / norm(d);
            normals[i][1] = d[0] / norm(d);
            prev = (prev + 1) % 4;
        }
        
    }

    virtual ~Polygon_geom(void) {
    }

    void rebuild(void) {
	    // compute normals
        int prev = nvertices - 1;
        for (int i=0; i < nvertices; i++) {
            cv::Vec2d d = bases[i] - bases[prev];
            normals[i][0] = -d[1] / norm(d);
            normals[i][1] = d[0] / norm(d);
            prev = (prev + 1) % nvertices;
        }
        
        precompute_point_test();

        compute_bounding_box();
        
        own_area = compute_area();
    }

    void construct_regular_polygon(double width, double height, double angle) {
        convex = true;
        
        bases   = vector<cv::Vec2d>(nvertices);
        normals = vector<cv::Vec2d>(nvertices);

        assert(nvertices >= 3);
        for (int i=0; i < nvertices; i++) {
            double phi = i*M_PI*2/double(nvertices);
            bases[i][0] = width/2*cos(angle+phi - M_PI/4.0) + cx;
            bases[i][1] = height/2*sin(angle+phi - M_PI/4.0) + cy;
        }
        // yes, I know I should be able to compute the normals directly, but this still works, ok :)
        int prev = nvertices - 1;
        for (int i=0; i < nvertices; i++) {
            cv::Vec2d d = bases[i] - bases[prev];
            normals[i][0] = -d[1] / norm(d);
            normals[i][1] = d[0] / norm(d);
            prev = i;
        }

        compute_bounding_box();
        
    }

    void compute_bounding_box(void) {
        // compute a bounding box
        bounds.min_y = 1e12;
        bounds.max_x = -1e12;
        bounds.max_y = -1e12;
        bounds.min_x = 1e12;
        for (size_t i=0; i < bases.size(); i++) {
            bounds.min_y = min(bases[i][1], bounds.min_y);
            bounds.max_x = max(bases[i][0], bounds.max_x);
            bounds.max_y = max(bases[i][1], bounds.max_y);
            bounds.min_x = min(bases[i][0], bounds.min_x);
        }
        bounds.area = (bounds.max_x - bounds.min_x) * (bounds.max_y - bounds.min_y);
    }
    
    void precompute_point_test(void) {
        constant = vector<double>(nvertices);
        multiple = vector<double>(nvertices);
        int j = (int)bases.size() - 1;

        for(int i=0; i < (int)bases.size(); i++) {
            if(bases[j][1] == bases[i][1]) {
                constant[i] = bases[i][0];
                multiple[i] = 0; 
            } else {
                constant[i] = bases[i][0] - (bases[i][1]*bases[j][0]) / (bases[j][1] - bases[i][1]) + 
                        (bases[i][1]*bases[i][0]) / (bases[j][1] - bases[i][1]);

                multiple[i] = (bases[j][0] - bases[i][0]) / (bases[j][1] - bases[i][1]); 
            }
            j = i; 
        }
    }
    
    // NOTE: indeterminate result if point falls exactly on boundary
    // the "classify" test can correctly detect boundary cases
    // but it is too slow for general point-in-poly tests
    bool is_inside(double x, double y) const {
        if (!bounds.is_inside(x, y)) {
            return false;
        }
    
        if (convex) {
            return convex_is_inside(x, y);
        }

        int j = bases.size() - 1;
        bool  oddNodes = false;

        for (int i=0; i < (int)bases.size(); i++) {
            if ( (bases[i][1] < y && bases[j][1] >= y) || 
                 (bases[j][1] < y && bases[i][1] >= y) ) {

                oddNodes ^= (y*multiple[i] + constant[i]) < x; 
            }
            j = i; 
        }

        return oddNodes; 
    }

    
    inline bool convex_is_inside(double x, double y) const {
        for (int i=0; i < nvertices; i++) {
        
            if ( (normals[i][0]*(x - bases[i][0]) + 
                 normals[i][1]*(y - bases[i][1])) < 0) {
                 
                return false;
            }
            
        }
        return true; // failed to elliminate point using any of the half-spaces
    }
    
    inline double det(int i, double x, double y, bool& edge) const {
        int next = (i+1)%bases.size();
        double d = (bases[i][0] - x) * (bases[next][1] - y) -
                   (bases[next][0] - x) * (bases[i][1] - y);
                   
        if (d == 0) {
            edge = true;
        } else {
            edge = false;
        }
        return d;
    }
    
    #define crossing ( (bases[i][1] < y) != (bases[next][1] < y) )
    #define right_crossing  (det(i, x, y, edge) > 0) == (bases[next][1] > bases[i][1])  
    #define modify_w ( w += ((bases[next][1] > bases[i][1]) ? 2 : 0) - 1  )
   
    inline Vertex_type classify(double x, double y) const {
        bool edge = false;
        if (x == bases[0][0] && y == bases[0][1]) {
            return ON; // on first vertex
        }
        int w = 0;
        for (int i=0; i < (int)bases.size(); i++) {
            int next = (i+1)%bases.size();
            if (bases[next][1] == y) {
                if (bases[next][0] == x) {
                    return ON; // on a vertex
                } else {
                    if ( bases[i][1] == y && ((bases[next][0] > x) == (bases[i][0] < x)) ) {
                        return ON; // on an edge
                    }
                }
            }
            if (crossing) {
                if (bases[i][0] >= x) {
                    if (bases[next][0] > x) {
                        modify_w;
                    } else {
                        if (right_crossing) {
                            modify_w;
                        }
                        if (edge) return ON; // non-horizontal edge?
                    }
                } else {
                    if (bases[next][0] > x) {
                        if (right_crossing) {
                            modify_w;
                        }
                        if (edge) return ON; // non-horizontal edge
                    }
                }
            }
        }
        return w == 0 ? OUTSIDE : INSIDE; 
    }
    
    bool is_convex(void) const {
        int flag = 0;
        int n = bases.size();

        if (n < 3) {
            return true; // define degenerates as convex
        }

        for (int i=0; i < n;i++) {
            int j = (i + 1) % n;
            int k = (i + 2) % n;
            
            double z  = (bases[j][0] - bases[i][0]) * (bases[k][1] - bases[j][1]);
            z -= (bases[j][1] - bases[i][1]) * (bases[k][0] - bases[j][0]);
            
            if (z < 0) {
                flag |= 1;
            } else {
                if (z > 0) {
                    flag |= 2;
                }
            }
            if (flag == 3) {
                return false;
            }
        }
        if (flag != 0) {
            return true;
        }
        return false;
    }
    
    double intersection_area(const Geometry& ib, double xoffset = 0, double yoffset = 0)  const {
    
        double points_x[max_verts_per_poly];
        double points_y[max_verts_per_poly];
        int points_len = 0;
        
        // TODO: this will probably throw an exception if you try to
        // pass a multipolygon as the photosite geometry. Don't do it!
        const Polygon_geom& b = dynamic_cast<const Polygon_geom&>(ib);

        // first check for overlap between the bounding boxes
        if (!b.bounds.bounds_overlap(bounds, xoffset, yoffset)) {
            return 0;
        }
        
        // if poly b is moderately complex, first check if the bounding box of
        // b intersects the real geometry (since we know the the bounds overlap),
        // then check if this is total overlap
        if (b.nvertices >= 6) { 
            b.intersect_bounds(points_x, points_y, points_len, *this, xoffset, yoffset);
            double i_bb_area = compute_area(points_x, points_y, points_len);
            
            if (fabs(i_bb_area) < 1e-11) { // bounds overlapped, but actual overlap is null
                return 0; 
            } else {
                if (fabs(i_bb_area - b.bounds.area) < 1e-11) { // full overlap, no need to check further
                    return b.own_area; // full photosite area
                } 
            }

            // partial intersection, so reset and perform 
            // intersection again using full geometry
            points_len = 0;
        }
        
        
        if (b.convex) {
            // if we know the clipper (photosite) is convex, use sutherland-hodgeman
            b.intersect(points_x, points_y, points_len, *this, xoffset, yoffset);  
            return compute_area(points_x, points_y, points_len);
        } else {
            // otherwise use Greiner-Horman, which is about 6x slower, but works for concave photosites
            return b.intersect_greiner_horman_area(*this, xoffset, yoffset);
        }
    }
    
    inline bool t_intersect(double& pix, double& piy, 
                     const double& v1x, const double& v1y,
                     const double& d1x, const double& d1y,
                     const double& v2x, const double& v2y,
                     const double& d2x, const double& d2y) const {
                   
        double denom = (d2y*d1x - d2x*d1y);
        
        if (fabs(denom) < 1e-11) {
           // this happens when the lines are parallel
           // the caller handles this correctly by not
           // adding an additional intersection point
           return false;
        }
        
        double u = (d2x*(v1y - v2y) - d2y*(v1x - v2x)) / denom;
        
        pix = v1x + u*d1x;
        piy = v1y + u*d1y;
                   
        return true;               
    }

    void print(void) const {
        for (int i=0; i < nvertices; i++) {
            printf("\t%lf %lf\n", bases[i][0], bases[i][1]);
        }
    }
    
    double compute_area(double* points_x, double* points_y, int points_len) const {
        double A = 0;
        for (int i=0; i < points_len; i++) {
            int ni = (i+1) % points_len;
            A += points_x[i]*points_y[ni] - points_x[ni]*points_y[i];
        }
        return 0.5 * fabs(A);
    }
    
    double compute_area(void) const {
        double A = 0;
        for (int i=0; i < nvertices; i++) {
            int ni = (i+1) % nvertices;
            A += bases[i][0]*bases[ni][1] - bases[ni][0]*bases[i][1];
        }
        return 0.5 * fabs(A);
    }
    
    bool has_ccw_winding(void) const {
        double A = 0;
        for (int i=0; i < nvertices; i++) {
            int ni = (i+1) % nvertices;
            A += bases[i][0]*bases[ni][1] - bases[ni][0]*bases[i][1];
        }
        return A > 0;
    }

    void intersect(double* points_x, double* points_y, int& points_len, 
        const Polygon_geom& b, double xoffset = 0, double yoffset = 0) const {

        for (int i=0; i < b.nvertices; i++) {
            points_x[i] = b.bases[i][0];
            points_y[i] = b.bases[i][1];
            points_len++;
        }
        
        for (int e=0; e < nvertices; e++) {
        
            int ne = (e + 1) % nvertices;
            double Px = bases[e][0] + xoffset;
            double Py = bases[e][1] + yoffset;
            double Dx = bases[ne][0] - bases[e][0];
            double Dy = bases[ne][1] - bases[e][1];
        
            intersect_core(points_x, points_y, points_len, Px, Py, Dx, Dy);
        }
        
    }
    
    
    void intersect_bounds(double* points_x, double* points_y, int& points_len, 
        const Polygon_geom& b, double xoffset = 0, double yoffset = 0) const { // special case for bounding boxes
        
        for (int i=0; i < b.nvertices; i++) { 
            points_x[i] = b.bases[i][0];
            points_y[i] = b.bases[i][1];
            points_len++;
        }
        
        intersect_core(points_x, points_y, points_len, 
            bounds.max_x + xoffset, bounds.max_y + yoffset, 
            -1, 0
        );
        
        intersect_core(points_x, points_y, points_len, 
            bounds.max_x + xoffset, bounds.min_y + yoffset, 
            0, 1
        );
        
        intersect_core(points_x, points_y, points_len, 
            bounds.min_x + xoffset, bounds.min_y + yoffset, 
            1, 0
        );
        
        intersect_core(points_x, points_y, points_len, 
            bounds.min_x + xoffset, bounds.max_y + yoffset, 
            0, -1
        );
    }
    

    // this looks like the Sutherland-Hodgman algorithm
    // the clipping polygon must be convex (req. by SH algo)
    // will produce overlapping edges if a concave point exists
    // outside of the clipping polygon. these overlapping
    // edges cause no harm, because they have zero area (which
    // seems to work out fine with compute_area())
    // see the Greiner-Hormann algorithm below if you
    // are interested in the correct geometry for
    // concave/concave or convex/concave cases
    inline void intersect_core(double* inpoints_x, double* inpoints_y, int& in_len, 
        double Px, double Py, double Dx, double Dy) const {
                        
        double Nx = -Dy;
        double Ny = Dx;

        double outpoints_x[max_verts_per_poly];
        double outpoints_y[max_verts_per_poly];
        int out_idx = 0;
        
        double Sx = inpoints_x[in_len - 1];
        double Sy = inpoints_y[in_len - 1];
        
        double dSx = Sx - Px;
        double dSy = Sy - Py;
        
        for (size_t i=0; i < size_t(in_len); i++) {
            const double& Ex = inpoints_x[i];
            const double& Ey = inpoints_y[i];
            
            double dEx = Ex - Px;
            double dEy = Ey - Py;
            
            if ( (Nx*dEx + Ny*dEy) >= 0) {
                if ( (Nx*dSx + Ny*dSy) < 0) {
                    t_intersect(
                        outpoints_x[out_idx], outpoints_y[out_idx],
                        Sx, Sy,
                        Ex - Sx, Ey - Sy,
                        Px, Py,
                        Dx, Dy
                    );
                    out_idx++;
                }
                outpoints_x[out_idx] = Ex;
                outpoints_y[out_idx] = Ey;
                out_idx++;
            } else {  
                if ( (Nx*dSx + Ny*dSy) >= 0) {
                    t_intersect(
                        outpoints_x[out_idx], outpoints_y[out_idx],
                        Sx, Sy,
                        Ex - Sx, Ey - Sy,
                        Px, Py,
                        Dx, Dy
                    );
                    out_idx++;
                }
            }
            Sx = Ex;
            Sy = Ey;
            dSx = dEx;
            dSy = dEy;
        }
    
        in_len = out_idx;
        memcpy(inpoints_x, outpoints_x, sizeof(double)*out_idx);
        memcpy(inpoints_y, outpoints_y, sizeof(double)*out_idx);
    }
    
    
    // slower polygon clipping algorithm, but this one should handle concave-concave
    // clipping, and it should also avoid creating degenerate parts
    vector<Polygon_geom> intersect_greiner_horman(const Polygon_geom& b) const {
        vector<Polygon_geom> polys;

        assert(has_ccw_winding() == b.has_ccw_winding());
        
        vector<GH_clipping::gh_vertex> verts(nvertices * b.nvertices*2 + nvertices + b.nvertices);
        
        // populate the verts vector with the two polys
        int poly1_start = GH_clipping::init_gh_list(verts, bases, 0, 1);
        int vs = GH_clipping::init_gh_list(verts, b.bases, poly1_start, -1);
        
        int vs_before_intersections = vs;
        bool found_intersections = GH_clipping::gh_phase_one(verts, vs, bases.size(), b.bases.size());
        
        if (vs == vs_before_intersections && !found_intersections) {
        
            bool all_on = true;
            for (size_t p=0; p < bases.size(); p++) {
                int cl = b.classify(bases[p][0], bases[p][1]);
                if (cl == OUTSIDE) {
                    all_on = false;
                }
            }
            
            if (all_on) {
                // *this must be entirely within b, so return *this
                polys.push_back(*this);
                return polys;
            } else {
                // maybe b is entirely inside *this?
                bool all_in = true;
                for (size_t p=0; p < b.bases.size(); p++) {
                    int cl = classify(b.bases[p][0], b.bases[p][1]);
                    if (cl == OUTSIDE) {
                        all_in = false;
                    }
                }
                
                if (all_in) {
                    polys.push_back(b);
                    return polys;
                }
            
                // *this is entirely outside b, so return empty list
                return polys;
            }
        }
        
        // first process C
        GH_clipping::gh_phase_two(verts, *this, poly1_start);
        // then process S with alternate version 
        GH_clipping::gh_phase_two_b(verts, b, poly1_start);
        
        GH_clipping::gh_phase_three(verts, vs, vs_before_intersections, polys);
        
        for (size_t k=0; k < polys.size(); k++){
            if (!polys[k].has_ccw_winding()) {
                polys[k] = Polygon_geom(vector<cv::Vec2d>(polys[k].bases.rbegin(), polys[k].bases.rend()));
            }
        }
        
        return polys;
    }
    
    double intersect_greiner_horman_area(const Polygon_geom& b, double xoffset=0, double yoffset=0) const {

        assert(has_ccw_winding() == b.has_ccw_winding());
        
        vector<GH_clipping::gh_vertex> verts(nvertices * b.nvertices*2 + nvertices + b.nvertices); // this malloc will slow things down a lot
        
        // populate the verts vector with the two polys
        // xoffset and yoffset are typically used on a photosite aperture polygon
        // so *this is by default considered to be 
        int poly1_start = GH_clipping::init_gh_list(verts, bases, 0, 1, xoffset, yoffset);
        int vs = GH_clipping::init_gh_list(verts, b.bases, poly1_start, -1);
        
        int vs_before_intersections = vs;
        GH_clipping::gh_phase_one(verts, vs, bases.size(), b.bases.size());
        
        if (vs == vs_before_intersections) {
        
            bool all_on = true;
            for (size_t p=0; p < bases.size(); p++) {
                int cl = b.classify(bases[p][0] + xoffset, bases[p][1] + yoffset);
                if (cl == OUTSIDE) {
                    all_on = false;
                }
            }
            
            if (all_on) {
                // *this must be entirely within b, so return *this
                return own_area;
            } else {
                // maybe b is entirely inside *this?
                bool all_in = true;
                for (size_t p=0; p < b.bases.size(); p++) {
                    int cl = classify(b.bases[p][0] - xoffset, b.bases[p][1] - yoffset);
                    if (cl == OUTSIDE) {
                        all_in = false;
                    }
                }
                
                if (all_in) {
                    return b.own_area;
                }
            
                // *this is entirely outside b, so return empty list
                return 0; 
            }
        }
        
        
        // first process C
        GH_clipping::gh_phase_two(verts, *this, poly1_start, xoffset, yoffset); // set_traversal modified
        // then process S with alternate version 
        GH_clipping::gh_phase_two_b(verts, b, poly1_start); // this way round is safe, because b has not been modified
        
        vector<Polygon_geom> polys; // is this slow?
        
        return GH_clipping::gh_phase_three(verts, vs, vs_before_intersections, polys, true);
    }

    vector<cv::Vec2d> normals;
    vector<cv::Vec2d> bases;

    vector<double> constant;
    vector<double> multiple;
    
    int nvertices;

    bool convex;
};

#endif // RENDER_H
