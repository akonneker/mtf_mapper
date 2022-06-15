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
#ifndef MULTIPOLYGON_GEOM_H
#define MULTIPOLYGON_GEOM_H

using namespace cv;

#include <vector>
using std::vector;

#include <algorithm>
using std::sort;

#include "geom.h"
#include "polygon_geom.h"

//==============================================================================
class Multipolygon_geom : public Geometry {
  public:

    Multipolygon_geom(double xoff, double yoff, const string& fname, double analogue_scale=1.0) 
    : Geometry(0, 0, 0, Aa_bb(0,0,1,1)), total_vertices(0) {

        FILE* fin = fopen(fname.c_str(), "rt");
        if (!fin) {
            fprintf(stderr, "Error. Could not open polygon geometry file %s. Aborting\n", fname.c_str());
            exit(1);
        }

        // TODO: maybe one day this can read SVG?
        while (!feof(fin)) {
            int nverts = 0;
            int nread = fscanf(fin, "%4d", &nverts);
            if (nread == 1) {
                if (nverts > Polygon_geom::max_verts_per_poly) {
                    printf("Error! Polygon input file [%s] contains a polygon with %d vertices, which\n", fname.c_str(), nverts);
                    printf("\t which exceeds the built-in limit of %d vertices\n", Polygon_geom::max_verts_per_poly);
                }
                vector<cv::Vec2d> verts(nverts);
                int j = 0; // out vertex index
                for (int i=0; i < nverts; i++) {
                    nread = fscanf(fin, "%22lf %22lf", &verts[j][0], &verts[j][1]);
                    verts[j] *= analogue_scale;
                    verts[j][0] += xoff; 
                    verts[j][1] += yoff;
                    // elliminate duplicates
                    if ( (j > 0 && (fabs(verts[j][0] - verts[j-1][0]) < 1e-11 && fabs(verts[j][1] - verts[j-1][1]) < 1e-11)) ||
                         (j > 1 && (fabs(verts[0][0] - verts[j][0]) < 1e-11 && fabs(verts[0][1] - verts[j][1]) < 1e-11)) ) {
                        
                        printf("erasing duplicate vertex (%lf, %lf)\n", verts[j][0], verts[j][1]);
                        verts.erase(--verts.end());
                    } else {
                        j++;
                    }
                }
                total_vertices += j;
                Polygon_geom pol(verts);
                if (!pol.has_ccw_winding()) {
                    pol = Polygon_geom(vector<cv::Vec2d>(pol.bases.rbegin(), pol.bases.rend()));
                }
                parts.push_back(pol);
            }
        }
        
        fclose(fin);

        compute_bounding_box();
        
        own_area = 1;
    }
    
    Multipolygon_geom(void) : total_vertices(0) {
        own_area = 1;
    }
    
    virtual ~Multipolygon_geom(void) {
    }
    
    void compute_bounding_box(void) {

        // compute a bounding box
        bounds.min_y = 1e12;
        bounds.max_x = -1e12;
        bounds.max_y = -1e12;
        bounds.min_x = 1e12;
        for (size_t p=0; p < parts.size(); p++) {
            for (int i=0; i < parts[p].nvertices; i++) {
                bounds.min_y = min(parts[p].bases[i][1], bounds.min_y);
                bounds.max_x = max(parts[p].bases[i][0], bounds.max_x);
                bounds.max_y = max(parts[p].bases[i][1], bounds.max_y);
                bounds.min_x = min(parts[p].bases[i][0], bounds.min_x);
            }
        }
        bounds.area = (bounds.max_x - bounds.min_x) * (bounds.max_y - bounds.min_y);
    }

    bool is_inside(double x, double y) const {
    
        if (!bounds.is_inside(x, y)) {
            return false;
        }
    
        for (size_t p=0; p < parts.size(); p++) {
            if (parts[p].is_inside(x, y)) {
                return true;
            }
        }

        return false;
    }
    
    double intersection_area(const Geometry& b, double xoffset = 0, double yoffset = 0)  const {

        if (!b.bounds.bounds_overlap(bounds, xoffset, yoffset)) {
            return 0;
        }
        
        double total_area = 0;
        for (size_t p=0; p < parts.size(); p++) {
            total_area += parts[p].intersection_area(b, xoffset, yoffset);
        }
        
        return total_area;
    }
    
    void print(void) const {
        for (size_t p=0; p < parts.size(); p++) {
            printf("part %d:\n", int(p));
            parts[p].print();
        }
    }
    

    vector<Polygon_geom> parts;
    int total_vertices;
};

#endif // MULTIPOLYGON_GEOM_H
