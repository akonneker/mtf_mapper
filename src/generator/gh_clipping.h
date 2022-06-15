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
#ifndef GH_CLIPPING_H
#define GH_CLIPPING_H


#include <opencv2/core/core.hpp>
using namespace cv;

#include <vector>
using std::vector;

#include "geom.h"

class Polygon_geom;

//==============================================================================
namespace GH_clipping {

    typedef enum {
        NONE = 0,
        EN   = 1,
        EX   = 2,
        EXEN = 4,
        ENEX = 8
    } traversal_flag;
    
    typedef enum {
        EDGE_OUT = 0,
        EDGE_IN = 1,
        EDGE_ON = 2
    } edge_status;
    
    typedef enum {
        D1 = 1,
        D2,
        D3,
        D4
    } traversal_edge;

    struct gh_vertex {
        double x;
        double y;
        int next;       // we use integers here to avoid multiple dynamic allocations
        int prev;       // and the corresponding pain of freeing up the linked list again
        bool isect;
        traversal_flag flag;
        double alpha;
        int neighbour;
        int next_poly;  // we also use a global vector for storing the polygons ... but there are only two?
        int couple;     // a coupled vertex pointer
        bool cross_change; // special case handled in fig 12 of Kim's paper
    };

                      
    int init_gh_list(vector<gh_vertex>& verts, const vector<cv::Vec2d>& in_verts, int vs, int next_poly, double xoffset=0, double yoffset=0);
    bool gh_phase_one(vector<gh_vertex>& verts, int& vs, int poly0_size, int poly1_size);
    void gh_phase_two(vector<gh_vertex>& verts, const Polygon_geom& b, int first_vert_index, double xoffset=0, double yoffset=0);
    void gh_phase_two_b(vector<gh_vertex>& verts, const Polygon_geom& b, int poly1_start);
    double gh_phase_three(vector<gh_vertex>& verts, int vs, int first_isect_index, vector<Polygon_geom>& polys, bool area_only=false);
    
}

#endif // GH_CLIPPING_H
