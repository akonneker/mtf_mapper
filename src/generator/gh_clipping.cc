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

#include "gh_clipping.h"
#include "polygon_geom.h"

namespace GH_clipping {

inline static double cross(const cv::Vec2d& v0, const cv::Vec2d& v1, const cv::Vec2d& v2) {
    double ds_x = v1[0] - v0[0];
    double ds_y = v1[1] - v0[1];
    double dc_x = v2[0] - v0[0];
    double dc_y = v2[1] - v0[1];
                                      
    return (dc_y*ds_x - dc_x*ds_y);
}

static inline bool t_intersect(const gh_vertex& s0, const gh_vertex& s1,
                        const gh_vertex& c0, const gh_vertex& c1,
                        gh_vertex& is, gh_vertex& ic) {
                  
    double ds_x = s1.x - s0.x;
    double ds_y = s1.y - s0.y;
    double dc_x = c1.x - c0.x;
    double dc_y = c1.y - c0.y;
                                      
    double denom = (dc_y*ds_x - dc_x*ds_y);
                                                      
    if (fabs(denom) < 1e-11) {
        return false;
    }
    
    denom = 1.0/denom;
                                                                                                                          
    is.alpha = (dc_x*(s0.y - c0.y) - dc_y*(s0.x - c0.x)) * denom;
    ic.alpha = -(ds_x*(c0.y - s0.y) - ds_y*(c0.x - s0.x)) * denom;
    
    
    if (is.alpha < -1e-12 || is.alpha > (1+1e-12) ||  
        ic.alpha < -1e-12 || ic.alpha > (1+1e-12)) {
    
        return false;
    }
    
    is.x = s0.x + is.alpha*ds_x;
    is.y = s0.y + is.alpha*ds_y;
    is.isect = true;
    
    ic.x = c0.x + ic.alpha*dc_x;
    ic.y = c0.y + ic.alpha*dc_y;
    ic.isect = true;
    
    
    return true;        
}


static void insert_vertex(vector<gh_vertex>& verts, int v1, int next, int vs) {
    int current = v1;
        
    if (verts[current].next != next) {
        do { 
            if (verts[vs].alpha > verts[verts[current].next].alpha) {
                current = verts[current].next;
            } else {
                break;
            }
        } while (verts[vs].alpha > verts[verts[current].next].alpha);
    }
    
    verts[vs].next = verts[current].next;
    verts[vs].prev = current;
    verts[verts[current].next].prev = vs;
    verts[current].next = vs;
    verts[vs].flag = NONE;
    verts[vs].couple = -1;
    verts[vs].cross_change = false;
}

static inline bool same(const gh_vertex& a, const gh_vertex& b) {
    return fabs(a.x - b.x) < 1e-11 && fabs(a.y - b.y) < 1e-11;
}

int check_degeneracy(vector<gh_vertex>& verts, int vs, int si, int ci, int nsi, int nci) {
    int& si_prev = verts[si].prev;
    int& si_next = nsi; 
    
    int& ci_prev = verts[ci].prev;
    int& ci_next = nci;
    
    if (same(verts[vs], verts[si_prev])) {
        return si_prev;
    }
    
    if (same(verts[vs], verts[ci_prev])) {
        return ci_prev;
    }
    
    if (same(verts[vs], verts[si])) {
        return si;
    }
    
    if (same(verts[vs], verts[ci])) {
        return ci;
    }
    
    if (same(verts[vs], verts[si_next])) {
        return si_next;
    }
    
    if (same(verts[vs], verts[ci_next])) {
        return ci_next;
    }
    
    return -1;
}


int init_gh_list(vector<gh_vertex>& verts, const vector<cv::Vec2d>& in_verts, int vs, int next_poly,
    double xoffset, double yoffset) {
    for (int v=0; v < (int)in_verts.size(); v++) {
        verts[v+vs].x = in_verts[v][0] + xoffset;
        verts[v+vs].y = in_verts[v][1] + yoffset;
        verts[v+vs].next = (v+1) % in_verts.size() + vs;
        verts[v+vs].prev = (v + in_verts.size() - 1) % in_verts.size() + vs;
        verts[v+vs].isect = false;   // this it not an intersection, so it is a vertex
        verts[v+vs].flag = NONE;     // thus, en is the inside/outside flag???
        verts[v+vs].neighbour = -1;  // no matching vertex in neighbour yet
        verts[v+vs].next_poly = next_poly; // no next poly yet
        verts[v+vs].alpha = 10;      // any value > 1 indicates no intersection, which is the default
        verts[v+vs].couple = -1;     // not coupled at first
        verts[v+vs].cross_change = false;  // and not a cross_change case either
        
    }
    return vs + in_verts.size(); // index of next free vertex entry
}



bool gh_phase_one(vector<gh_vertex>& verts, int& vs, int poly0_size, int poly1_size) {
    bool some_intersections = false;
    // phase one
    for (int si=0; si < poly0_size; si++) {
        int nsi = (si + 1) % poly0_size;
        for (int ci=0; ci < poly1_size; ci++) {
            int nci = (ci + 1) % poly1_size;
            
            bool isect = t_intersect(
                verts[si], verts[nsi], 
                verts[ci+poly0_size], verts[nci+poly0_size], 
                verts[vs], verts[vs+1]
            );
            
            // if this is identical to the previous intersection, just skip over it
            if (isect &&
                vs - 2 > (poly1_size + poly0_size) &&
                same(verts[vs], verts[vs-2])) {
                isect = false;
            }
            
            some_intersections |= isect;
            
            if (isect) { 
                verts[vs].neighbour = vs+1;
                verts[vs].next_poly = 1;
                verts[vs+1].neighbour = vs;
                verts[vs+1].next_poly = -7;
                
                int deg = check_degeneracy(verts, vs, si, ci+poly0_size, nsi, nci+poly0_size);
                
                if (deg >= 0) { // intersection already exists
                    if (same(verts[si], verts[ci+poly0_size])) {
                        verts[si].neighbour = ci+poly0_size;
                        verts[ci+poly0_size].neighbour = si;
                        verts[si].isect = 1;
                        verts[ci+poly0_size].isect = 1;
                        continue;
                    }
                    if (verts[deg].neighbour == -1) { // if vertex already has a neighbour, just skip
                        if (deg < poly0_size) { // match comes from S, so add intersection to C
                        
                            bool found = false;
                            
                            int found_idx = poly0_size;
                            for (; found_idx < (poly0_size+poly1_size); found_idx++) {
                                if (same(verts[found_idx], verts[deg])) {
                                    found = true;
                                    break;
                                }
                            }
                            
                            int dest = vs+1;
                            
                            if (found) {
                                dest = found_idx;
                            } else {
                                insert_vertex(verts, ci+poly0_size, nci+poly0_size, vs+1);
                            }
                            
                            verts[dest].neighbour = deg;
                            verts[dest].isect = 1;
                            
                            verts[deg].neighbour = dest;
                            verts[deg].isect = 1;
                            
                            verts[vs].isect = 0; // disable redundant vertex
                            verts[vs].next_poly = -20;
                            
                            if (!found) vs += 2;
                        } else { // match comes from C, so add intersection to S
                            
                            bool found = false;
                            int found_idx = 0;
                            for (; found_idx < poly0_size; found_idx++) {
                                if (same(verts[found_idx], verts[deg])) {
                                    found = true;
                                    break;
                                }
                            }
                            
                            int dest = vs;
                            
                            if (found) {
                                dest = found_idx;
                            } else {
                                insert_vertex(verts, si, nsi, vs);
                            }
                            
                            verts[deg].neighbour = dest;
                            verts[deg].isect = 1;
                            
                            verts[dest].neighbour = deg;
                            verts[dest].isect = 1;
                              
                            verts[vs+1].isect = 0; // disable redundant vertex
                            verts[vs+1].next_poly = -20;
                            
                            if (!found) vs += 2;
                            
                        }
                    } 
                } else {
                    insert_vertex(verts, si, nsi, vs);
                    insert_vertex(verts, ci+poly0_size, nci+poly0_size, vs+1);
                    vs += 2;
                }
                
            }
        }
    }
    return some_intersections;
}

inline static double triangle_cross(const gh_vertex& v0, const gh_vertex& v1, const gh_vertex& v2) {
    double ds_x = v1.x - v0.x;
    double ds_y = v1.y - v0.y;
    double dc_x = v2.x - v0.x;
    double dc_y = v2.y - v0.y;
                                      
    return (dc_y*ds_x - dc_x*ds_y);
}

inline bool edge_present(vector<gh_vertex>& verts, int other, int current) {

    int& neighbour = verts[current].neighbour;

    if (neighbour >= 0) { // we have a neighbour
        int& nc = verts[current].neighbour;
        if (verts[nc].prev == verts[other].neighbour ||
            verts[nc].next == verts[other].neighbour) {
            
            return true;
        }
    } 
    
    return false;
}

static void set_traversal_flag(vector<gh_vertex>& verts, const Polygon_geom& b, int current,
    double xoffset=0, double yoffset=0) {
    
    const int& prev = verts[current].prev;
    const int& next = verts[current].next;
    
    const int& neigh = verts[current].neighbour;
    const int& nprev = verts[neigh].prev;
    const int& nnext = verts[neigh].next;
    
    
    // classify the prev->current edge
    edge_status prev_status = edge_present(verts, prev, current) ? EDGE_ON : EDGE_OUT;
    if (prev_status != EDGE_ON) {
        switch (
            b.classify( 
                0.5*(verts[prev].x + verts[current].x) - xoffset, 
                0.5*(verts[prev].y + verts[current].y) - yoffset
            ) 
        ) {
        case ON:
            printf("Warning! (prev) This should not be possible\n");
            break;
        case INSIDE:
            prev_status = EDGE_IN;
            break;
        case OUTSIDE:
            prev_status = EDGE_OUT;
        }
    }
    
    // classify the current->next edge
    edge_status next_status = edge_present(verts, next, current) ? EDGE_ON : EDGE_OUT;
    if (next_status != EDGE_ON) {
        switch (
            b.classify( 
                0.5*(verts[next].x + verts[current].x) - xoffset, 
                0.5*(verts[next].y + verts[current].y) - yoffset
            ) 
        ) {
        case ON:
            printf("Warning! (next) This should not be possible\n");
            break;
        case INSIDE:
            next_status = EDGE_IN;
            break;
        case OUTSIDE:
            next_status = EDGE_OUT;
        }
    }
    
    traversal_flag& trav = verts[current].flag;
    
    int combined_edge_status = prev_status * 3 + next_status;
    switch (combined_edge_status) {
    case 8: // on/on
        trav = NONE;
        break;
        
    case 6: // on/out
    case 5: // in/on
    case 3: // in/out
        trav = EX;
        break;
    
    case 7: // on/in
    case 2: // out/on
    case 1: // out/in
        trav = EN;
        break;
        
    case 4: // in/in
        trav = EXEN;
        break;
    
    case 0: // out/out
        trav = ENEX;
        break;
    }
    
    // compute cross_change flag. probably not needed for every vertex ...
    if (trav == EXEN || trav == ENEX) {
        // assume this function is called with C (poly 1), so S is neigh, C is self
        double tr1_cross = triangle_cross(verts[nprev], verts[prev], verts[nnext]);
        double tr2_cross = triangle_cross(verts[nprev], verts[prev], verts[next]);
        verts[current].cross_change = tr1_cross * tr2_cross < 0;
    }
}

void gh_phase_two(vector<gh_vertex>& verts, const Polygon_geom& b, int first_vert_index, double xoffset, double yoffset) {

    // set the traversal flags for C (poly 1)
    int current = first_vert_index;
    do {
        if (verts[current].isect) {
            set_traversal_flag(verts, b, current, xoffset, yoffset);
        }
        
        current = verts[current].next;
    } while (current != first_vert_index); 
    
    // now set the "couple" field
    current = first_vert_index;
    do {
        if ( verts[current].isect && (verts[current].flag == EN || verts[current].flag == EX) ) {
            int next = verts[current].next;
                
            while (verts[next].isect && verts[next].flag != verts[current].flag &&
                   verts[next].couple == -1) {
                next = verts[next].next;
            }
            if (current != next &&
                verts[next].isect && verts[next].flag == verts[current].flag && 
                verts[current].couple == -1 && verts[next].couple == -1) { // only couple uncoupled vertices?
                verts[current].couple = next;
                verts[next].couple = current;
                
                // should we couple the same vertices in the other polygon?
                verts[verts[current].neighbour].couple = verts[next].neighbour;
                verts[verts[next].neighbour].couple = verts[current].neighbour;
            } 
        }
        
        current = verts[current].next;
    } while (current != first_vert_index); 
}

void gh_phase_two_b(vector<gh_vertex>& verts, const Polygon_geom& b, int poly1_start) {

    // since we have already set the flags for C, we can quickly
    // scan its intersections for a suitable candidate intersection
    const int first_vert_index = 0;
    int current = poly1_start;
    int start = -1;
    
    do {
        // look for an intersection with flag that tells us something
        // about the orientation (reverse or not)
        if (verts[current].isect && verts[current].flag != 0) { 
            start = current;
            break;
        }
        current = verts[current].next;
    } while (current != poly1_start); 
    
    if (start != -1) {
        start = verts[start].neighbour;
    } else {
        // it appears that all of the intersections had their
        // traversal flags set to "NONE". this probably means that
        // we can start with any of them 
        current = first_vert_index;
        do {
            if (verts[current].isect) {
                start = current;
                break;
            }
            current = verts[current].next;
        } while (current != first_vert_index); 
        printf("no intersections with flags != NONE found. picking first intersection instead ...\n");
    }

    bool reverse = false;
    if (start >= 0) {
        set_traversal_flag(verts, b, start);
    
        if ( (verts[start].flag == EN && verts[verts[start].neighbour].flag == EX) ||
             (verts[start].flag == EX && verts[verts[start].neighbour].flag == EN) || 
             (verts[start].flag == ENEX && verts[verts[start].neighbour].flag == EXEN) ||
             (verts[start].flag == EXEN && verts[verts[start].neighbour].flag == ENEX)) { 
            reverse = true;
        } else {
            if (verts[start].flag == verts[verts[start].neighbour].flag) {
                reverse = false;
            } else {
                // everything else, assume reverse is true, but we would still
                // like to know about it ...
                printf("flags in undefined state (%d, %d)->%d. help!\n",
                    verts[start].flag, verts[verts[start].neighbour].flag, true);
                reverse = true;
            }
        }
        
        current = start;
        do {
            if (verts[current].isect) {
                if (reverse) {
                    traversal_flag nflag = NONE;
                    switch (verts[verts[current].neighbour].flag) {
                    case NONE:
                        nflag = NONE; // NONE -> NONE, I think ...
                        break;
                    case EN:
                        nflag = EX;
                        break;
                    case EX:
                        nflag = EN;
                        break;
                    case EXEN:
                        nflag = ENEX;
                        break;
                    case ENEX:
                        nflag = EXEN;
                        break;
                    }
                    verts[current].flag = nflag;
                } else {
                    verts[current].flag = verts[verts[current].neighbour].flag;
                }
            }
            current = verts[current].next;
        } while (current != start);
        
    } else {
        printf("no intersections in subject poly. weird.\n");
    }
}

static int proceed(int current, traversal_edge status, vector<gh_vertex>& verts, vector<cv::Vec2d>& vertices) {
    if (verts[current].neighbour != -1) { // TODO: This is a bit extreme? What about ENEX or EXEN nodes?
        verts[current].next_poly = -51;
        verts[verts[current].neighbour].next_poly = -50; 
    }
    
    if (status == D1) {
        // check for degeneracy first
        cv::Vec2d next(verts[current].x, verts[current].y);
        size_t v_size = vertices.size();
        if (v_size > 1) {
            if ( fabs(cross(vertices[v_size-1], vertices[v_size-2], next)) < 1e-11 ) {
                vertices[v_size-1] = next;
                return verts[current].next;
            }
        }
        vertices.push_back(next);
        return verts[current].next;
    }
    if (status == D2) {
        // check for degeneracy first
        cv::Vec2d next(verts[current].x, verts[current].y);
        size_t v_size = vertices.size();
        if (v_size > 1) {
            if ( fabs(cross(vertices[v_size-1], vertices[v_size-2], next)) < 1e-11 ) {
                vertices[v_size-1] = next;
                return verts[current].prev;
            }
        }
        vertices.push_back(next);
        return verts[current].prev;
    }
    // touch both intersections
    verts[current].next_poly = -2;
    verts[verts[current].neighbour].next_poly = -2; 
    
    return verts[current].neighbour;
}

static traversal_edge delete_flag1(int current, traversal_edge status, vector<gh_vertex>& verts) {
    switch (verts[current].flag) {
    case ENEX:
        verts[current].flag = NONE;
        if (verts[current].cross_change) {
            return (status == D3) ? D3 : D4;
        }
        return (status == D3) ? D4 : D3;
        break;
    case EXEN:
        verts[current].flag = (status == D3) ? EN : EX;
        return (status == D3) ? D2 : D1;
        break;
    case EN:
        verts[current].flag = NONE;
        return D1;
        break;
    case EX:
        verts[current].flag = NONE;
        return D2;
        break;
    case NONE: // never reached, parent skips this case
        break;
    }
    //assert(false);
    return D1; // not usually reachable
}

static traversal_edge delete_flag2(int current, int prev, traversal_edge status, vector<gh_vertex>& verts) {
    switch (verts[current].flag) {
    case ENEX:
        verts[current].flag = (status == D1) ? EX : EN;
        if (verts[current].cross_change) {
            return (status == D1) ? D4 : D3;
        }
        return (status == D1) ? D3 : D4;
        break;
    case EXEN:
        verts[current].flag = (status == D1) ? EN : EX;
        if (verts[current].cross_change) {
            return (status == D1) ? D4 : D3;
        }
        return (status == D1) ? D3 : D4;
        break;
    case EN:
        verts[current].flag = NONE;
        if (status == D1 && 
            verts[current].couple != -1 && 
            verts[prev].couple == current) {
            
            // mark both vertices in couple to prevent them from being selected
            verts[current].next_poly = -4;
            verts[verts[current].couple].next_poly = -5; //?
            
            return D1; // keep going in the same direction when encountering a coupled vertex
        }
        return (status == D1) ? D3 : D4;
        break;
    case EX:
        verts[current].flag = NONE;
        if (status != D1 && 
            verts[current].couple != -1 &&
            verts[prev].couple == current) {
            
            // mark both vertices in couple to prevent them from being selected
            verts[current].next_poly = -4;
            verts[verts[current].couple].next_poly = -5; //?
            
            return D2;
        }
        return (status == D1) ? D3 : D4;
        break;
    case NONE: // never reached, parent skips this case
        break;
    }
    //assert(false);
    return D1; // not usually reachable
}

int select_vertex(vector<gh_vertex>& verts, int vs) {
    
    for (int si=0; si < vs; si++) {
        if (verts[si].next_poly == 1 && verts[si].isect && verts[si].flag != NONE) {
        
            // we have an unmarked intersection
            
            if (verts[si].couple != -1) {  // si is part of an unmarked couple
                if (verts[si].flag == EN && verts[verts[si].couple].flag == EN) {
                    return verts[si].couple;
                } else {
                    if (verts[si].flag == EX && verts[verts[si].couple].flag == EX)  {
                        return si;
                    } else {
                        // we reach this state if one of the couple was accessed through the neighbour pointer
                        verts[si].next_poly = -11;
                        continue; // try next si?
                    }
                }
            } 
            
            return si;
        }
    }
    
    return -1;
}

static double compute_area(const vector<cv::Vec2d>& points) {
    double A = 0;
    for (size_t i=0; i < points.size(); i++) {
        int ni = (i+1) % points.size();
        A += points[i][0]*points[ni][1] - points[ni][0]*points[i][1];
    }
    return 0.5 * fabs(A);
}

double gh_phase_three(vector<gh_vertex>& verts, int vs, [[maybe_unused]] int first_isect_index, vector<Polygon_geom>& polys, bool area_only) {
    double combined_area = 0;

    bool done = false;
    while (!done) {
        vector<cv::Vec2d> vertices;
        vertices.reserve(vs);
        
        int current = select_vertex(verts, vs);
        done = current < 0;
        if (!done) {
            verts[current].next_poly = -10;
            
            int prev = current;
            int start = current;
            
            traversal_edge status = delete_flag1(current, D3, verts);
            current = proceed(current, status, verts, vertices);
            
            while (current != start) {
                if (verts[current].flag) {
                    if (status == D1 || status == D2) {
                        status = delete_flag2(current, prev, status, verts); 
                    } else {
                        status = delete_flag1(current, status, verts); 
                    }
                    prev = current;
                }
                current = proceed(current, status, verts, vertices);
            }
            if (vertices.size() > 2) {
                
                if (area_only) {
                    combined_area += compute_area(vertices);
                } else {
            
                    // quickly check for degeneracy with last point wrapping around
                    if ( fabs(cross(vertices[vertices.size()-1], vertices[0], vertices[vertices.size()-2])) < 1e-11 ) {
                        vertices.erase(--vertices.end());
                    }
                    polys.push_back(Polygon_geom(vertices));
                }
            }
            vertices.clear();
        }
        
    }
    
    return combined_area;
}

} // namespace GH_clipping

