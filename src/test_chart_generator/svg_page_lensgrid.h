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
#ifndef SVG_PAGE_LENSGRID_H
#define SVG_PAGE_LENSGRID_H

#include "svg_page.h"
#include "include/mtf50_edge_quality_rating.h"
#include "lightweight_polygon.h"

#include <vector>
#include <cmath>
using std::vector;
using std::pair;
using std::make_pair;

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

class Svg_page_lensgrid : public Svg_page {
  public:
    Svg_page_lensgrid(const string& page_spec, const string& fname) 
    : Svg_page(page_spec, fname), centre(0.5, 0.5*double(height_mm)/double(width_mm)),
      off(0), scale(10.0/100) {
      
        fiducial_scale_index = fiducial_mapping_index(int(fiducial_scale_index) + int(fiducial_mapping_index::A0L));
    }
    
    void render(void) {
        vector<Point2f> centers = coded_sector_circles(10*fiducial_scale);
        grid(0.04, 17, (width_mm == height_mm) ? 17 : 11, centers); 
    }
    
  protected:
    virtual iPoint project(double x, double y, [[maybe_unused]] bool silent=false) {
        
        x = floor(x*width);
        y = floor(y*width);
        
        return iPoint(int(x), int(y));
    }

    double quant_angle(double x) {
        double quad1 = fabs(fmod(x, M_PI/2.0));
        if (quad1 > M_PI/4.0) {
            quad1 = M_PI/2.0 - quad1;
        }
        quad1 = quad1 / M_PI * 180;
        return quad1;
    }
    
    void triangle(dPoint a, dPoint b, dPoint c) {
       
        iPoint p = project(a.x, a.y);
        fprintf(fout, "  <polygon points=\"%d,%d ", p.x, p.y);
        p = project(b.x, b.y);
        fprintf(fout, "%d,%d ", p.x, p.y);
        p = project(c.x, c.y);
        fprintf(fout, "%d,%d\" style=\"%s\"/>\n", p.x, p.y, style.c_str());
    }
    
    bool place_trapezoid(double xpos, double ypos, double width, double angle_thresh, const vector<Point2f>& circles) {
        vector<iPoint> coords(4);
        
        static double direction = 1;
        
        //direction *= -1;
        
        double ewidth = 0;
        
        bool overlap;
        bool done = false;
        int tries = 0;
        int last_change = 0;
        do {
            dPoint delta = dPoint(xpos, ypos) - centre;
            
            //#define SQUARES
            
            double norm = sqrt(delta.x*delta.x + delta.y*delta.y);
            delta.x /= norm;
            delta.y /= norm;
            
            dPoint omx = dPoint(xpos + delta.x*width*sqrt(0.5), ypos + delta.y*width*sqrt(0.5));
            dPoint tang = dPoint(-delta.y, delta.x);
            
            dPoint tl(omx.x - 0.5*width*tang.x, omx.y - 0.5*width*tang.y);
            iPoint p = project(tl.x, tl.y);
            coords[0] = p;
            
            dPoint rad = tl - centre;
            norm = sqrt(rad.x*rad.x + rad.y*rad.y);
            rad.x /= norm;
            rad.y /= norm;
            double radang1 = atan2(fabs(rad.y), fabs(rad.x)) * 180.0/M_PI;
            if (radang1 > 45) radang1 = 90 - radang1;
            #ifdef SQUARES
            rad = delta;
            #endif
            
            dPoint bl(tl.x - rad.x*width + width*ewidth*tang.x, tl.y - rad.y*width + width*ewidth*tang.y);
            p = project(bl.x, bl.y);
            coords[1] = p;
            
            dPoint tr(omx.x + 0.5*width*tang.x, omx.y + 0.5*width*tang.y);
            
            rad = tr - centre;
            norm = sqrt(rad.x*rad.x + rad.y*rad.y);
            rad.x /= norm;
            rad.y /= norm;
            double radang2 = atan2(fabs(rad.y), fabs(rad.x)) * 180.0/M_PI;
            if (radang2 > 45) radang2 = 90 - radang2;
            #ifdef SQUARES
            rad = delta;
            #endif
            
            dPoint br(tr.x - rad.x*width - width*ewidth*tang.x, tr.y - rad.y*width - width*ewidth*tang.y);
            p = project(br.x, br.y);
            coords[2] = p;
            
            p = project(tr.x, tr.y);
            coords[3] = p;
            
            double angles[4];
            for (int i=0; i < 4; i++) {
                int next = (i + 1) % 4;
                double angle = atan2(fabs(double(coords[i].y - coords[next].y)), fabs(double(coords[i].x - coords[next].x)));
                angle *= 180.0/M_PI;
                if (angle > 45) {
                    angle = 90 - angle;
                }
                // identify which edge this is
                int j = -1;
                double l2 = sqrt(SQR(coords[i].x - coords[next].x) + SQR(coords[i].y - coords[next].y));
                double dot = (coords[i].x - coords[next].x)*delta.x + (coords[i].y - coords[next].y)*delta.y;
                dot /= l2;
                if (fabs(dot) < 0.707) {
                    j = (dot < 0) ? 0 : 1;
                } else {
                    j = (dot < 0) ? 2 : 3;
                }
                if (j < 0) {
                    printf("error. not a valid orientation\n");
                    exit(1);
                }
                angles[j] = angle;
            }
            
            if (tries % 2 == 0) {
                if (angles[0] < angle_thresh || angles[0] > 43.5 || fabs(angles[0] - 26.565) < 2 ||
                    angles[1] < angle_thresh || angles[1] > 43.5 || fabs(angles[1] - 26.565) < 2 ) { 
                    xpos -= 0.001*delta.y*direction;
                    ypos += 0.001*delta.x*direction;
                    last_change = tries;
                }
                
            } else {
                if (angles[2] < 2 || angles[2] > 43 || angles[3] < 2 || angles[3] > 43 ||
                    fabs(angles[2] - 26.565) < 2 || fabs(angles[3] - 26.565) < 2) {
                    ewidth += 0.01;
                    last_change = tries;
                }
            }
            
            overlap = false;
            double mindist = 1e50;
            for (int i=0; i < 4; i++) {
                Point2f fp = Point2f(coords[i].x, coords[i].y);
                Point2f fpn = Point2f(coords[(i+1)%4].x, coords[(i+1)%4].y);
                for (double t=0; t <= 1.0; t += 0.1) {
                    Point2f lp = t*fp + (1-t)*fpn;
                    for (size_t j=0; j < circles.size(); j++) {
                        double dist = cv::norm(lp - circles[j]);
                        mindist = std::min(dist, mindist);
                    }
                }
            }
            
            if (mindist < 7.5*fiducial_scale*sscale) {
                overlap = true;
            }
            
            vector<cv::Vec2d> verts;
            double pc_x = 0;
            double pc_y = 0;
            for (int i=0; i < 4; i++) {
                pc_x += coords[i].x*0.25;
                pc_y += coords[i].y*0.25;
            }
            for (int i=0; i < 4; i++) {
                verts.push_back(
                    cv::Vec2d(
                        coords[3-i].x + 0.5*(coords[3-i].x - pc_x), 
                        coords[3-i].y + 0.5*(coords[3-i].y - pc_y)
                    )
                );
            }
            Convex_polygon pg(verts);
            
            bool all_clear = true;
            for (size_t pi=0; all_clear && pi < existing_trapezoids.size(); pi++) {
                if (pg.intersection_area(existing_trapezoids[pi]) > 0) {
                    all_clear = false;
                }
            }
            
            if (!all_clear) {
                overlap = true;
            }
            
            tries++;
            if (tries > last_change+3) done = true;
            
            // estimate radial error
            //printf("%lf %lf\n", 
            //    std::min(fabs(radang1-angles[2]), fabs(radang2-angles[2])),
            //    std::min(fabs(radang1-angles[3]), fabs(radang2-angles[3]))
            //);
            
        } while (!done);
        
        
        if (!overlap) {
            fprintf(fout, "  <polygon points=\"%d,%d ", coords[0].x, coords[0].y);
            //fprintf(stderr, "4\n%lf %lf\n", coords[0].y*scale + off, coords[0].x*scale + off);
            fprintf(fout, "%d,%d ", coords[1].x, coords[1].y);
            //fprintf(stderr, "%lf %lf\n", coords[1].y*scale + off, coords[1].x*scale + off);
            fprintf(fout, "%d,%d ", coords[2].x, coords[2].y);
            //fprintf(stderr, "%lf %lf\n", coords[2].y*scale + off, coords[2].x*scale + off);
            fprintf(fout, "%d,%d\" style=\"%s\"/>\n", coords[3].x, coords[3].y, style.c_str());
            //fprintf(stderr, "%lf %lf\n", coords[3].y*scale + off, coords[3].x*scale + off);
            
            vector<cv::Vec2d> verts;
            for (int i=0; i < 4; i++) {
                verts.push_back(cv::Vec2d(coords[3-i].x, coords[3-i].y));
            }
            existing_trapezoids.push_back(verts);
        }
        
        return overlap;
    }
    
    vector<Point2f> coded_sector_circles(double swidth) {
        vector<Point2f> centers;
        
        for (int i=0; i < n_fiducials; i++) { // defined in "include/fiducial_positions.h"
            double ypos = main_fiducials[i].rcoords.y*fiducial_scale*fiducial_position_scale_factor[fiducial_scale_index].second;
            double xpos = main_fiducials[i].rcoords.x*fiducial_scale*fiducial_position_scale_factor[fiducial_scale_index].first;
            sector_circle(xpos, ypos, swidth*0.5, fiducial_code_mapping[fiducial_scale_index][i]); 
            Point2f p(xpos*sscale+ 0.5*width, ypos*sscale + 0.5*height);
            centers.push_back(p);
        }
        return centers;
    }

    void grid(double swidth, size_t nrows, size_t ncols, const vector<Point2f>& circles) {
        const double aspect = double(height_mm) / double(width_mm);
        
        /*
        // If the chart is square, place an invisible rectangle down the centre
        // this forces the trapesiums to avoid this central strip
        // allowing the chart to be split in two.
        if (height_mm == width_mm) {
            vector<cv::Vec2d> verts;
            const double sliver = 0.005;
            iPoint p = project(0.5 - 2*sliver, 0.0);
            verts.push_back(cv::Vec2d(p.x, p.y));
            p = project(0.5 + sliver, 0.0);
            verts.push_back(cv::Vec2d(p.x, p.y));
            p = project(0.5 + sliver, aspect);
            verts.push_back(cv::Vec2d(p.x, p.y));
            p = project(0.5 - 2*sliver, 0.0);
            verts.push_back(cv::Vec2d(p.x, p.y));
            existing_trapezoids.push_back(verts);
        }
        */
        
        const double phi = 1.9/180.0*M_PI;
        for (int row=0; row < (int)nrows; row++) {
            for (int col=0; col < (int)ncols; col++) {
                if (row == (int)nrows/2 && col == (int)ncols/2) continue;
                double sx = (col-int(ncols)/2)/double(ncols)*0.93;
                double sy = aspect*(row-int(nrows)/2)/double(nrows)*0.93;
                double rx = sx*cos(phi) - sy*sin(phi);
                double ry = sx*sin(phi) + sy*cos(phi);
                bool overlapped = place_trapezoid(
                    0.5 + rx,
                    0.5*aspect + ry,
                    swidth*sqrt(0.5)*aspect,
                    (col == (int)ncols/2) || (row == (int)nrows/2) ? 3 : 1.5,
                    circles
                );
                
                if (overlapped) {
                    double offset = 3.0;
                    if (abs(row - (int)nrows/2) <= 1 || abs(col - (int)ncols/2) <= 1) {
                        offset = 15.0;
                    }
                    for (double dp=0; overlapped && dp < offset/180.0*M_PI; dp += 0.25/180.0*M_PI) {
                        rx = sx*cos(phi + dp) - sy*sin(phi + dp);
                        ry = sx*sin(phi + dp) + sy*cos(phi + dp);
                        
                        overlapped = place_trapezoid(
                            0.5 + rx,
                            0.5*aspect + ry,
                            swidth*sqrt(0.5)*aspect,
                            (col == (int)ncols/2) || (row == (int)nrows/2) ? 3 : 1.5,
                            circles
                        );
                        
                        if (overlapped) {
                            rx = sx*cos(phi - dp) - sy*sin(phi - dp);
                            ry = sx*sin(phi - dp) + sy*cos(phi - dp);
                            
                            overlapped = place_trapezoid(
                                0.5 + rx,
                                0.5*aspect + ry,
                                swidth*sqrt(0.5)*aspect,
                                (col == (int)ncols/2) || (row == (int)nrows/2) ? 3 : 1.5,
                                circles
                            );
                        }
                    }
                }
            }
        }
        
        
        triangle(
          dPoint(centre.x, centre.y),
          dPoint(centre.x - 0.04, centre.y + 0.04),
          dPoint(centre.x - 0.04, centre.y - 0.04)
        );
        
        
        triangle(
          dPoint(centre.x, centre.y),
          dPoint(centre.x + 0.04, centre.y + 0.04),
          dPoint(centre.x + 0.04, centre.y - 0.04)
        );
        
        triangle(
          dPoint(1, 0),
          dPoint(1, 0.05),
          dPoint(0.5,  0)
        );
        
        triangle(
          dPoint(1-0.04, 0.05 - 0.04*0.1),
          dPoint(1, 0.05),
          dPoint(1,  0.1)
        );
        
        triangle(
          dPoint(0, aspect),
          dPoint(0, aspect-0.05),
          dPoint(0.5,  aspect)
        );
        
        triangle(
          dPoint(0.04, aspect - 0.05 + 0.04*0.1),
          dPoint(0,  aspect - 0.05),
          dPoint(0,  aspect - 0.1)
        );
        
        triangle(
          dPoint(0.05, 0),
          dPoint(0, 0.5),
          dPoint(0, 0.0)
        );
        
        triangle(
          dPoint(0.05, 0),
          dPoint(0.1, 0),
          dPoint(0.05 - 0.04*0.1, 0.04)
        );
        
        triangle(
          dPoint(1 - 0.05, aspect),
          dPoint(1, aspect - 0.5),
          dPoint(1, aspect)
        );
        
        triangle(
          dPoint(1 - 0.05, aspect),
          dPoint(1 - 0.1, aspect),
          dPoint(1 - 0.05 + 0.04*0.1, aspect - 0.04)
        );
    }
    
    dPoint centre;
    double off;
    double scale;
    
    vector<Convex_polygon> existing_trapezoids;
};

#endif
