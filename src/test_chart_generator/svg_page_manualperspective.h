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
#ifndef SVG_PAGE_MANUALPERSPECTIVE_H
#define SVG_PAGE_MANUALPERSPECTIVE_H

#include "svg_page.h"
#include "include/fiducial_positions.h"

#include <string>

using std::string;
using std::ostringstream;

class Svg_page_manualperspective : public Svg_page {
  public:
    Svg_page_manualperspective(const string& page_spec, const string& fname) 
      : Svg_page(page_spec, fname, 10000), cop(0,0,0), pl(0,0,2000), 
        normal(0, 1.0/sqrt(2.0), -1.0/sqrt(2.0)), fl(50), 
        sensor(15.6, 23.6), outside_limits(false), width_scale(1) {
        
        // must call set_viewing_parameters before render
    }
    
    void render(void) {
        double ang=5.0/180.0*M_PI;
        double bsize=0.005; // block size ....
        
        const double cw = 25;
        const int nc = 9;
        const double cs = 1.1;
        for (int col=0; col <= nc; col++) {
            double offset = -0.75*cw*bsize;
            if (col > nc/2) offset = -offset - 0.6*cw*bsize;
            double cstart = -(nc/2.0+0.25)*(cw+5)*bsize*cs + col*(cw+5)*bsize*cs;
            strip(offset + cstart, -0.8 + (col%2)*5*bsize, cw*bsize, 5*bsize, 35, 1, 
                -ang * ((col > nc/2) ? -1 : 1)
            ); // for perspective blocks
        }
        
        codedcircles(5);
        
        chevrons();
        
        // print out chart specification
        fprintf(fout, "\n");
        text(chart_description, 0.1, 0.97, 30);
        
        if (outside_limits) {
            printf("Warning: chart size exceeds the dimensions of the page you have specified.\n");
            printf("\tThe most likely cause of this is that the specified distance is too small,\n");
            printf("\ti.e., you are too close to the chart; increase distance with -d.\n");
        }
    }
    
    void set_viewing_parameters(double distance_to_target, double angle=-45.0/180*M_PI) {
        pl[2] = distance_to_target;
        normal[1] = cos(angle);
        normal[2] = sin(angle);
        
        double effective_angle = 2*atan(width_mm/(2*distance_to_target));
        
        width_scale = fl * tan(effective_angle/2);
        double f = 15.6/(2*tan(effective_angle/2));
        printf("maximum focal length (15.6 mm sensor height): %.2lf mm\n", f);
        printf("maximum focal length (24 mm sensor height): %.2lf mm\n", 24/(2*tan(effective_angle/2)));
        printf("distance-to-focal-length ratio: %lfx\n", pl[2]/f);
        
        ostringstream ss;
        ss << "Chart size: " << page_size << ", ";
        ss << "design viewing distance: " << pl[2] << " mm, ";
        ss << "design viewing angle: " << lrint(angle/M_PI*1800)/10.0 << " degrees, ";
        ss << "maximum focal length (sensor height = " << 15.6 << " mm): " << lrint(f*10)/10.0 << " mm, ";
        ss << "ratio at maximum focal length: " << lrint(pl[2]/f*10)/10.0;
        chart_description = ss.str();
    }
    
  protected:  
  
    dPoint project_core(double x, double y) {
        Vec3d d(x*width_scale, y*width_scale, fl);
        
        d = d * (1.0/norm(d));
        double t = ((pl - cop).ddot(normal)) / (d.dot(normal));
        Vec3d pi = cop + d * t;
        
        Vec3d cr = normal.cross(Vec3d(1,0,0)); // basis vector in plane, along y direction
        dPoint p;
        p.y = (pi - pl).dot(cr);
        p.x = (pi - pl).dot(Vec3d(1,0,0));
        
        return p;
    }
  
    virtual iPoint project(double x, double y, [[maybe_unused]] bool silent=false) {
    
        dPoint p = project_core(x,y);
        
        p.x = floor(p.x*sscale + 0.5*width);
        p.y = floor(p.y*sscale + 0.5*height);
        
        if (p.x < 0 || p.x > width ||
            p.y < 0 || p.y > height) {
            
            outside_limits = true;
        }
        
        return iPoint(int(p.x), int(p.y));
    }
    
    
    void perspective_rectangle(double tlx, double tly, double width, double height, double angle, bool right=false) {
    
        vector<iPoint> points;
        if (right) {
            points.push_back(project(tlx, tly));
            points.push_back(project(tlx - cos(angle)*width, tly + sin(angle)*height));
            points.push_back(project(tlx - sin(angle)*width - cos(angle)*width, tly - cos(angle)*height + sin(angle)*height));
            points.push_back(project(tlx - sin(angle)*width, tly - cos(angle)*height));
        } else {
            points.push_back(project(tlx, tly));
            points.push_back(project(tlx + cos(angle)*width, tly + sin(angle)*height));
            points.push_back(project(tlx + sin(angle)*width + cos(angle)*width, tly - cos(angle)*height + sin(angle)*height));
            points.push_back(project(tlx + sin(angle)*width, tly - cos(angle)*height));
        }
        
        if (collision_free(points)) {
            fprintf(fout, "  <polygon points=\"%d,%d ", points[0].x, points[0].y);
            fprintf(fout, "%d,%d ", points[1].x, points[1].y);
            fprintf(fout, "%d,%d ", points[2].x, points[2].y);
            fprintf(fout, "%d,%d\" style=\"%s\"/>", points[3].x, points[3].y, style.c_str());
        }
    }
    
    void strip(double tlx, double tly, double swidth, double sheight, size_t nrows, size_t ncols, double ang) {
        
        for (size_t ry=0; ry < nrows; ry++) {
            for (size_t rx=0; rx < ncols; rx++) {
                double ypos = ry * 1.8 * sheight + tly; // perspective
                double xpos = tlx - rx * 2.0 * swidth;
          
                perspective_rectangle(xpos, ypos, swidth, sheight, -ang);
                fprintf(fout, "\n");
      
            } // columns
        } // rows
    }
    
    bool collision_free(vector<iPoint>& points) {
        bool clean = true;
        const double rad = 7;
        for (int i=0; i < n_fiducials && clean; i++) { // defined in "include/fiducial_positions.h"
            iPoint rtl(points[0]);
            iPoint rbr(points[0]);
            for (size_t j=1; j < points.size(); j++) {
                rtl.x = min(rtl.x, points[j].x);
                rtl.y = min(rtl.y, points[j].y);
                rbr.x = max(rbr.x, points[j].x);
                rbr.y = max(rbr.y, points[j].y);
            }
            for (size_t j=0; j < points.size() && clean; j++) {
                iPoint tl = scale(main_fiducials[i].rcoords.x - rad, main_fiducials[i].rcoords.y - rad*0.9);
                iPoint br = scale(main_fiducials[i].rcoords.x + rad, main_fiducials[i].rcoords.y + rad*0.9);
                
                bool lclean=false;
                
                if ((br.x) < rtl.x ||    // *this is left of b
                    (tl.x) > rbr.x) {    // *this is right of b
                      
                    lclean = true;
                }
                  
                if ((br.y) < rtl.y ||  // *this is above b
                    (tl.y) > rbr.y) {  // *this is below b
                  
                    lclean |= true;
                }
                
                // crop it to one side. but which side?
                int nminx = rbr.x;
                int nmaxx = rtl.x;
                if (!lclean) {
                    // is the circle centroid left or right of the rectangle centroid?
                    int midx = (tl.x + br.x)/2;
                    int rmidx = (rtl.x + rbr.x)/2;
                    if (midx < rmidx) { // circle is left
                        for (size_t k=0; k < points.size(); k++) {
                            points[k].x = max(points[k].x, br.x);
                            nminx = min(nminx, points[k].x);
                            nmaxx = max(nmaxx, points[k].x);
                        }
                    } else { // circle is right
                        for (size_t k=0; k < points.size(); k++) {
                            points[k].x = min(points[k].x, tl.x);
                            nminx = min(nminx, points[k].x);
                            nmaxx = max(nmaxx, points[k].x);
                        }
                    }
                    
                    // drop it if too small...
                    if ( fabs(double(nmaxx - nminx)) < 0.2*fabs(double(rbr.x - rtl.x)) ) {
                        clean = false;
                    }
                }
            }
        }
        return clean;
    }
    
    void codedcircles(double swidth) {
        
        for (int i=0; i < n_fiducials; i++) { // defined in "include/fiducial_positions.h"
            sector_circle(main_fiducials[i].rcoords.x, main_fiducials[i].rcoords.y, swidth, main_fiducials[i].code);
            fprintf(fout, "\n");
        }
        
        const double shift = -2;
        
        // strip of small boxes
        for (int y=-80; y < 80; y++) {
            iPoint c1 = scale(swidth-1 + shift, 2*y);
            iPoint c2 = scale(swidth+1 + shift, 2*y+1);
            fprintf(fout, "\n");
            rect(c1.x, c1.y, c2.x - c1.x, c2.y - c1.y);
        }
        
        // focusing aids
        iPoint c1 = scale(0 + shift - 1.5, -0.25);
        iPoint c2 = scale(+1 + shift - 1.5,  0.25);
        fprintf(fout, "\n");
        rect(c1.x, c1.y, c2.x - c1.x, c2.y - c1.y);
        c1 = scale(-1 + shift + 1.5, -0.25);
        c2 = scale(0 + shift + 1.5,  0.25);
        fprintf(fout, "\n");
        rect(c1.x, c1.y, c2.x - c1.x, c2.y - c1.y);
    }
    
    void direct_rectangle(double cx, double cy, double width, double height, double angle) {
            fprintf(fout, "  <polygon points=\"");
            for (int i=0; i < 4; i++) {
                double theta = 2*M_PI*i/double(4) + M_PI/4.0; 
                double urx = 0.5*width*cos(theta);
                double ury = 0.5*height*sin(theta);
                double rx = urx*cos(angle) - ury*sin(angle);
                double ry = urx*sin(angle) + ury*cos(angle);
                iPoint p = scale(cx + rx, cy + ry);
                fprintf(fout, "%d,%d ", p.x, p.y);
            }
            fprintf(fout, "\" style=\"%s\"/>", style.c_str());
          
    }
    
    void chevrons(void) {
        const double shift = 0;
        const int fine_rows = 320;
        for (int row=0; row < fine_rows; row++) {
            direct_rectangle(-7 + shift, row*1 - fine_rows/2, 3.5, 0.5, 45.0/180*M_PI);
            fprintf(fout, "\n");
        }
        for (int row=0; row < fine_rows/4; row++) {
            for (int col=0; col < 3; col ++) {
                direct_rectangle(col - 7 + shift, row*4 - fine_rows/2, 3.5, 0.5, -45.0/180*M_PI);
                fprintf(fout, "\n");
            }
        }
    }
    
    Vec3d cop;      // centre of projection
    Vec3d pl;       // position of focal plane
    Vec3d normal;     // normal vector of target
    double fl;          // focal length of lens
    Vec2d sensor;
    bool outside_limits;
    double width_scale;
    std::string chart_description;
};

#endif
