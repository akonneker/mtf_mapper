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
#ifndef SVG_PAGE_PERSPECTIVE_H
#define SVG_PAGE_PERSPECTIVE_H

#include "svg_page.h"

#include <string>

using std::string;
using std::ostringstream;

class Svg_page_perspective : public Svg_page {
  public:
    Svg_page_perspective(const string& page_spec, const string& fname) 
      : Svg_page(page_spec, fname), cop(0,0,0), pl(0,0,2000), 
        normal(0, 1.0/sqrt(2.0), -1.0/sqrt(2.0)), fl(50), 
        sensor(15.6, 23.6), outside_limits(false), width_scale(1) {
        
        // must call set_viewing_parameters before render
    }
    
    void render(void) {
        double ang=10.0/180.0*M_PI;
        double bsize=0.0214; // block size ....
        
        // TODO: work out the correct number of blocks, given the block size
        
        strip(-6*2*1.5*2*bsize, -30*bsize, 2.5*bsize, 2*bsize, 15, 5, -ang, 1);
        strip( 6*2*1.5*2*bsize, -30*bsize, 2.5*bsize, 2*bsize, 15, 5, -ang, 0);
        perspective_rectangle(-0.2, 0.5, 0.4, 0.5, 0*M_PI/180.0, 0);
        
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
        
        p.x = floor(p.x*100 + 0.5*width);
        p.y = floor(p.y*100 + 0.5*height);
        
        if (p.x < 0 || p.x > width ||
            p.y < 0 || p.y > height) {
            
            outside_limits = true;
        }
        
        return iPoint(int(p.x), int(p.y));
    }
    
    
    
    void perspective_rectangle(double tlx, double tly, double width, double height, double angle, bool right) {
        if (right) {
            iPoint p = project(tlx, tly);
            fprintf(fout, "  <polygon points=\"%d,%d ", p.x, p.y);
            p = project(tlx - cos(angle)*width, tly + sin(angle)*height);
            fprintf(fout, "%d,%d ", p.x, p.y);
            p = project(tlx - sin(angle)*width - cos(angle)*width, tly - cos(angle)*height + sin(angle)*height);
            fprintf(fout, "%d,%d ", p.x, p.y);
            p = project(tlx - sin(angle)*width, tly - cos(angle)*height);
            fprintf(fout, "%d,%d\" style=\"%s\"/>", p.x, p.y, style.c_str());
          
        } else {
            iPoint p = project(tlx, tly);
            fprintf(fout, "  <polygon points=\"%d,%d ", p.x, p.y);
            p = project(tlx + cos(angle)*width, tly + sin(angle)*height);
            fprintf(fout, "%d,%d ", p.x, p.y);
            p = project(tlx + sin(angle)*width + cos(angle)*width, tly - cos(angle)*height + sin(angle)*height);
            fprintf(fout, "%d,%d ", p.x, p.y);
            p = project(tlx + sin(angle)*width, tly - cos(angle)*height);
            fprintf(fout, "%d,%d\" style=\"%s\"/>", p.x, p.y, style.c_str());
          
        }
    }
    
    void strip(double tlx, double tly, double swidth, double sheight, size_t nrows, size_t ncols, double ang, bool left) {
  
        for (size_t ry=0; ry < nrows; ry++) {
            for (size_t rx=0; rx < ncols; rx++) {
      
                if (left) {
                    double ypos = ry * 2.0 * sheight + tly;
                    double xpos = rx * 2.0 * swidth  + tlx;
              
              
                    perspective_rectangle(xpos, ypos, swidth*0.8, sheight*0.8, ang, 0);
                    fprintf(fout, "\n");
                    if (ry < nrows-1) {
                        xpos += swidth;
                        ypos += sheight;
                        perspective_rectangle(xpos, ypos, swidth*0.8, sheight*0.8, ang, 0);
                        fprintf(fout, "\n");
                    }
                } else {
                    double ypos = ry * 2.0 * sheight + tly;
                    double xpos = tlx - rx * 2.0 * swidth;
              
                    perspective_rectangle(xpos, ypos, swidth*0.8, sheight*0.8, ang, 1);
                    fprintf(fout, "\n");
                    if (ry < nrows-1) {
                        xpos -= swidth;
                        ypos += sheight;
                        perspective_rectangle(xpos, ypos, swidth*0.8, sheight*0.8, ang, 1);
                        fprintf(fout, "\n");
                    }
                }
      
            } // columns
        } // rows
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
