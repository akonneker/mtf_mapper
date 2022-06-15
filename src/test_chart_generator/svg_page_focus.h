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
#ifndef SVG_PAGE_FOCUS_H
#define SVG_PAGE_FOCUS_H

#include "svg_page.h"

#include <string>

#include "include/fiducial_positions.h"
#include "alpa_correction.h"

using std::string;
using std::ostringstream;

class Svg_page_focus : public Svg_page {
  public:
    Svg_page_focus(const string& page_spec, const string& fname, bool alpa_scale=false) 
      : Svg_page(page_spec, fname, 10000), cop(0,0,0), pl(0,0,2000), 
        normal(0, 1.0/sqrt(2.0), -1.0/sqrt(2.0)), fl(50), 
        sensor(15.6, 23.6), width_scale(1), clipid(1), alpa_scale(alpa_scale) {
        
        // must call set_viewing_parameters before render
    }
    
    void render(void) {
        double ang=5.0/180.0*M_PI;
        double bsize=0.005; // block size ....
        
        strip(-5*bsize, -0.8, 45*bsize, 5*bsize, 44, 1, -ang); // for perspective blocks
        
        coded_sector_circles(10*fiducial_scale);
        
        chevrons();
        
        if (alpa_scale) {
            alpa_corrections();
        }
        print_checks();
        
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
  
    virtual iPoint project(double x, double y, bool silent=false) {
    
        dPoint p = project_core(x,y);
        
        p.x = floor(p.x*sscale + 0.5*width);
        p.y = floor(p.y*sscale + 0.5*height);
        
        if (p.x < 0 || p.x > width ||
            p.y < 0 || p.y > height) {
            
            outside_limits = !silent;
        }
        
        return iPoint(int(p.x), int(p.y));
    }
    
    
    void perspective_rectangle(double cx, double cy, double width, double height, double angle, [[maybe_unused]] bool print_centroid=false) {
        fprintf(fout, "  <polygon points=\"");
        for (int i=0; i < 4; i++) {
            double theta = 2*M_PI*i/double(4) + M_PI/4.0; 
            double urx = 0.5*width*cos(theta);
            double ury = 0.5*height*sin(theta);
            double rx = urx*cos(angle) - ury*sin(angle);
            double ry = urx*sin(angle) + ury*cos(angle);
            iPoint p = project(cx + rx, cy + ry);
            fprintf(fout, "%d,%d ", p.x, p.y);
            
            
        }
        fprintf(fout, "\" style=\"%s\"/>", style.c_str());
        #if 0
        if (print_centroid) {
            iPoint p2 = project(cx, cy);
            fprintf(stderr, "{%lg,%lg}, ", double(p2.x)/sscale, double(p2.y)/sscale);
        }
        #endif
    }
    
    void strip(double tlx, double tly, double swidth, double sheight, size_t nrows, size_t ncols, double ang) {
        
        for (size_t ry=0; ry < nrows; ry++) {
            for (size_t rx=0; rx < ncols; rx++) {
      
            
                double ypos = ry * 1.4 * sheight + tly; // perspective
                double xpos = tlx - rx * 2.0 * swidth;
          
                perspective_rectangle(xpos, ypos, swidth, sheight, -ang, true);
                fprintf(fout, "\n");
      
            } // columns
        } // rows
    }
    
    void perspective_poly(double cx, double cy, double rad, int sides=64) {
            fprintf(fout, "  <polygon points=\"");
            for (int i=0; i < sides; i++) {
                double theta = 2*M_PI*i/double(sides+1);
                iPoint p = project(cx + rad*cos(theta), cy + rad*sin(theta));
                fprintf(fout, "%d,%d ", p.x, p.y);
            }
            fprintf(fout, "\" style=\"%s\"/>", style.c_str());
          
    }
    
    void coded_sector_circles(double swidth) {
        const double hshift = 20 * fiducial_scale;
        
        #if 0
        
        const double cdist = 30 * fiducial_scale;
        sector_circle(cdist, 0, swidth/2.0, 0);
        fprintf(fout, "\n");
        printf("{0, 0, %lf, %lf, 0, 0},\n", cdist, 0.0);
        sector_circle(-cdist, 0, swidth/2.0, 0);
        fprintf(fout, "\n");
        printf("{0, 0, %lf, %lf, 0, 0},\n", -cdist, 0.0);
        
        const double theta = 35/180.0*M_PI;
        const double radius = 90*sqrt(0.5);
        int code = 2;
        for (double phi=0; phi > -2*M_PI; phi -= M_PI/2.0) {
            double cx = radius*cos(phi + theta);
            double cy = radius*sin(phi + theta);
            
            sector_circle(cx, cy, swidth/2.0, code);
            fprintf(fout, "\n");
            
            int quad = ((cx < 0) ? 1 : 0) | ((cy > 0) ? 2 : 0);
            if (quad <= 1) { // ugly, but swap quadrants 0 and 1
                quad = 1 - quad;
            }
            
            printf("{0, 0, %lf, %lf, %d, %d}%c\n", cx, cy, code, quad, quad == 3 ? ' ' : ','); // don't actually need quadrants for these four
            
            quad++;
            code += 2;
        }
        
        vector<double> base_rad{100, 120, 130, 140};
        vector<double> offset_angle{-35, 20, 70, 45};
        for (size_t i=0; i < base_rad.size(); i++) {
        
            double rad=0;
            for (double phi=0; phi > -2*M_PI; phi -= M_PI/2.0) {
                double cx = (base_rad[i]+rad)*cos(phi + offset_angle[i]/180.0*M_PI);
                double cy = (base_rad[i]+rad)*sin(phi + offset_angle[i]/180.0*M_PI);
                
                sector_circle(cx, cy, swidth/2.0, code);
                
                int quad = ((cx < 0) ? 1 : 0) | ((cy > 0) ? 2 : 0);
                if (quad <= 1) { // ugly, but swap quadrants 0 and 1
                    quad = 1 - quad;
                }
                
                printf("{0, 0, %lf, %lf, %d, %d}%c\n", cx, cy, code, quad, quad == 3 ? ' ' : ','); // don't actually need quadrants for these four
                
                quad++;
                rad += 2.5;
            }
            
            code += 2;
        }
        
        #else 
        for (int i=0; i < n_fiducials; i++) { // defined in "include/fiducial_positions.h"
            sector_circle(main_fiducials[i].rcoords.x*fiducial_scale, main_fiducials[i].rcoords.y*fiducial_scale, swidth*0.5, fiducial_code_mapping[fiducial_scale_index][i]);
            printf("fiducial %d (%lf, %lf) assigned code %d\n", 
                i, main_fiducials[i].rcoords.x, main_fiducials[i].rcoords.y,
                fiducial_code_mapping[fiducial_scale_index][i]
            );
            fprintf(fout, "\n");
        }
        #endif
        
        iPoint vprobe = project(1.07,1, true);
        const int fine_rows = floor(vprobe.x/sscale);
        
        
        for (int y=-fine_rows/(4*fiducial_scale); y < fine_rows/(4*fiducial_scale); y++) {
            iPoint c1 = scale(swidth-fiducial_scale + hshift/2, fiducial_scale*(2*y));
            iPoint c2 = scale(swidth+fiducial_scale + hshift/2, fiducial_scale*(2*y+1));
            fprintf(fout, "\n");
            rect(c1.x, c1.y, c2.x - c1.x, c2.y - c1.y);
        }
        
        // these are probably the horizontal lines ....
        iPoint hprobe = project(0,-1, true);
        const double char_scale = floor(hprobe.x/sscale/sqrt(2.0));
        printf("horizontal scale: %lf\n", char_scale);
        const double rwidth = 2;
        for (int x=0; x < int(char_scale/4); x++) {
            iPoint c1 = scale( (4*x-1)*rwidth*0.5 + 2*hshift, -0.25);
            iPoint c2 = scale( (4*x+1)*rwidth*0.5 + 2*hshift,  0.25);
            fprintf(fout, "\n");
            rect(c1.x, c1.y, c2.x - c1.x, c2.y - c1.y);
        }
        for (int x=0; x < int(char_scale/4); x++) {
            iPoint c1 = scale( (4*(-x)-1)*rwidth*0.5 - 2*hshift, -0.25);
            iPoint c2 = scale( (4*(-x)+1)*rwidth*0.5 - 2*hshift, 0.25);
            fprintf(fout, "\n");
            rect(c1.x, c1.y, c2.x - c1.x, c2.y - c1.y);
        }
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
        iPoint proj = project(1.07,1, true);
        const int fine_rows = floor(proj.x/sscale/fiducial_scale);
        const double hshift = 20;
        for (int row=0; row < fine_rows; row++) {
            direct_rectangle((-7 + hshift)*fiducial_scale, (row*1 - fine_rows/2)*fiducial_scale, 3.5*fiducial_scale, 0.5*fiducial_scale, 45.0/180*M_PI);
            fprintf(fout, "\n");
        }
        for (int row=0; row < fine_rows/4; row++) {
            for (int col=0; col < 3; col ++) {
                direct_rectangle((col - 7 + hshift)*fiducial_scale, (row*4 - fine_rows/2)*fiducial_scale, 3.5*fiducial_scale, 0.5*fiducial_scale, -45.0/180*M_PI);
                fprintf(fout, "\n");
            }
        }
    }
    
    void alpa_corrections(void) {
        bool first = true;
        for (size_t i=0; i < alpa_correction_table.size(); i++) {
            
            double pos  = alpa_correction_table[i].second * (first ? -1 : 1);
            double shim = alpa_correction_table[i].first * (first ? 1 : -1);
            
            if (fabs(pos) < 1e-6) {
                first = false;
            } else {
                direct_rectangle(-25, pos, 10, 0.25, 0);
                iPoint p = scale(-28, pos + (!first ? -1 : 1.75)*1.5);
                fprintf(fout, "  <text x=\"%d\" y=\"%d\" font-family=\"Verdana\" font-size=\"%dmm\" fill=\"black\" > %c%d </text>\n",
                    p.x, p.y, int(0.6*sscale), first ? '+' : ' ', int(shim)
                );   
            }
        }
    }
    
    void print_checks(void) {
        // markers of known dimension to check printed scale
        
        vector<Point2d> corners = {
            {floor(-120 * fiducial_scale), floor(-150 * fiducial_scale)},
            {floor( 120 * fiducial_scale), floor( 150 * fiducial_scale)},
            {floor(-120 * fiducial_scale), floor( 150 * fiducial_scale)},
            {floor( 120 * fiducial_scale), floor(-150 * fiducial_scale)}
            
        };
        
        for (auto p: corners) {
            direct_rectangle(p.x, p.y, 10, 0.25, 0);
            direct_rectangle(p.x, p.y, 0.25, 10, 0);
        }
        
        Point2d ip = corners[0];
        
        Point2d pleft(ip.x + (ip.x < 0 ? - 10 : 10), ip.y);
        
        direct_rectangle(pleft.x, pleft.y, 10, 0.25, 0);
        direct_rectangle(pleft.x-3, pleft.y, 2, 0.25, M_PI/2);
        direct_rectangle(pleft.x+3, pleft.y, 2, 0.25, M_PI/2);
        iPoint p = scale(pleft.x - 2, pleft.y + 1.25);
        fprintf(fout, "  <text x=\"%d\" y=\"%d\" font-family=\"Verdana\" font-size=\"%dmm\" fill=\"black\" > %.0lf mm </text>\n",
            p.x, p.y, int(0.25*sscale), 2*fabs(ip.x)
        );
        
        Point2d ptop(ip.x, ip.y + (ip.y < 0 ? - 10 : 10));
        direct_rectangle(ptop.x, ptop.y, 10, 0.25, M_PI/2);
        direct_rectangle(ptop.x, ptop.y-3, 2, 0.25, 0);
        direct_rectangle(ptop.x, ptop.y+3, 2, 0.25, 0);
        p = scale(ptop.x + 0.5, ptop.y - 2);
        fprintf(fout, "  <text x=\"%d\" y=\"%d\" font-family=\"Verdana\" font-size=\"%dmm\" fill=\"black\" transform=\"rotate(90 %d,%d)\"> %.0lf mm </text>\n",
            p.x, p.y, int(0.25*sscale), 
            p.x, p.y, 
            2*fabs(ip.y)
        );
        
    }
    
    Vec3d cop;      // centre of projection
    Vec3d pl;       // position of focal plane
    Vec3d normal;     // normal vector of target
    double fl;          // focal length of lens
    Vec2d sensor;
    
    double width_scale;
    std::string chart_description;
    int clipid;
    bool alpa_scale;
};

#endif
