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
#ifndef SVG_PAGE_H
#define SVG_PAGE_H

#include <stdio.h>

#include <cmath>
#include <string>
using std::string;

#ifdef _MSC_VER
	typedef unsigned short int uint16_t;

	#define M_PI 3.14159265358979
	#define lrint(x) ( (x < 0) ? int(floor(x-0.5)) : int(floor(x+0.5)) )
#endif

#include <opencv2/core/core.hpp>
using namespace cv;

#include "include/fiducial_positions.h"

typedef cv::Point_<double> dPoint;
typedef cv::Point_<int> iPoint;

class Svg_page {
  public:
    Svg_page(const string& page_spec, const string& fname, double pscale=100) 
      : fname(fname), page_size(page_spec), sscale(pscale), outside_limits(false), clipid(1)  {
        
        if (page_spec == "A4" || page_spec == "a4") {
            width_mm = 210;
            height_mm = 297;
            fiducial_scale_index = fiducial_mapping_index::A4;
        } 
        if (page_spec == "A3" || page_spec == "a3") {
            width_mm = 297;
            height_mm = 420;
            fiducial_scale_index = fiducial_mapping_index::A3;
        }
        if (page_spec == "A3S" || page_spec == "a3s") {
            width_mm = 420;
            height_mm = 420;
            fiducial_scale_index = fiducial_mapping_index::A3;
        }
        if (page_spec == "A2" || page_spec == "a2") {
            width_mm = 420;
            height_mm = 594;
            fiducial_scale_index = fiducial_mapping_index::A2;
        }
        if (page_spec == "A1" || page_spec == "a1") {
            width_mm = 594;
            height_mm = 841;
            fiducial_scale_index = fiducial_mapping_index::A1;
        }
        if (page_spec == "A0" || page_spec == "a0") {
            width_mm = 841;
            height_mm = 1189;
            fiducial_scale_index = fiducial_mapping_index::A0;
        }
        
        fiducial_scale = fiducial_scale_factor[fiducial_scale_index];
        
        width  = lrint(width_mm * sscale);
        height = lrint(height_mm * sscale);
        
        fout = fopen(fname.c_str(), "wt");
        emit_header();
        
        // set default style
        style="fill:black;stroke:black;stroke-width:1";
    }
    
    ~Svg_page(void) {
        fprintf(fout, "\n\n</svg>\n");
        fclose(fout);
    }
    
    
    virtual void render(void) = 0;
    
    
  protected:  
    void emit_header(void) {
        fprintf(fout, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        fprintf(fout, "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.0\" ");
        fprintf(fout, "width=\"%dmm\" height=\"%dmm\" viewBox=\"0 0 %d %d\">\n", int(width_mm), int(height_mm), int(width), int(height));
        
        // draw a white rectangle to cover the background
        fprintf(fout, "  <rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" style=\"fill:white;stroke:none\"/>\n", 
            int(0), int(0), int(width), int(height)
        );

        // draw a black border rectangle
        fprintf(fout, "  <rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" style=\"fill:none;stroke:black;stroke-width:1\"/>\n", 
            int(4), int(4), int(width-4), int(height-4)
        );
    }
    
    void square(size_t tlx, size_t tly, size_t width) {
        fprintf(fout, "  <rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" style=\"%s\"/>", int(tlx), int(tly), int(width), int(width), style.c_str());
    }

    void rect(size_t tlx, size_t tly, size_t width, size_t height) {
        fprintf(fout, "  <rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" style=\"%s\"/>", int(tlx), int(tly), int(width), int(height), style.c_str());
    }
    
    virtual iPoint project(double x, double y, [[maybe_unused]] bool silent=false) {
        
        x = floor(x*width);
        y = floor(y*height);
        
        return iPoint(int(x), int(y));
    }
    
    virtual iPoint scale(double x, double y) {
    
        dPoint p(x,y);
        
        p.x = floor(p.x*sscale + 0.5*width);
        p.y = floor(p.y*sscale + 0.5*height);
        
        if (p.x < 0 || p.x > width ||
            p.y < 0 || p.y > height) {
            
            outside_limits = true;
        }
        
        return iPoint(int(p.x), int(p.y));
    }
    
    
    void rotated_square(double tlx, double tly, double bwidth, double angle) {
        iPoint p = project(tlx + bwidth*cos(angle), tly + bwidth*sin(angle));
        fprintf(fout, "  <polygon points=\"%d,%d ", p.x, p.y);
        p = project(tlx + bwidth*cos(angle+M_PI/2.0), tly + bwidth*sin(angle+M_PI/2.0));
        fprintf(fout, "%d,%d ", p.x, p.y);
        p = project(tlx + bwidth*cos(angle+M_PI), tly + bwidth*sin(angle+M_PI));
        fprintf(fout, "%d,%d ", p.x, p.y);
        p = project(tlx + bwidth*cos(angle+1.5*M_PI), tly + bwidth*sin(angle+1.5*M_PI));
        fprintf(fout, "%d,%d\" style=\"%s\"/>\n", p.x, p.y, style.c_str());
    }
    
    void text(const std::string& s, double x, double y, int font_size=20) { // x, y are fractional positions
        fprintf(fout, "  <text x=\"%ld\" y=\"%ld\" font-family=\"Verdana\" font-size=\"%ldmm\" fill=\"black\" > %s </text>\n",
            lrint(x*width), lrint(y*height), lrint(font_size), s.c_str()
        );   
    }
    
    void wedge_hole(double cx, double cy, double inner_rad, double outer_rad, double xangle, double fraction, bool double_hole=false) {
        
        if (fraction >= 0.999999) {
            if (!double_hole) {
                hole(cx, cy, outer_rad);
            } else {
                hole(cx, cy, outer_rad);
                hole(cx, cy, inner_rad);
            }
        } else {
            int lf = fraction > 0.5 ? 1 : 0;
            double x = cx - outer_rad * cos((xangle)/180.0*M_PI);
            double y = cy - outer_rad * sin((xangle)/180.0*M_PI);
            iPoint p = scale(x, y);
            fprintf(fout, "M %d,%d ", p.x, p.y);
            x = cx - outer_rad * cos((xangle + 360*fraction)/180.0*M_PI);
            y = cy - outer_rad * sin((xangle + 360*fraction)/180.0*M_PI);
            p = scale(x, y);
            fprintf(fout, "A %ld %ld 0 %d,1 %d,%d ", lrint(sscale*(outer_rad)), lrint(sscale*(outer_rad)), lf, p.x, p.y);
            x = cx - inner_rad * cos((xangle + 360*fraction)/180.0*M_PI);
            y = cy - inner_rad * sin((xangle + 360*fraction)/180.0*M_PI);
            p = scale(x, y);
            fprintf(fout, "L %d,%d ", p.x, p.y);
            x = cx - inner_rad * cos((xangle)/180.0*M_PI);
            y = cy - inner_rad * sin((xangle)/180.0*M_PI);
            p = scale(x, y);
            fprintf(fout, "A %ld,%ld 0 %d,0 %d,%d z\n", lrint(sscale*inner_rad), lrint(sscale*inner_rad), lf, p.x, p.y);
        }
    }
    
    void ring_open(double cx, double cy, double outer_rad) {
        int ior = lrint(sscale*outer_rad);
        iPoint p;
        p = scale(cx + outer_rad, cy);
        fprintf(fout, " <path \n");
        fprintf(fout, "\td=\"M %d,%d ", p.x, p.y);
        p = scale(cx - outer_rad, cy);
        fprintf(fout, "A %d %d 0 0,1 %d,%d ", ior, ior, p.x, p.y);
        p = scale(cx + outer_rad, cy);
        fprintf(fout, "A %d %d 0 0,1 %d,%d z", ior, ior, p.x, p.y);
    }
    
    void hole(double cx, double cy, double outer_rad) {
        int ior = lrint(sscale*outer_rad);
        iPoint p;
        p = scale(cx + outer_rad, cy);
        fprintf(fout, "M %d,%d ", p.x, p.y);
        p = scale(cx - outer_rad, cy);
        fprintf(fout, "A %d %d 0 0,1 %d,%d ", ior, ior, p.x, p.y);
        p = scale(cx + outer_rad, cy);
        fprintf(fout, "A %d %d 0 0,1 %d,%d z", ior, ior, p.x, p.y);
    }

    void ring(double cx, double cy, double inner_rad, double outer_rad) {
        int ior = lrint(sscale*outer_rad);
        int iir = lrint(sscale*inner_rad);
        iPoint p;
        p = scale(cx + outer_rad, cy);
        fprintf(fout, " <path \n");
        fprintf(fout, "\td=\"M %d,%d ", p.x, p.y);
        p = scale(cx - outer_rad, cy);
        fprintf(fout, "A %d %d 0 0,1 %d,%d ", ior, ior, p.x, p.y);
        p = scale(cx + outer_rad, cy);
        fprintf(fout, "A %d %d 0 0,1 %d,%d z", ior, ior, p.x, p.y);
        p = scale(cx + inner_rad, cy);
        fprintf(fout, "M %d,%d \n", p.x, p.y);
        p = scale(cx - inner_rad, cy);
        fprintf(fout, "A %d %d 0 0,1 %d,%d ", iir, iir, p.x, p.y);
        p = scale(cx + inner_rad, cy);
        fprintf(fout, "A %d %d 0 0,1 %d,%d z", iir, iir, p.x, p.y);
        fprintf(fout, "\" fill-rule=\"evenodd\" />\n");
    }
    
    void sector_circle(double cx, double cy, double rad, int code) {
        
        if (code == 0) {
            ring(cx, cy, 0.12*rad, rad);
            return;
        }
        
        const int outer_step = 5;
        const int inner_step = 3;
        
        const int outer_vals = outer_step + 1;
        const int inner_vals = inner_step + 1;
        
        int inner_code = inner_step - code % inner_vals;
        int outer_code = outer_step - (code / inner_vals) % outer_vals;
        
        
        ring_open(cx, cy, rad);
        if (outer_code >= 0) {
            wedge_hole(cx, cy, 0.4*rad, 0.6*rad, 90, 1 - outer_code*1.0/double(outer_step), true);
        }
        if (inner_code >= 0) {
            wedge_hole(cx, cy, 0.12*rad, 0.4*rad, 270, 1 - inner_code*1.0/double(inner_step));
        }
        if (inner_code > 0) {
            hole(cx, cy, 0.12*rad);
        }
        fprintf(fout, "\" fill-rule=\"evenodd\" />\n");
    }
    
    
    string fname;
    string style;
    FILE*  fout;
    
    size_t width_mm;   // in mm
    size_t height_mm;  // in mm
    size_t width;   // in svg units
    size_t height;  // in svg units
    std::string page_size;
    
    double sscale;
    bool outside_limits;
    int clipid;
    
    fiducial_mapping_index fiducial_scale_index = fiducial_mapping_index::A3;
    double fiducial_scale = 1.0;
};

#endif
