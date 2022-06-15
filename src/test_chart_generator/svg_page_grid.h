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
#ifndef SVG_PAGE_GRID_H
#define SVG_PAGE_GRID_H

#include "svg_page.h"
#include "include/mtf50_edge_quality_rating.h"

typedef enum {INVALIDGRID, FULLGRID, HALFGRID, THIRDGRID} grid_scale;

class Svg_page_grid : public Svg_page {
  public:
    Svg_page_grid(const string& page_spec, const string& fname, grid_scale scale=FULLGRID) 
    : Svg_page(page_spec, fname), centre(0.5, 0.5*sqrt(2.0)), scale(scale) {
    
    }
    
    
    
    void render(void) {
        switch(scale) {
        case FULLGRID: 
            grid(0.04, 17, 11); 
            break;
        case HALFGRID:
            grid(0.08, 9, 6);
            break;
        case THIRDGRID:
            grid(0.12, 6, 4);
            break;
        default:
            printf("Invalid grid type. Refusing to render\n");
            break;
        };
    }
    
  protected:
    
    virtual iPoint project(double x, double y, bool) {
        
        x = floor(x*width);
        y = floor(y*height*sqrt(0.5));
        
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
    
    void place_square(double xpos, double ypos, double width) {
        dPoint delta = dPoint(xpos, ypos) - centre;
        
        double eang = atan2(fabs(delta.y), fabs(delta.x));
        double degrees = eang / M_PI * 180;
        double orient = degrees;
        if (orient > 45) {
            orient = 90 - orient;
        }
        
        const double tolerance = 3.0;
        const int ncrit = 4;
        double crit[ncrit] = {0, 14.036, 26.565, 45};
        double correction = 0;
        
        for (int k=0; k < ncrit; k++) {
            double del = orient - crit[k];
            if (fabs(del) < tolerance) {
                if (del < 0) {
                    correction = orient - (crit[k] + tolerance);
                } else {
                    correction = orient - (crit[k] - tolerance);
                }
            }
        }
        
        degrees -= correction;
        eang = degrees / 180.0 * M_PI + M_PI/4;
        
        rotated_square(xpos, ypos, width*0.5, eang);
    }

    void grid(double swidth, size_t nrows, size_t ncols) {
        int hr = nrows/2;
        int hc = ncols/2;
        for (int col=-hc; col <= hc; col++) {
            
            for (int row=-hr; row <= hr; row++) {
                place_square((col+hc+0.5)/double(2*hc+1), sqrt(2.0)*(row+hr+0.5)/double(2*hr+1), swidth);
            }
            
            // in-between columns
            if (col < hc) {
                for (int row=-hr; row < hr; row++) {
                    place_square((col+hc+1)/double(2*hc+1), sqrt(2.0)*(row+hr+1)/double(2*hr+1), swidth);
                }
            }
        }
    }
    
    dPoint centre;
    grid_scale scale;
};

#endif
