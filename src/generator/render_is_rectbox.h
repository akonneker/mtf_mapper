/*
Copyright 2021 Frans van den Bergh. All rights reserved.

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
#ifndef RENDER_IMPORTANCE_SAMPLING_RECTBOX_H
#define RENDER_IMPORTANCE_SAMPLING_RECTBOX_H

#include "include/common_types.h"

#include "airy_sampler.h"
#include "render.h"
#include "polygon_geom.h"
#include "render_is_rect_base.h"

//==============================================================================
class Render_polygon_is_rectbox : public Render_polygon_is_rect_base {
  public:
    Render_polygon_is_rectbox(Geometry& target, Geometry& photosite, 
        double in_aperture=8, double in_pitch=4.73, double in_lambda=0.55, double in_aperture2=8, int hs=40) 
        : Render_polygon_is_rect_base(
            target, photosite, 
            RECT_PLUS_BOX, in_aperture, in_pitch, in_lambda, in_aperture2,
            hs // #half-samples
          )
          {
          
          initialise();
    }
    
    virtual ~Render_polygon_is_rectbox(void) {
    }
    
    virtual string get_mtf_curve(void) const {
        char buffer[1024];
        double scale = (lambda/pitch) * aperture;
        sprintf(buffer, "(x < %lf) ? abs(x > 0 ? sin(x*pi)/(x*pi) : 1)*(1 - x*%lf) : 0", 
            1.0/scale, scale
        );
        return string(buffer);
    }
    
    virtual string get_psf_curve(void) const {
        return string("not implemented (no simple analytical form)");
    }
    
    virtual double get_mtf50_value(void) const {
        return  bisect_airy(&rect_box_mtf, 0); // TODO: update
    }
      
  protected:
    virtual inline double sample_core(const double& ex, const double& ey, const double& x, const double& y,
        const double& object_value, const double& background_value) const {
    
        double area = target.intersection_area(photosite, ex + x, ey + y) / photosite.own_area;
        return object_value * area + background_value * (1 - area);
    }
    
    double static rect_box_mtf(double x, double s, double)  {
        double cut = 1.0 / s;
        return (x >= cut ? 0 : fabs(x > 0 ? sin(x*M_PI)/(x*M_PI) : 1) * (1 - x*s)) - 0.5;
    }
};

#endif // RENDER_H
