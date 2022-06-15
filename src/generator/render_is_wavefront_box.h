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
#ifndef RENDER_IMPORTANCE_SAMPLING_WAVEFRONT_BOX_H
#define RENDER_IMPORTANCE_SAMPLING_WAVEFRONT_BOX_H

#include "include/common_types.h"
#include "render_is_wavefront_base.h"

//==============================================================================
class Render_polygon_is_wavefront_box : public Render_polygon_is_wavefront_base {
  public:
    Render_polygon_is_wavefront_box(Geometry& target, Geometry& photosite,
        double in_aperture=8, double in_pitch=4.73, double in_lambda=0.55, int hs=60,
        double w020=0.0, double w040=0.0) 
        : Render_polygon_is_wavefront_base(
            target, photosite, 
            WAVEFRONT_PLUS_BOX, in_aperture, in_pitch, in_lambda,
            hs, // hs=60 for Wavefront PSF
            w020, w040
          ) {
          
    }
    
    virtual ~Render_polygon_is_wavefront_box(void) {
    }
    
  protected:
    virtual double mtf_modifier(double f) const {
        return f == 0 ? 1.0 : fabs(sin(f*M_PI)/(f*M_PI));
    }
    
    virtual inline double sample_core(const double& ex, const double& ey, const double& x, const double& y,
        const double& object_value, const double& background_value) const {
        
        double area = target.intersection_area(photosite, ex + x, ey + y) / photosite.own_area;
        return object_value * area + background_value * (1 - area);
    }
};

#endif // RENDER_H
