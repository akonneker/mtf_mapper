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
#ifndef RENDER_IMPORTANCE_SAMPLING_AIRY_H
#define RENDER_IMPORTANCE_SAMPLING_AIRY_H

#include "include/common_types.h"

#include "airy_sampler.h"
#include "render.h"
#include "polygon_geom.h"
#include "render_is_airy_base.h"

//==============================================================================
class Render_polygon_is_airy : public Render_polygon_is_airy_base {
  public:
    Render_polygon_is_airy(Geometry& target, Geometry& photosite,
        double in_aperture=8, double in_pitch=4.73, double in_lambda=0.55, int hs=60) 
        : Render_polygon_is_airy_base(
            target, photosite, 
            AIRY, in_aperture, in_pitch, in_lambda,
            hs // hs=60 for Airy PSF
          )
          {
          
          initialise();
    }
    
    virtual ~Render_polygon_is_airy(void) {
    }
    
    virtual string get_mtf_curve(void) const {
        char buffer[1024];
        double scale = (lambda/pitch) * aperture;
        sprintf(buffer, "2.0/pi * (acos(x*%lg) - (x*%lg)*sqrt(1-(x*%lg)*(x*%lg)))", 
            scale, scale, scale, scale
        );
        return string(buffer);
    }
    
    virtual string get_psf_curve(void) const {
        char buffer[1024];
        double scale = 1 / ((lambda/pitch) * aperture);
        sprintf(buffer, "(2*besj1((x*%lg))/(x*%lg))**2", scale, scale);
        return string(buffer);
    }
    
    virtual double get_mtf50_value(void) const {
        double scale = (lambda/pitch) * aperture;
        return 0.4039728 / scale;
    }
      
  protected:
    virtual inline double sample_core(const double& ex, const double& ey, const double& x, const double& y,
        const double& object_value, const double& background_value) const {
        
        return t_geom.is_inside(ex + x, ey + y) ? object_value : background_value;
    }
};

#endif // RENDER_H
