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
#ifndef RENDER_IMPORTANCE_SAMPLING_AIRYOLPF_H
#define RENDER_IMPORTANCE_SAMPLING_AIRYOLPF_H

#include "include/common_types.h"

#include "airy_sampler.h"
#include "render.h"
#include "polygon_geom.h"
#include "render_is_airy_base.h"

//==============================================================================
class Render_polygon_is_airyolpf : public Render_polygon_is_airy_base {
  public:
    Render_polygon_is_airyolpf(Geometry& target, Geometry& photosite, 
        double in_aperture=8, double in_pitch=4.73, double in_lambda=0.55, double olpf_split=0.375, 
        int hs=30) 
        : Render_polygon_is_airy_base(
            target, photosite, 
            AIRY_PLUS_4DOT_OLPF, in_aperture, in_pitch, in_lambda,
            hs // half-samples
          ),
		  olpf_split(olpf_split)
		  {
          
          initialise();
          
          // construct a new bounding box around all four the supports
          bounds = photosite.bounds;
          bounds.max_x += olpf_split;
          bounds.min_x -= olpf_split;
          bounds.max_y += olpf_split;
          bounds.min_y -= olpf_split;
    }

    virtual ~Render_polygon_is_airyolpf(void) {
    }
    
    virtual string get_mtf_curve(void) const {
        char buffer[1024];
        double scale = (lambda/pitch) * aperture;
        sprintf(buffer, "2.0/pi*abs( %g*sin(%g*x*pi)/(%g*x*pi) + %g*sin(%g*x*pi)/(%g*x*pi) )*(acos(x*%lg) - (x*%lg)*sqrt(1-(x*%lg)*(x*%lg)))", 
			0.5 + olpf_split, 2.0*(0.5+olpf_split), 2.0*(0.5+olpf_split),
			0.5 - olpf_split, 2.0*(0.5-olpf_split), 2.0*(0.5-olpf_split),
            scale, scale, scale, scale
        );
        return string(buffer);
    }
    
    virtual string get_psf_curve(void) const {
        return(string("not implemented (no simple analytical form)"));
    }
    
    virtual double get_mtf50_value(void) const {
        return bisect_airy(&airy_olpf_mtf, olpf_split);
    }
      
  protected:

	double olpf_split;
	Aa_bb  bounds;
	
    virtual inline double sample_core(const double& ex, const double& ey, const double& x, const double& y,
        const double& object_value, const double& background_value) const {
        
        // this is a very minor performance boost, since the indivdual
        // support bounds checks are quite efficient
        if (!bounds.bounds_overlap(target.bounds, ex + x, ey + y)) {
            return background_value;
        }
    
        double sample = 0;
        
		double olpf_x[4] = {-olpf_split, -olpf_split,  olpf_split, olpf_split};
		double olpf_y[4] = {-olpf_split,  olpf_split, -olpf_split, olpf_split};
		        
        for (int k=0; k < 4; k++) {
            double area = target.intersection_area(photosite, ex + x + olpf_x[k], ey + y + olpf_y[k]) / photosite.own_area;
            sample += 0.25*(object_value * area + background_value * (1 - area));
        }
                
        return sample;
    }
    
    double static airy_olpf_mtf(double x, double s, double p_olpf_split)  {
        double wide = 0.5 + p_olpf_split;
        double narrow = 0.5 - p_olpf_split;
        return 2.0/M_PI * 
          fabs( wide*sin(2*wide*x*M_PI)/(2*wide*x*M_PI) + 
                narrow*sin(2*narrow*x*M_PI)/(2*narrow*x*M_PI) ) * 
                (acos(x*s) - (x*s)*sqrt(1-(x*s)*(x*s))
              ) - 0.5;
    }
};

#endif // RENDER_H
