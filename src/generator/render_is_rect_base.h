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
#ifndef RENDER_IMPORTANCE_SAMPLING_RECT_BASE_H
#define RENDER_IMPORTANCE_SAMPLING_RECT_BASE_H

#include "include/common_types.h"

#include "sinc_sampler.h"
#include "render.h"
#include "polygon_geom.h"
#include "render_importance_sampling.h"

#include <algorithm>

//==============================================================================
class Render_polygon_is_rect_base : public Render_polygon_is {
  public:
    Render_polygon_is_rect_base(Geometry& target, Geometry& photosite,
        Render_type render_type=RECT_PLUS_BOX,
        double in_aperture=8, double in_pitch=4.73, double in_lambda=0.55, double in_aperture2=8, int hs=60) 
        : Render_polygon_is(
            target, photosite, 
            render_type, in_aperture, in_pitch, in_lambda,
            hs // hs=60 for Airy PSF
          ), aperture2(in_aperture2)
          {
          
    }
    
    virtual ~Render_polygon_is_rect_base(void) {
    }
    
    virtual string get_mtf_curve(void) const {
        return string("not defined");
    }
    
    virtual string get_psf_curve(void) const {
        return string("not defined");
    }
    
    virtual double get_mtf50_value(void) const {
        return 0;
    }
      
  protected:
    virtual void initialise(void) {
    
        nsamples = SQR(hs*2 + 1);
        pos_x   = vector<double>(nsamples);
        pos_y   = vector<double>(nsamples);
        weights = vector<double>(nsamples);
        
        normal_sampler sampler;
        Sinc_sampler sinc(lambda, pitch, aperture, aperture2);
        
        vector< pair<double, int> > radius_list;
        
        for (int sidx=0; sidx < nsamples; sidx++) {
            double& ex = pos_x[sidx];
            double& ey = pos_y[sidx];
            
            double sample_prob = sinc.rairy2d(ex, ey, sampler);
            double rad = sqrt(ex*ex + ey*ey);
            
            rad /= (lambda/pitch)*aperture;
            
            // collect the indices of the outermost points
            // this ensures that the early stopping criterion works well
            if (rad > sinc.diam*0.15) {
                radius_list.push_back(make_pair(fabs(rad - sinc.diam*0.8), sidx));
            }
            
            weights[sidx] = sample_prob / sinc.weight(ex, ey);
        } // supersamples
        
        // now swap out the outermost points with the first points
        sort(radius_list.begin(), radius_list.end());
        for (size_t i=0; i < radius_list.size(); i++) {
            swap(pos_x[i], pos_x[radius_list[i].second]);
            swap(pos_y[i], pos_y[radius_list[i].second]);
            swap(weights[i], weights[radius_list[i].second]);
        }
        
        
        printf("using IS renderer with %d samples per pixel (RECT)\n", nsamples);
        #if 0
        FILE* fout = fopen("rect_points.txt", "wt");
        for (int sidx=0; sidx < nsamples; sidx++) {
            fprintf(fout, "%lf %lf %lf\n", pos_x[sidx], pos_y[sidx], weights[sidx]);
        }
        fclose(fout);
        #endif
    }
    
    virtual inline double sample_core(const double&, const double&, const double&, const double&,
        const double&, const double&) const {
        fprintf(stderr, "Render_polygon_is_rect_base::sample_core should never be called\n");
        exit(-1);
        return 0;
    }
    
  public:
    double aperture2;
};

#endif // RENDER_H
