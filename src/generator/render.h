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
#ifndef RENDER_H
#define RENDER_H

#include "include/common_types.h"

const double transparent = -1.0;

#include "normal_sampler.h"
#include <string>
using std::string;

#include "geom.h"
#include "polygon_geom.h"

#include <random>

//==============================================================================
class Render_target {
  public:
      Render_target(void) {}
      virtual double evaluate(double x, double y, double oject_value, double background_value) const = 0;
};


//==============================================================================
class Render_polygon : public Render_target {
  public:
    typedef enum {
        GAUSSIAN,
        GAUSSIAN_SAMPLED,
        AIRY,
        AIRY_PLUS_BOX,
        AIRY_PLUS_4DOT_OLPF,
        WAVEFRONT,
        WAVEFRONT_PLUS_BOX,
        RECT,
        RECT_PLUS_BOX
    } Render_type;
  
    Render_polygon(Geometry& target, 
        double in_sigma=6.0, double minor_sigma=6.0, double theta=0, bool init=true) 
        : sigma(in_sigma), 
          t_geom(target),
          psf_ratio(-1),
          img_width(0),
          img_height(0) {
        
        
        if (init) {
            initialize_samples(22, fabs(minor_sigma - 6.0) < 1e-6 ? sigma : minor_sigma, theta); // seems like enough samples for up to sigma=6, at least
        }
    }
    
    virtual ~Render_polygon(void) {
    }
    
    void initialize_samples(int hs=22, double minor_sigma=0, double theta=0) {
        this->hs = hs;
        int nsamples = SQR(hs*2 + 1);
        pos_x   = vector<double>(nsamples);
        pos_y   = vector<double>(nsamples);
    
        normal_sampler sampler;
        for (int sidx=0; sidx < nsamples; sidx++) {
            sampler.rnorm2d(
                pos_x[sidx], pos_y[sidx], 
                sigma, 
                minor_sigma == 0 ? sigma : minor_sigma, 
                theta
            );
        } 
    }
    
    double evaluate(double x, double y, double object_value, double background_value) const {
   
        double accumulator = 0;
        for (size_t sidx=0; sidx < pos_x.size(); sidx++) {
                if ( t_geom.is_inside(pos_x[sidx] + x, pos_y[sidx] + y) ) {
                    accumulator += object_value;
                } else {
                    accumulator += background_value;
                }
        } // supersamples
         
        double value = accumulator / ((2*hs+1)*(2*hs+1));
        return value;
    }
    
    virtual string get_mtf_curve(void) const {
        char buffer[1024];
        double a = 1.0/(sigma*sigma);
        sprintf(buffer, "exp(%.8lg*x*x)", -2*M_PI*M_PI/a);
        return string(buffer);
    }
    
    virtual string get_psf_curve(void) const {
        char buffer[1024];
        sprintf(buffer, "exp(-x*x/%lg)", 2*sigma*sigma);
        return string(buffer);
    }
    
    virtual double get_mtf50_value(void) const {
        return sqrt( log(0.5)/(-2*M_PI*M_PI*sigma*sigma) );
    }
    
    void set_img_dimensions(int rows, int cols, bool random_angle=true, double angle=0) {
        img_width = cols;
        img_height = rows;
        
        
        std::mt19937 gen(10);
        std::uniform_real_distribution<float> dis(0, 2*M_PI);
        
        uint64_t isize = uint64_t(img_height)*uint64_t(img_width);
        rseed = vector<float>(isize);
        for (uint64_t i=0; i < isize; i++) {
            rseed[i] = random_angle ? dis(gen) : angle;
        }
    }
    
    void set_psf_ratio(double ratio) {
        psf_ratio = ratio;
    }
    
  public:
      
    double sigma;
    int    hs;
    
    vector<double> pos_x;
    vector<double> pos_y;

    Geometry& t_geom; // target geometry
    double psf_ratio;
    int img_width;
    int img_height;
    vector<float> rseed;
};

#endif // RENDER_H
