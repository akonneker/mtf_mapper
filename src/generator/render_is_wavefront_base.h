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
#ifndef RENDER_IMPORTANCE_SAMPLING_WAVEFRONT_BASE_H
#define RENDER_IMPORTANCE_SAMPLING_WAVEFRONT_BASE_H

#include "include/common_types.h"

#include "wavefront_sampler.h"
#include "render.h"
#include "polygon_geom.h"

#include "opencv2/opencv.hpp"

//==============================================================================
class Render_polygon_is_wavefront_base : public Render_polygon_is {
  public:
    Render_polygon_is_wavefront_base(Geometry& target, Geometry& photosite,
        Render_type render_type=WAVEFRONT,
        double in_aperture=8, double in_pitch=4.73, double in_lambda=0.55, int hs=60,
        double w020=0.0, double w040=0.0) 
        : Render_polygon_is(
            target, photosite, 
            render_type, in_aperture, in_pitch, in_lambda,
            hs // hs=60 for Wavefront PSF
          ), w020(w020), w040(w040)
          {
          
          initialise();
    }
    
    virtual ~Render_polygon_is_wavefront_base(void) {
    }
    
    virtual string get_mtf_curve(void) const {
        char buffer[1024];
        sprintf(buffer, "no analytical form\n");
        return string(buffer);
    }
    
    virtual string get_psf_curve(void) const {
        char buffer[1024];
        sprintf(buffer, "no analytical form\n");
        return string(buffer);
    }
    
    // #define DEBUG_PSF
    
    virtual double get_mtf50_value(void) const {
        
        #ifdef DEBUG_PSF
        printf("get_mtf50_value() called on wavefront_based\n");
        
        FILE* fout = fopen("rpsf.txt", "wt");
        for (double r=0; r < Wavefront_sampler::diam; r += 0.01) {
            fprintf(fout, "%le %le\n", r, wave.get_psf(r)); 
        }
        fclose(fout);
        #endif
        
        const double psf_width = 200*aperture;
        constexpr int fft_width = 1024;
        
        cv::Mat psf(fft_width, fft_width, CV_32FC1);
        const double radnorm = psf.rows/2.0 * (1.0/psf_width);
        
        // find PSF peak in the range that we will be using
        // because peak is not always at radius = 0
        double maxpsf = 0;
        for (int c=0; c < psf.cols/2+1; c++) {
            double rad = fabs(c - psf.cols/2) / radnorm;
            maxpsf = std::max(maxpsf, wave.get_psf(rad/aperture));
        }
        
        for (int r=0; r < psf.rows; r++) {
            for (int c=0; c < psf.cols; c++) {
                double rad = sqrt( (r - psf.rows/2)*(r - psf.rows/2) + (c - psf.cols/2)*(c - psf.cols/2)) / radnorm;
                float val = wave.get_psf(rad/aperture) / maxpsf;
                psf.at<float>(r, c) = val;
            }
        }
        
        
        cv::Mat mtf;
        cv::dft(psf, mtf);
        
        #ifdef DEBUG_PSF
        for (int r=0; r < psf.rows; r++) {
            for (int c=0; c < psf.cols; c++) {
                psf.at<float>(r, c) = sqrt(sqrt(psf.at<float>(r, c)));
            }
        }
        
        cv::Mat psf8(psf.rows, psf.cols, CV_8UC1);
        psf.convertTo(psf8, CV_8U, 255.0, 0);
        imwrite("psf8.png", psf8);
        #endif
        
        double mnorm = fabs(mtf.at<float>(0, 0));
        double prev_mag = 1.0;
        double prev_f = 0.0;
        double mtf50 = 0;
        for (int c=1; c < mtf.cols/2; c += 2) {
            double f = (0.5*(c+1))*pitch / (2*psf_width*lambda);
            double mod = sqrt(sqr(mtf.at<float>(0, c)) + sqr(mtf.at<float>(0, c+1)));
            double mag = mtf_modifier(f) * mod / mnorm;
            
            if (mtf50 == 0 && mag < 0.5 && prev_mag >= 0.5) {
                double dy = f - prev_f;
                double dx = mag - prev_mag;
                double slope = dy/dx;
                mtf50 = prev_f + slope*(0.5 - prev_mag);
            }
            prev_mag = mag;
            prev_f = f;
        }
        
        return mtf50;
    }
      
  protected:
    virtual double mtf_modifier([[maybe_unused]] double f) const {
        return 1.0;
    }
  
    template<class T> inline T sqr(T x) const {
        return x*x;
    }
    
    virtual void initialise(void) {
    
        nsamples = SQR(hs*2 + 1);
        pos_x   = vector<double>(nsamples);
        pos_y   = vector<double>(nsamples);
        weights = vector<double>(nsamples);
        
        normal_sampler sampler;
        wave = Wavefront_sampler(lambda, pitch, aperture, w020, w040);
        
        vector< pair<double, int> > radius_list;
        
        for (int sidx=0; sidx < nsamples; sidx++) {
            double& ex = pos_x[sidx];
            double& ey = pos_y[sidx];
            
            [[maybe_unused]] double sample_prob = wave.rairy2d(ex, ey, sampler);
            double rad = sqrt(ex*ex + ey*ey);
            
            rad /= (lambda/pitch)*aperture;
            
            // collect the indices of the outermost points
            // this ensures that the early stopping criterion works well
            if (rad > Wavefront_sampler::diam*0.15) {
                radius_list.push_back(make_pair(fabs(rad - Wavefront_sampler::diam*0.8), sidx));
            }
            
            weights[sidx] = 1.0; 
        } // supersamples
        
        // now swap out the outermost points with the first points
        sort(radius_list.begin(), radius_list.end());
        for (size_t i=0; i < radius_list.size(); i++) {
            swap(pos_x[i], pos_x[radius_list[i].second]);
            swap(pos_y[i], pos_y[radius_list[i].second]);
            swap(weights[i], weights[radius_list[i].second]);
        }
        
        printf("using IS renderer with %d samples per pixel\n", nsamples);
    }
  
    virtual inline double sample_core(const double& ex, const double& ey, const double& x, const double& y,
        const double& object_value, const double& background_value) const {
        
        return t_geom.is_inside(ex + x, ey + y) ? object_value : background_value;
    }
    
    Wavefront_sampler wave;
    double w020;
    double w040;
};

#endif // RENDER_H
