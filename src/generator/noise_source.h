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

#ifndef NOISE_SOURCE_H
#define NOISE_SOURCE_H

class Noise_source {
  public:
    Noise_source(int size) : size(size) {
    }

    virtual ~Noise_source() {

    }

    virtual double sample(double s, int idx) {
        return clip(s + 0.0*idx);
    }
    
    virtual double clip(double x) {
        clip_count += x < 0 || x > 1.0;
        return std::max(0.0, std::min(1.0, x));
    }

    double randdu(void) {
        return double(rand())/double(RAND_MAX);
    }
    
    double rand_norm(double m, double s) {
        double x = normal_sampler::Moro_norm_inv(randdu());
        return( m + x * s );
    }

    size_t size;
    size_t clip_count = 0;
};

class Additive_gaussian_noise : public Noise_source {
  public:
    Additive_gaussian_noise(size_t size, double sigma) 
      : Noise_source(size), noise(size) {
        for (size_t i=0; i < size; i++) {
            noise[i] = rand_norm(0, sigma);
        }
    }

    virtual double sample(double s, int idx) {
        return clip(s + noise[idx]);
    }

    vector<double> noise;
};

// Sensor noise model described in chapter 5 of
// James R. Janesick, "Photon Transfer", SPIE Press, 2007
class Sensor_model_noise : public Noise_source {
  public:
    Sensor_model_noise(size_t size, 
        double read_noise=3, // readout noise, in electrons
        double Pn=0.02,      // Fixed Pattern Noise quality factor (percentage)
        double adc=1.5,
        size_t adc_bits=14)      // ADC e/DN factor?

      : Noise_source(size), noise(3, vector<double>(size)),
        read_noise(read_noise), Pn(Pn), adc(adc),
        adc_bits(adc_bits), DN_scale(1 << adc_bits) {

        for (size_t j=0; j < 3; j++) {
            for (size_t i=0; i < size; i++) {
                noise[j][i] = rand_norm(0, 1);
            }
        }
    }

    virtual double sample(double s, int idx) {
        double signal = s * DN_scale;

        double total_noise  = 
            noise[0][idx]*read_noise/adc +
            noise[1][idx]*sqrt(signal/adc) +
            noise[2][idx]*Pn*signal;

        return clip((signal + total_noise) / 65535.0);
    }
    
    virtual double clip(double x) {
        clip_count += x < 0 || x > (DN_scale-1.0) / 65535.0;
        return std::max(0.0, std::min((DN_scale-1.0) / 65535.0, x));
    }

    vector< vector<double> > noise;
    double read_noise;
    double Pn;
    double adc;
    size_t adc_bits;
    double DN_scale;
};

        

#endif

