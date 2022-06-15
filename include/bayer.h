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
#ifndef BAYER_H
#define BAYER_H

#include <string>
using std::string;

class Bayer {
  public:
    typedef enum {
        NONE,
        RED,
        GREEN,
        BLUE
    } bayer_t;
    
    // mask layout, in pixels:
    // A | B
    // C | D
    // => A<<3 | B << 2 | C << 1 | D
    typedef enum {
      DEFAULT = 0x1f, // placeholder used as default argument to samplers
      ALL = 0xf,
      
      RGGB_RED = 1 << 3,
      RGGB_GREEN = (1 << 2) | (1 << 1),
      RGGB_BLUE = 1,
      
      BGGR_RED = 1,
      BGGR_GREEN = (1 << 2) | (1 << 1),
      BGGR_BLUE = 1 << 3,
      
      GRBG_RED = 1 << 2,
      GRBG_GREEN = (1 << 3) | 1,
      GRBG_BLUE = 1 << 1,
      
      GBRG_RED = 1 << 1,
      GBRG_GREEN = (1 << 3) | 1,
      GBRG_BLUE = 1 << 2
      
    } cfa_mask_t;
    
    typedef enum {  // TODO: add some more cfa pattern types
        RGGB,
        BGGR,
        GRBG,
        GBRG
    } cfa_pattern_t;
    
    Bayer(void) {}
    
    static bayer_t from_string(const string& bayer_subset) {
        if (bayer_subset.compare("none") == 0) {
            return NONE;
        }
        if (bayer_subset.compare("red") == 0) {
            return RED;
        }
        if (bayer_subset.compare("blue") == 0) {
            return BLUE;
        }
        if (bayer_subset.compare("green") == 0) {
            return GREEN;
        }
        return NONE; // undefined strings too
    }
    
    static cfa_pattern_t from_cfa_string(const string& cfa_pattern) {
        if (cfa_pattern.compare("rggb") == 0) {
            return cfa_pattern_t::RGGB;
        }
        if (cfa_pattern.compare("bggr") == 0) {
            return cfa_pattern_t::BGGR;
        }
        if (cfa_pattern.compare("grbg") == 0) {
            return cfa_pattern_t::GRBG;
        }
        if (cfa_pattern.compare("gbrg") == 0) {
            return cfa_pattern_t::GBRG;
        }
        return cfa_pattern_t::RGGB; // undefined strings too
    }
    
    static cfa_mask_t to_cfa_mask(bayer_t subset, cfa_pattern_t cfa_pattern = RGGB) {
        switch(cfa_pattern) {
        case RGGB:
            switch (subset) {
            case NONE:  return cfa_mask_t::ALL; break;
            case RED:   return cfa_mask_t::RGGB_RED; break;
            case GREEN: return cfa_mask_t::RGGB_GREEN; break;
            case BLUE:  return cfa_mask_t::RGGB_BLUE; break;
            }
            break;
        case BGGR:
            switch (subset) {
            case NONE:  return cfa_mask_t::ALL; break;
            case RED:   return cfa_mask_t::BGGR_RED; break;
            case GREEN: return cfa_mask_t::BGGR_GREEN; break;
            case BLUE:  return cfa_mask_t::BGGR_BLUE; break;
            }
            break;
        case GRBG:
            switch (subset) {
            case NONE:  return cfa_mask_t::ALL; break;
            case RED:   return cfa_mask_t::GRBG_RED; break;
            case GREEN: return cfa_mask_t::GRBG_GREEN; break;
            case BLUE:  return cfa_mask_t::GRBG_BLUE; break;
            }
            break;
        case GBRG:
            switch (subset) {
            case NONE:  return cfa_mask_t::ALL; break;
            case RED:   return cfa_mask_t::GBRG_RED; break;
            case GREEN: return cfa_mask_t::GBRG_GREEN; break;
            case BLUE:  return cfa_mask_t::GBRG_BLUE; break;
            }
            break;
        }
        return cfa_mask_t::ALL;
    }
    
    static string to_string(bayer_t bayer) {
        string bname="";
        switch (bayer) {
            case RED: bname = "red"; break;
            case BLUE: bname = "blue"; break;
            case GREEN: bname = "green"; break;
            default: bname = ""; break;
        }
        return bname;
    }

};

#endif
