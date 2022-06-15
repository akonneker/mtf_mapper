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
#ifndef TIFFSNIFF_H
#define TIFFSNIFF_H

#include <stdio.h>
#include <string>
using std::string;

#include <memory>
using std::shared_ptr;
using std::make_shared;
#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
using std::vector;
using std::pair;
using std::make_pair;

#include "include/display_profile.h"

enum class jpeg_app_t {
    EXIF=1,
    ICC=2,
    NONE=15
};

class Tiffsniff {
  public:
    Tiffsniff(const string& fname, bool is_8bit = false);
    bool profile_found(void) const { return has_profile; }
    Display_profile profile(void);
    
  private:
    enum class profile_t {
        sRGB,
        adobeRGB,
        UNKNOWN,
        CUSTOM
    };
    
    enum class ifd_t {
        ICC,
        TIFF,
        EXIF,
        EXIF_INTEROP
    };
   
    void parse_tiff(off_t offset);
    void read_ifd(off_t offset, off_t base_offset = 0, ifd_t ifd_type = ifd_t::TIFF);
    void read_icc_profile(off_t offset);
    void read_trc_entry(off_t offset, uint32_t size);
    vector<double> read_xyztype_entry(off_t offset, uint32_t size);
    void read_curv_trc(off_t offset, uint32_t size);
    void read_para_trc(off_t offset);
    vector< pair<jpeg_app_t, off_t> > scan_jpeg_app_blocks(void);
    double read_exif_gamma(off_t offset);
    void parse_png(off_t offset);
    
    uint32_t read_uint32(void);
    uint16_t read_uint16(void);
    
    shared_ptr< std::iostream > fin;
    
    bool big_endian = false;
    bool has_profile = false;
    
    int exif_cs = -1;
    bool exif_gamma_found = false;
    bool exif_interop_r03 = false;
    
    profile_t inferred_profile = profile_t::UNKNOWN;
    off_t file_size;
    
    vector<double> gparm {1, 1, 0, 0, 0, 0, 0}; // linear gamma as default
    vector< pair<uint16_t, uint16_t> > gtable;
    vector<double> luminance_weights {0.2225045, 0.7168786, 0.0606169}; // sRGB RGB->Y adapted to D50 by default
};

typedef struct {
    uint16_t tag_id;
    uint16_t data_type;
    uint32_t data_count;
    uint32_t data_offset;
} tiff_field;

typedef struct {
    uint32_t tag_signature;
    uint32_t data_offset;
    uint32_t element_size;
    
    static double read_fixed8_8(shared_ptr< std::iostream > fin) {
        unsigned char b[2];
        fin->read((char*)b, 2);
        return double(b[0]) + double(b[1])/256.0;
    }
    
    static uint16_t read_uint16(shared_ptr< std::iostream > fin) {
        unsigned char b[2];
        fin->read((char*)b, 2);
        return (uint16_t(b[0]) << 8) | uint16_t(b[1]);
    }
    
    static double read_fixed15_16(shared_ptr< std::iostream > fin) {
        unsigned char b[4];
        fin->read((char*)b, 4);
        char sb0 = static_cast<char>(b[0]);
        return double((int16_t(sb0) << 8) | b[1]) + 
            double(((uint16_t(b[2]) << 8) | b[3]))/65536.0;
    }
    
    static uint32_t read_uint32(shared_ptr< std::iostream > fin) {
        unsigned char b[4];
        fin->read((char*)b, 4);
        return (uint32_t(b[0]) << 24) | (uint32_t(b[1]) << 16) |
            (uint32_t(b[2]) << 8) | uint32_t(b[3]);
    }
    
} icc_tag;

#endif
