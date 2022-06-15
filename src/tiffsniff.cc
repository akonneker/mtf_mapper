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

#include "include/tiffsniff.h"
#include "include/logger.h"
#include "include/common_types.h"
#include <string.h>

#include "config.h"
#if mtfmapper_ZLIB_FOUND == 1
    #include <zlib.h>
#endif

Tiffsniff::Tiffsniff(const string& fname, bool is_8bit) {

    fin = shared_ptr<std::iostream>(new std::fstream(fname, std::ios_base::binary | std::ios_base::in));
    if (fin->good()) {
        if (fin->seekg(0, fin->end).good()) {
            file_size = off_t(fin->tellg());
        }
        fin->seekg(0);
        
        try {
        
            unsigned char magic[4];
            fin->read((char*)magic, 4);
            
            if (!fin->good()) {
                throw -1;
            }
        
            if ((magic[0] == 0x4d && magic[1] == 0x4d) ||
                (magic[0] == 0x49 && magic[1] == 0x49)) { // possible TIFF file
                
                parse_tiff(0);
            }
            
            if (magic[0] == 0xff && magic[1] == 0xd8) { // possible JPEG file
                // we could check the end of the file for 0xff 0xf9, but Sony pads the back
                // of the file with zeros, so we'll just wing it
                has_profile = true;
                // now build a list of APP blocks
                auto blocks = scan_jpeg_app_blocks();
                
                if (blocks.size() > 0) {
                    bool icc_found = false;
                    // if we have an ICC block, use that
                    for (size_t i=0; i < blocks.size() && !icc_found; i++) {
                        if (blocks[i].first == jpeg_app_t::ICC) {
                            icc_found = true;
                            read_icc_profile(blocks[i].second);
                        }
                    }
                    
                    bool exif_found = false;
                    if (!icc_found) { // no ICC block, settle for EXIF
                        for (size_t i=0; i < blocks.size() && !exif_found; i++) {
                            if (blocks[i].first == jpeg_app_t::EXIF) {
                                exif_found = true;
                                parse_tiff(blocks[i].second);
                            }
                        }
                    }
                }
            }
            
            if (memcmp(magic, "\x89\x50\x4e\x47", 4) == 0) {
                parse_png(0);
            }
            
        } catch (int) {
            logger.error("%s\n", "Error while parsing input image, unable to extract colourspace information.");
        }
    }
    
    if (inferred_profile == profile_t::UNKNOWN) {
        logger.debug("No ICC profile found, attempting to infer profile from Exif information, Exif CS: %d\n", exif_cs);
        if (exif_cs < 0) {
            // for 16-bit images, force linear colour space if no Exif Colour Space tag was found
            // for 8-bit images, assume sRGB
            inferred_profile = is_8bit ? profile_t::sRGB : profile_t::CUSTOM; 
        } else {
            if (exif_cs == 1) {
                // sRGB confirmed
                inferred_profile = profile_t::sRGB;
            } else {
                // otherwise assume this means aRGB
                inferred_profile = profile_t::adobeRGB;
                if (!exif_gamma_found) {
                    gparm[0] = 563/256.0;
                }
            }
        }
    }
    
    if (inferred_profile == profile_t::sRGB) {
        gparm[0] = 2.4;
        gparm[1] = 1.0/1.055;
        gparm[2] = 0.055/1.055;
        gparm[3] = 1.0/12.92;
        gparm[4] = 0.04045;
        
        luminance_weights[0] = 0.2225045;
        luminance_weights[1] = 0.7168786;
        luminance_weights[2] = 0.0606169;
    }
    if (inferred_profile == profile_t::adobeRGB) {
        luminance_weights[0] = 0.3111242;
        luminance_weights[1] = 0.6256560;
        luminance_weights[2] = 0.0632197;
    }
}

uint32_t Tiffsniff::read_uint32(void) {
    if (big_endian) {
        return icc_tag::read_uint32(fin);
    } 
    unsigned char b[4];
    fin->read((char*)b, 4);
    return uint32_t(b[0]) | (uint32_t(b[1]) << 8) |
        (uint32_t(b[2]) << 16) | (uint32_t(b[3]) << 24);
}

uint16_t Tiffsniff::read_uint16(void) {
    if (big_endian) {
        return icc_tag::read_uint16(fin);
    } 
    unsigned char b[2];
    fin->read((char*)b, 2);
    return uint32_t(b[0]) | (uint32_t(b[1]) << 8);
}

Display_profile Tiffsniff::profile(void) {
    if (gtable.size() == 0) {
        return Display_profile(gparm, luminance_weights);
    }
    return Display_profile(gtable, luminance_weights);
}

void Tiffsniff::parse_tiff(off_t offset) {
    if (fin->seekg(offset).good()) {
        const char be_id[4] = {0x4D, 0x4D, 0x00, 0x2A};
        const char le_id[4] = {0x49, 0x49, 0x2A, 0x00};
        
        unsigned char magic[4];
        
        if (fin->read((char*)magic, 4).good()) {
            // TIFF
            bool is_valid_tiff = false;
            if (memcmp(le_id, magic, 4) == 0) {
                if (offset == 0) {
                    logger.debug("%s\n", "Little endian TIFF detected.");
                } else {
                    logger.debug("%s\n", "EXIF / Little endian TIFF detected.");
                }
                is_valid_tiff = true;
            }
            if (memcmp(be_id, magic, 4) == 0) {
                if (offset == 0) {
                    logger.debug("%s\n", "Big endian TIFF detected.");
                } else {
                    logger.debug("%s\n", "EXIF / Big endian TIFF detected.");
                }
                is_valid_tiff = true;
                big_endian = true;
            }
            
            if (is_valid_tiff) {
                uint32_t ifd_offset = read_uint32();
                read_ifd(ifd_offset + offset, offset, ifd_t::TIFF);
                has_profile = true;
            } 
        } else {
            throw -1;
        }
    } else {
        throw -1;
    }
}

void Tiffsniff::read_ifd(off_t offset, off_t base_offset, ifd_t ifd_type) {
    if (fin->seekg(offset).good()) {
        uint16_t cur_ifd_entries = read_uint16();
        
        // check for weird values; realistic tiff files will not have thousands of entries
        if (cur_ifd_entries == 0 || cur_ifd_entries > 32765 || !fin->good()) {
            throw -1;
        }
        
        tiff_field field;
        for (uint16_t i=0; i < cur_ifd_entries; i++) {
            field.tag_id = read_uint16();
            field.data_type = read_uint16();
            field.data_count = read_uint32();
            field.data_offset = read_uint32();;
            
            if (i < (cur_ifd_entries-1) && !fin->good()) {
                throw -1;
            }
            
            
            auto fpos = fin->tellg();
            
            if (ifd_type == ifd_t::TIFF && field.tag_id == 0x8773) { // ICC profile 
                read_icc_profile(field.data_offset);
            }
            
            if (ifd_type == ifd_t::TIFF && field.tag_id == 0x8769) { // EXIF field
                read_ifd(base_offset + field.data_offset, base_offset, ifd_t::EXIF);
            }
            
            if (ifd_type == ifd_t::EXIF && field.tag_id == 0xa001) { // EXIF colourspace tag
                if (fin->seekg(-4, fin->cur).good()) {
                    exif_cs = read_uint16();
                }
            }
            
            if (ifd_type == ifd_t::EXIF && field.tag_id == 0xa500) { // EXIF gamma tag
                double gamma = read_exif_gamma(base_offset + field.data_offset);
                gparm[0] = gamma;
                exif_gamma_found = true;
            }
            
            if (ifd_type == ifd_t::EXIF && field.tag_id == 0xa005) { // EXIF interoperability IFD pointer
                read_ifd(base_offset + field.data_offset, base_offset, ifd_t::EXIF_INTEROP);
            }
            
            if (ifd_type == ifd_t::EXIF_INTEROP && field.tag_id == 0x0001) { // EXIF interoperability tag
                char ids[4];
                if (fin->seekg(-4, fin->cur).good() && fin->read((char*)ids, 4).good()) {
                    exif_interop_r03 = strncmp(ids, "R03", 3) == 0;
                }
            }
            
            fin->seekg(fpos);
        }
        
        // read next IDF offset
        uint32_t next_offset = read_uint32();
        if (next_offset) {
            if (next_offset > file_size || !fin->good()) {
                throw -1;
            }
            read_ifd(base_offset + next_offset, base_offset, ifd_t::TIFF);
        }
    } else {
        throw -1;
    }
}

void Tiffsniff::read_icc_profile(off_t offset) {
    char icc_header[128];
    if (fin->seekg(offset).good()) {
        if (!fin->read((char*)icc_header, 128).good()) {
            throw -1;
        }
        logger.debug("Device: %c%c%c%c\n", icc_header[12], icc_header[13], icc_header[14], icc_header[15]);
        logger.debug("Colour space: %c%c%c%c\n", icc_header[16], icc_header[17], icc_header[18], icc_header[19]);
        
        // looks like ICC profiles are stored in big endian format
        uint32_t icc_entries = icc_tag::read_uint32(fin);
        if (icc_entries == 0 || icc_entries > 32765 || !fin->good()) {
            throw -1;
        }
        
        bool found_trc = false;
        
        icc_tag tag;
        for (uint32_t i=0; i < icc_entries; i++) {
            tag.tag_signature = icc_tag::read_uint32(fin);
            tag.data_offset = icc_tag::read_uint32(fin);
            tag.element_size = icc_tag::read_uint32(fin);
            
            const off_t icc_tag_offset = offset + tag.data_offset;
            if (icc_tag_offset > file_size || !fin->good()) {
                throw -1;
            }
            
            auto fpos = fin->tellg(); // before we jump to an individual entry
            
            // Just grab the first TRC get find, since they are probably all the same anyway
            if (!found_trc && (tag.tag_signature == 0x67545243 || tag.tag_signature == 0x6B545243 ||
                tag.tag_signature == 0x62545243 || tag.tag_signature == 0x72545243)) {
                read_trc_entry(icc_tag_offset, tag.element_size);
                found_trc = true;
            }
            
            if (tag.tag_signature == 0x7258595A) { // red matrix column
                vector<double> col = read_xyztype_entry(icc_tag_offset, tag.element_size);
                luminance_weights[0] = col[1];
            }
            
            if (tag.tag_signature == 0x6758595a) { // green matrix column
                vector<double> col = read_xyztype_entry(icc_tag_offset, tag.element_size);
                luminance_weights[1] = col[1];
            }
            
            if (tag.tag_signature == 0x6258595A) { // blue matrix column
                vector<double> col = read_xyztype_entry(icc_tag_offset, tag.element_size);
                luminance_weights[2] = col[1];
            }
            
            if (tag.tag_signature == 0x77747074) { // media white point
                vector<double> wp = read_xyztype_entry(icc_tag_offset, tag.element_size);
                if (fabs(wp[1] - 1.0) < 1e-6) { // only if it looks like a valid illuminant
                    logger.debug("white point is %lg, %lg, %lg\n", wp[0], wp[1], wp[2]);
                }
            }
            
            fin->seekg(fpos); // return to ICC IFD
        }
        if (!found_trc) {
            logger.error("%s\n", "Warning: image contains ICC profile, but no TRC curve was found.");
            inferred_profile = profile_t::sRGB;
        } else {
            inferred_profile = profile_t::CUSTOM;
        }
    } else {
        throw -1;
    }
}

void Tiffsniff::read_trc_entry(off_t offset, uint32_t size) {
    if (fin->seekg(offset).good()) {
        char trc_type[4];
        if (fin->read((char*)trc_type, 4).good()) {
            if (strncmp(trc_type, "curv", 4) == 0) {
                read_curv_trc(offset, size);
            } else {
                if (strncmp(trc_type, "para", 4) == 0) {
                    read_para_trc(offset);
                }
            }
        } else {
            throw -1;
        }
    } else {
        throw -1;
    }
}

void Tiffsniff::read_curv_trc(off_t offset, uint32_t size) {
    if (fin->seekg(offset + 8).good()) {
        uint32_t ecount = icc_tag::read_uint32(fin);
        
        if ((ecount*2+12) > size || !fin->good()) {
            throw -1;
        }
        
        switch(ecount) {
        case 0: // nothing to do here
            break;
        case 1:
            gparm[0] = icc_tag::read_fixed8_8(fin);
            logger.debug("ICC 'curv' gamma is %lf\n", gparm[0]);
            break;
        default: // with two or more entries, use the table
            gtable = vector< pair<uint16_t, uint16_t> >(ecount);
            for (uint64_t i=0; i < ecount; i++) {
                gtable[i].first = i*65535/(ecount-1);
                gtable[i].second = icc_tag::read_uint16(fin);
            }
            logger.debug("ICC 'curv' gamma table with %ld entries\n", gtable.size());
            break;
        }
    } else {
        throw -1;
    }
}

void Tiffsniff::read_para_trc(off_t offset) {
    if (fin->seekg(offset + 8).good()) {
        uint16_t ftype = icc_tag::read_uint16(fin);
        icc_tag::read_uint16(fin); // dump the reserved bytes
        
        // See Table 70 of ICC 1:2001-12 (ICC profile v4)
        switch(ftype) {
        case 0:
            gparm[0] = icc_tag::read_fixed15_16(fin);
            gparm[1] = 1.0;
            break;
        case 1:
            gparm[0] = icc_tag::read_fixed15_16(fin);
            gparm[1] = icc_tag::read_fixed15_16(fin);
            gparm[2] = icc_tag::read_fixed15_16(fin);
            gparm[3] = 0;
            gparm[4] = -gparm[2]/gparm[1];
            break;
        case 2:
            gparm[0] = icc_tag::read_fixed15_16(fin);
            gparm[1] = icc_tag::read_fixed15_16(fin);
            gparm[2] = icc_tag::read_fixed15_16(fin);
            gparm[3] = 0;
            gparm[4] = -gparm[2]/gparm[1];
            gparm[6] = icc_tag::read_fixed15_16(fin);
            break;
        case 3: // sRGB, supposedly
            gparm[0] = icc_tag::read_fixed15_16(fin);
            gparm[1] = icc_tag::read_fixed15_16(fin);
            gparm[2] = icc_tag::read_fixed15_16(fin);
            gparm[3] = icc_tag::read_fixed15_16(fin);
            gparm[4] = icc_tag::read_fixed15_16(fin);
            break;
        case 4:
            gparm[0] = icc_tag::read_fixed15_16(fin);
            gparm[1] = icc_tag::read_fixed15_16(fin);
            gparm[2] = icc_tag::read_fixed15_16(fin);
            gparm[3] = icc_tag::read_fixed15_16(fin);
            gparm[4] = icc_tag::read_fixed15_16(fin);
            gparm[5] = icc_tag::read_fixed15_16(fin);
            gparm[6] = icc_tag::read_fixed15_16(fin);
            break;
        default: // treat this as linear
            logger.debug("unknown ICC parametricCurveType %d\n", ftype);
            throw -1;
            break;
        }
    } else {
        throw -1;
    }
}

vector<double> Tiffsniff::read_xyztype_entry(off_t offset, uint32_t) {
    if (fin->seekg(offset).good()) {
        uint32_t sig = icc_tag::read_uint32(fin);
        
        if (sig != 0x58595A20) {
            logger.debug("Expected XYZ signature, found %x instead.\n", sig);
            throw -1;
        }
        
        icc_tag::read_uint32(fin); // drop reserved field
        
        if (!fin->good()) {
            throw -1;
        }
        
        double x = icc_tag::read_fixed15_16(fin);
        double y = icc_tag::read_fixed15_16(fin);
        double z = icc_tag::read_fixed15_16(fin);
        
        return vector<double>{x,y,z};
    } else {
        throw -1;
    }
    return vector<double>{0,0,0};
}

double Tiffsniff::read_exif_gamma(off_t offset) {
    if (fin->seekg(offset).good()) {
        uint32_t top = icc_tag::read_uint32(fin);
        uint32_t bot = icc_tag::read_uint32(fin);
        
        return double(top)/double(bot);
    } else {
        throw -1;
    }
}

vector< pair<jpeg_app_t, off_t> > Tiffsniff::scan_jpeg_app_blocks(void) {
    vector< pair<jpeg_app_t, off_t> > blocks;
    fin->seekg(2);
    bool done = false;
    
    while (fin->good() && !done) {
        uint16_t app_id = icc_tag::read_uint16(fin);
        uint16_t bsize  = icc_tag::read_uint16(fin); // this size excludes the size bytes themselves
        
        if ((app_id & 0xff00) != 0xff00) {
            throw -1;
        }
        
        done = (app_id & 0xffe0) != 0xffe0;
        if (!done) {
            std::streamoff fpos = fin->tellg();
            
            unsigned char sig[128];
            int sn = 0;
            while (sn < 127 && (sig[sn] = fin->get()) != 0) sn++;
            
            if (strncasecmp((char*)sig, "ICC_PROFILE", 11) == 0) {
                blocks.push_back(make_pair(jpeg_app_t::ICC, off_t(fin->tellg()) + 2)); // skip over chunk numbers
            }
            if (strncasecmp((char*)sig, "EXIF", 4) == 0) {
                int i;
                for (i=0; i < 4 && fin->get() == 0; i++);
                blocks.push_back(make_pair(jpeg_app_t::EXIF, off_t(fin->tellg()) - 1));
            }
            
            if (!fin->seekg(fpos + off_t(bsize - 2)).good()) {
                throw -1;
            }
        }
    }
    return blocks;
}

void Tiffsniff::parse_png(off_t offset) {
    if (fin->seekg(offset).good()) {
        unsigned char magic[4];
        
        fin->seekg(8); // skip over rest of PNG signature
        uint32_t chunk_size = icc_tag::read_uint32(fin);
        if (!fin->read((char*)magic, 4).good()) throw -1;
        fin->seekg(chunk_size + 4, fin->cur);
        bool has_icc_profile = false;
        do {
            uint32_t chunk_size = icc_tag::read_uint32(fin);
            if (!fin->good() || !fin->read((char*)magic, 4).good()) {
                // not necessarily an error, maybe we already grabbed the profile ...
                return;
            }
            if (memcmp(magic,"gAMA", 4) == 0 && !has_icc_profile) {
                gparm[0] = 100000.0 / double(icc_tag::read_uint32(fin));
                inferred_profile = profile_t::CUSTOM;
                has_profile = true;
            }
            if (memcmp(magic,"sRGB", 4) == 0 && !has_icc_profile) { // not tested yet
                int psrgb = fin->get();
                if (psrgb >= 0 && psrgb <= 3) {
                    inferred_profile = profile_t::sRGB;
                    has_profile = true;
                }
            }
            #if mtfmapper_ZLIB_FOUND == 1
            if (memcmp(magic,"iCCP", 4) == 0) {
                unsigned char pfname[81];
                int name_len = 0;
                bool done = false;
                while (name_len < 80 && !done) {
                    pfname[name_len] = fin->get();
                    done = pfname[name_len] == 0;
                    name_len++;
                }
                if (fin->get() != 0) throw -1; // PNG compression method not deflate!
                
                size_t c_size = chunk_size - (name_len + 1);
                
                vector<unsigned char> c_chunk(chunk_size, 0);
                vector<unsigned char> u_chunk(chunk_size*10, 0);
                if (!fin->read((char*)c_chunk.data(), c_size).good()) throw -1;
                
                z_stream strm;
                strm.zalloc = Z_NULL;
                strm.zfree = Z_NULL;
                strm.opaque = Z_NULL;
                strm.avail_in = 0;
                strm.next_in = Z_NULL;
                int ret = inflateInit(&strm);
                if (ret != Z_OK) {
                    logger.error("%s\n", "PNG/ICC profile error: zlib init failed.");
                    throw -1;
                }
                
                strm.avail_in = c_size;
                strm.next_in = c_chunk.data();
                strm.avail_out = u_chunk.size();
                strm.next_out = u_chunk.data();
                ret = inflate(&strm, Z_NO_FLUSH);
                assert(ret != Z_STREAM_ERROR); 
                switch (ret) {
                case Z_NEED_DICT:
                    ret = Z_DATA_ERROR;    
                    [[fallthrough]];
                case Z_DATA_ERROR:
                case Z_MEM_ERROR:
                    (void)inflateEnd(&strm);
                    logger.error("%s\n", "PNG/ICC profile error: zlib inflate error.");
                    throw -1;
                }
                int u_size = u_chunk.size() - strm.avail_out;
                shared_ptr<std::iostream> fout(new std::stringstream());
                if (fout->good()) {
                    fout->write((char*)u_chunk.data(), u_size);
                    
                    shared_ptr<std::iostream> img_file = fin;
                    fin = fout;
                    read_icc_profile(0);
                    has_icc_profile = has_profile = true;
                    fin = img_file;
                }
            }
            #else
            logger.debug("%s\n", "No PNG/ICC support because zlib not found.");
            #endif
            fin->seekg(chunk_size + 4, fin->cur);
        } while (fin->good() && memcmp(magic, "IDAT", 4) == 0);
    } else {
        throw -1;
    }
}
