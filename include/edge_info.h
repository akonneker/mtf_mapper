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
#ifndef EDGE_INFO_H
#define EDGE_INFO_H

#include <cstring>
#include <string>
using std::string;

#include <opencv2/core/core.hpp>

#include "include/job_metadata.h"
#include "include/logger.h"
#include <memory>

class Edge_info {
  public:
    Edge_info(void) {
    }
    
    static bool serialize_header(FILE* fout, size_t valid_count, const Job_metadata& metadata) {
        vector<char> buffer(5*sizeof(uint32_t) + 2*sizeof(double));
        
        char* obuf = buffer.data();
        write_uint32(&obuf, valid_count);
        write_double(&obuf, metadata.pixel_pitch);
        write_double(&obuf, metadata.mtf_contrast);
        write_uint32(&obuf, uint32_t(metadata.bayer));
        write_uint32(&obuf, uint32_t(metadata.channels));
        write_uint32(&obuf, uint32_t(metadata.image_width));
        write_uint32(&obuf, uint32_t(metadata.image_height));
        
        size_t nwritten = fwrite(buffer.data(), 1, buffer.size(), fout);
        
        if (nwritten != buffer.size()) {
            logger.error("Could not write header to edge info serialization file, tried %ld bytes, only wrote %ld\n", buffer.size(), nwritten);
            return false;
        }
        
        return true;
    }
    
    static bool deserialize_header(FILE* fin, size_t& valid_count, Job_metadata& metadata) {
        vector<char> buffer(5*sizeof(uint32_t) + 2*sizeof(double));
        
        size_t nread = fread(buffer.data(), 1, buffer.size(), fin);
        
        if (nread != buffer.size()) {
            logger.error("Could not read header from edge info serialization file, tried %ld bytes, only wrote %ld\n", buffer.size(), nread);
            return false;
        }
        
        char* ibuf = buffer.data();
        valid_count = read_uint32(&ibuf);
        metadata.pixel_pitch = read_double(&ibuf);
        metadata.mtf_contrast = read_double(&ibuf);
        metadata.bayer = Bayer::bayer_t(read_uint32(&ibuf));
        metadata.channels = read_uint32(&ibuf);
        metadata.image_width = read_uint32(&ibuf);
        metadata.image_height = read_uint32(&ibuf);

        return true;
    }
    
    // serialize multiple edges into a string
    static bool serialize(
        FILE* fout,
        const vector<cv::Point2d>& in_centroid,
        const vector<double>& in_angle,
        const vector<double>& in_mtf50,
        const vector<double>& /*in_quality*/,
        const vector<std::shared_ptr<vector<double>>>& in_sfr,
        const vector<std::shared_ptr<vector<double>>>& in_esf,
        const vector<cv::Point2d> in_snr,
        const vector<cv::Point2d> in_chromatic_aberration,
        const vector<bool> valid_edge,
        const vector<double>& in_edge_length,
        const vector<double>& in_geometric_length
    ) {
        // there is no need to keep data grouped by block at this point, so each record is independent
        vector<char> buffer(20*1024);
        for (size_t k=0; k < in_centroid.size(); k++) {
            if (!valid_edge[k]) continue;
            
            char* buf = buffer.data();
            
            write_double(&buf, in_centroid[k].x);
            write_double(&buf, in_centroid[k].y);
            write_double(&buf, in_angle[k]); // whole angle, not relative
            write_double(&buf, in_mtf50[k]); // MTF-XX, in c/p
            write_double(&buf, in_snr[k].x); // mean CNR
            write_double(&buf, in_snr[k].y); // oversampling factor
            write_double(&buf, in_chromatic_aberration[k].x); // CA red-green
            write_double(&buf, in_chromatic_aberration[k].y); // CA blue-green
            write_double(&buf, in_edge_length[k]); // edge length (effective ROI), in pixels
            write_double(&buf, in_geometric_length[k]); // geometric edge length, in pixels

            write_uint32(&buf, in_sfr[k]->size());
            for (size_t j=0; j < in_sfr[k]->size(); j++) {
                write_float(&buf, (float)(*in_sfr[k])[j]);
            }
            
            write_uint32(&buf, in_esf[k]->size());
            for (size_t j=0; j < in_esf[k]->size(); j++) {
                write_float(&buf, (float)(*in_esf[k])[j]);
            }
            
            size_t size_in_bytes = buf - buffer.data();
            size_t nwritten = fwrite(buffer.data(), 1, size_in_bytes, fout);
            if (nwritten != size_in_bytes) {
                logger.error("Could not write to edge info serialization file, tried %ld bytes, only wrote %ld\n", buffer.size(), size_in_bytes);
                return false;
            }
        }
        return true;
    }
    
    // deserialize a single edge from string s
    static Edge_info deserialize(FILE* fin, bool& valid) {
        Edge_info b;
        valid = false;
        
        vector<char> buffer(20*1024);
        
        char* ibuf = buffer.data();
        
        size_t nread = fread(ibuf, 1, 10*sizeof(double), fin);
        if (nread != 10*sizeof(double)) {
            logger.error("Could not read from edge info serialization file, tried %ld bytes, only read %ld\n", 9*sizeof(double), nread);
            return b;
        }
        
        b.centroid.x = read_double(&ibuf);
        b.centroid.y = read_double(&ibuf);
        b.angle = read_double(&ibuf);
        b.mtf50 = read_double(&ibuf);
        b.snr.x = read_double(&ibuf);
        b.snr.y = read_double(&ibuf);
        b.chromatic_aberration.x = read_double(&ibuf);
        b.chromatic_aberration.y = read_double(&ibuf);
        b.edge_length = read_double(&ibuf);
        b.geometric_edge_length = read_double(&ibuf);

        nread = fread(ibuf, 1, sizeof(uint32_t), fin);
        if (nread != sizeof(uint32_t)) {
            logger.error("%s\n", "Could not read from edge info serialization file, SFR size read failed");
            return b;
        }
        
        size_t nsfr = read_uint32(&ibuf);
        if (nsfr > 128) {
            logger.error("Unexpected SFR length (%ld) in edge deserialization, aborting\n", nsfr);
            return b;
        }
        
        nread = fread(ibuf, 1, sizeof(float)*nsfr, fin);
        if (nread != sizeof(float)*nsfr) {
            logger.error("%s\n", "Could not read from edge info serialization file, SFR read failed");
            return b;
        }
        
        b.sfr = std::shared_ptr<vector<double>>(new vector<double>(nsfr, 0.0));
        for (size_t j=0; j < nsfr; j++) {
            (*b.sfr)[j] = read_float(&ibuf);
        }
        
        nread = fread(ibuf, 1, sizeof(uint32_t), fin);
        if (nread != sizeof(uint32_t)) {
            logger.error("%s\n", "Could not read from edge info serialization file, ESF size read failed");
            return b;
        }
        
        size_t nesf = read_uint32(&ibuf);
        if (nsfr > 256) {
            logger.error("%s\n", "Unexpected ESF length in edge deserialization, aborting");
            return b;
        }
        
        nread = fread(ibuf, 1, sizeof(float)*nesf, fin);
        if (nread != sizeof(float)*nesf) {
            logger.error("%s\n", "Could not read from edge info serialization file, ESF read failed");
            return b;
        }
        
        b.esf = std::shared_ptr<vector<double>>(new vector<double>(nesf, 0.0));
        for (size_t j=0; j < nesf; j++) {
            (*b.esf)[j] = read_float(&ibuf);
        }
        
        valid = true;
        return b;
    }
    
    void set_metadata(const Job_metadata& md) {
        metadata = md;
    }
    
    cv::Point2d centroid;
    double angle;
    double mtf50;
    double quality;
    std::shared_ptr<vector<double>> sfr;
    std::shared_ptr<vector<double>> esf;
    cv::Point2d snr;
    cv::Point2d chromatic_aberration;
    Job_metadata metadata;
    double edge_length = 0;
    double geometric_edge_length = 0;
    
    static constexpr double nodata = -1073741824;
    
  private:
    void static write_uint32(char** buffer, uint32_t val) {
        memcpy(*buffer, (char*)&val, sizeof(uint32_t));
        *buffer += sizeof(uint32_t);
    }
    
    void static write_double(char** buffer, double val) {
        memcpy(*buffer, (char*)&val, sizeof(double));
        *buffer += sizeof(double);
    }
    
    void static write_float(char** buffer, float val) {
        memcpy(*buffer, (char*)&val, sizeof(float));
        *buffer += sizeof(float);
    }
    
    uint32_t static read_uint32(char** buffer) {
        uint32_t val = 0;
        memcpy((char*)&val, *buffer, sizeof(uint32_t));
        *buffer += sizeof(uint32_t);
        return val;
    }
    
    double static read_double(char** buffer) {
        double val = 0;
        memcpy((char*)&val, *buffer, sizeof(double));
        *buffer += sizeof(double);
        return val;
    }
    
    float static read_float(char** buffer) {
        float val = 0;
        memcpy((char*)&val, *buffer, sizeof(float));
        *buffer += sizeof(float);
        return val;
    }
    
};

#endif
