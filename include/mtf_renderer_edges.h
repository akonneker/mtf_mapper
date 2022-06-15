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
#ifndef MTF_RENDERER_EDGES_H
#define MTF_RENDERER_EDGES_H

#include "mtf_renderer.h"
#include "common_types.h"
#include "include/ordered_point.h"
#include "include/output_version.h"
#include "include/job_metadata.h"

class Mtf_renderer_edges : public Mtf_renderer {
  public:
    Mtf_renderer_edges(const std::string& fname, 
      const std::string& sfrname,
      const std::string& devname,
      const std::string& serialname,
      Output_version::type output_version,
      const Job_metadata& metadata,
      bool lpmm_mode=false, double pixel_size=1.0) 
      :  ofname(fname), sfrname(sfrname), devname(devname), serialname(serialname),
         lpmm_mode(lpmm_mode), pixel_size(pixel_size),
         output_version(output_version), metadata(metadata) {
         
         // TODO: this is a bit hacky; we should have a separate flag for units ?
         this->metadata.pixel_pitch = fabs(pixel_size - 1) < 1e-6 ? 1000.0 : 1000.0/pixel_size;
      
    }
    
    void render(const vector<Block>& blocks) {
        Point2d centr(0,0);
        for (size_t i=0; i < blocks.size(); i++) {
            centr += blocks[i].get_centroid();
        }
        centr = centr*(1.0/double(blocks.size()));
        
    
    
        FILE* fout = fopen(ofname.c_str(), "wt");
        FILE* sfrout = fopen(sfrname.c_str(), "wt");
        FILE* devout = fopen(devname.c_str(), "wt");
        
        FILE* serout = NULL;
        if (output_version >= Output_version::V2) {
            serout = fopen(serialname.c_str(), "wb");
            
            if (serout) {
                size_t valid_count = 0;
                for (size_t i=0; i < blocks.size(); i++) {
                    for (size_t k=0; k < 4; k++) {
                        valid_count += blocks[i].get_edge_valid(k);
                    }
                }
                
                Edge_info::serialize_header(
                    serout, 
                    valid_count,
                    metadata
                );
            }
        }
        
        if (output_version >= Output_version::V2) {
            fprintf(sfrout, "# MTF Mapper SFR, output format version %d \n", int(output_version));
            fprintf(sfrout, "# column  1: block_id\n");
            fprintf(sfrout, "# column  2: edge centroid x (pixels)\n");
            fprintf(sfrout, "# column  3: edge centroid y (pixels)\n");
            fprintf(sfrout, "# column  4: slanted edge orientation (degrees, modulo 45 degrees)\n");
            fprintf(sfrout, "# column  5: edge orientation relative to radial line to centre (degrees)\n");
            fprintf(sfrout, "# column  6: mean CNR (50%% weight to each of dark and bright sides)\n");
            fprintf(sfrout, "# column  7: dark side CNR\n");
            fprintf(sfrout, "# column  8: bright side CNR\n");
            fprintf(sfrout, "# column  9: dark side SNR\n");
            fprintf(sfrout, "# column 10: bright side SNR\n");
            fprintf(sfrout, "# column 11: contrast\n");
            fprintf(sfrout, "# column 12: effective oversampling factor (maximum is 8x, below 4x is considered poor)\n");
            fprintf(sfrout, "# column 13: effective edge length (a value below 25 pixels is considered poor)\n");
            
            fprintf(fout, "# MTF Mapper MTF summary, output format version %d \n", int(output_version));
            fprintf(fout, "# column  1: block_id\n");
            fprintf(fout, "# column  2: edge centroid x (pixels)\n");
            fprintf(fout, "# column  3: edge centroid y (pixels)\n");
            fprintf(fout, "# column  4: MTF-%d value (%s)\n", (int)lrint(metadata.mtf_contrast*100), lpmm_mode ? "lp/mm" : "c/p");
            fprintf(fout, "# column  5: nearby corner x (pixels)\n");
            fprintf(fout, "# column  6: nearby corner y (pixels)\n");
            fprintf(fout, "# column  7: mean CNR (50%% weight to each of dark and bright sides)\n");
            fprintf(fout, "# column  8: effective oversampling factor (maximum is 8x, below 4x is considered poor)\n");
            fprintf(fout, "# column  9: effective edge length (a value below 25 pixels is considered poor)\n");
            
            fprintf(devout, "# MTF Mapper line deviation, output format version %d \n", int(output_version));
            fprintf(devout, "# column  1: block_id\n");
            fprintf(devout, "# column  2: edge centroid x (pixels)\n");
            fprintf(devout, "# column  3: edge centroid y (pixels)\n");
            fprintf(devout, "# column  4: slope, i.e., rise/run\n");
            fprintf(devout, "# column  5: rise\n");
            fprintf(devout, "# column  6: run\n");
        }
        
        vector<int> corder(4);
        vector<int> eorder(4);
        bool serialization_errors = false;
        for (size_t i=0; i < blocks.size(); i++) {
        
            vector<Ordered_point> corners(4);
            // start with first corner that is left of 12:00
            for (size_t k=0; k < 4; k++) {
                Point2d dir = blocks[i].get_corner(k) - blocks[i].get_centroid();
                dir *= 1.0/norm(dir);
                corners[k].second = k;
                corners[k].first = fmod(atan2(dir.y, dir.x) + M_PI*0.5, 2.0*M_PI);
                if (corners[k].first < 0) {
                    corners[k].first += 2.0*M_PI;
                }
            }
            sort(corners.begin(), corners.end());
            for (size_t k=0; k < 4; k++) {
                corder[k] = int(corners[k].second);
            }
            
            // just a sanity check for debugging
            bool all_found = true;
            for (int k=0; k < 4; k++) {
                bool found_k = false;
                for (int j=0; j < 4; j++) {
                    found_k |= corder[j] == k;
                }
                all_found &= found_k;
            }
            if (!all_found) {
                logger.debug("Warning: not all corners are present in corder: [%d %d %d %d]\n", corder[0], corder[1], corder[2], corder[3]);
            }
            
            if (output_version >= Output_version::V2) {
                // assign edges to corners so that centroid falls on edge between corner and next
                for (size_t k=0; k < 4; k++) {
                    Point2d corn_dir = blocks[i].get_corner(corder[k]) - blocks[i].get_corner(corder[(k+1)%4]);
                    corn_dir *= 1.0/norm(corn_dir);
                    Point2d corn_n(-corn_dir.y, corn_dir.x);
                    double min_dist = std::numeric_limits<double>::max();
                    size_t min_edge = 0;
                    for (size_t j=0; j < 4; j++) {
                        double perp_dist = fabs((blocks[i].get_edge_centroid(j) - blocks[i].get_corner(corder[k])).dot(corn_n));
                        if (perp_dist < min_dist) {
                            min_edge = j;
                            min_dist = perp_dist;
                        }
                    }
                    eorder[k] = min_edge;
                }
            } else {
                // the old way of associating an edge with a corner; not really sure what I intended here
                // but I am keeping it for backwards compatibility
                for (size_t k=0; k < 4; k++) {
                    Point2d cdir = blocks[i].get_corner(corder[k]) - blocks[i].get_centroid();
                    cdir *= 1.0/norm(cdir);
                    corners[k].second = k;
                    corners[k].first = fmod(atan2(cdir.y, cdir.x) + M_PI*0.5, 2.0*M_PI);
                    if (corners[k].first < 0) {
                        corners[k].first += 2.0*M_PI;
                    }
                }
                sort(corners.begin(), corners.end());
                for (size_t k=0; k < 4; k++) {
                    eorder[k] = int(corners[k].second);
                }
            }
            
            // just a sanity check for debugging
            all_found = true;
            for (int k=0; k < 4; k++) {
                bool found_k = false;
                for (int j=0; j < 4; j++) {
                    found_k |= eorder[j] == k;
                }
                all_found &= found_k;
            }
            if (!all_found) {
                logger.debug("Warning: not all edges are present in eorder: [%d %d %d %d]\n", eorder[0], eorder[1], eorder[2], eorder[3]);
            }
            
            for (size_t k=0; k < 4; k++) {
                int j = corder[k];
                int l = eorder[k];
                double val = blocks[i].get_mtf50_value(l);
                Point2d ec = blocks[i].get_edge_centroid(l);
                Point2d cr = blocks[i].get_corner(j);
                const Snr& snr = blocks[i].get_snr(l);
                double edge_length = blocks[i].get_edge_length(l);
                
                // in later output format versions, skip the empty rows
                if (output_version >= Output_version::V2 &&
                    blocks[i].get_sfr(l)[0] == 0) {
                    continue;
                }
                
                if (output_version == Output_version::V1) {
                    fprintf(fout, "%d %lf %lf %lf %lf %lf\n",
                        int(i),
                        ec.x, ec.y,
                        lpmm_mode ? val*pixel_size : val,
                        cr.x, cr.y
                    );
                }
                
                if (output_version >= Output_version::V2) {
                    fprintf(fout, "%d %lf %lf %lf %lf %lf %.3lf %.1lf %.1lf\n",
                        int(i),
                        ec.x, ec.y,
                        lpmm_mode ? val*pixel_size : val,
                        cr.x, cr.y, snr.mean_cnr(), snr.oversampling(),
                        edge_length
                    );
                }
                
                fprintf(sfrout, "%d %lf %lf ", int(i), ec.x, ec.y);
                
                double edge_angle = atan2(-blocks[i].get_normal(l).x, blocks[i].get_normal(l).y);
                fprintf(sfrout, "%lf ", angle_reduce(edge_angle));
                
                Point2d cent = blocks[i].get_edge_centroid(l);

                Point2d dir = cent - centr;
                if (fabs(norm(dir)) > 1e-12) {
                    dir = dir * (1.0 / norm(dir));
                }
                
                Point2d norm = blocks[i].get_normal(l);

                double delta = dir.x*norm.x + dir.y*norm.y;
                fprintf(sfrout, "%lf ", acos(fabs(delta))/M_PI*180.0);
                
                if (output_version >= Output_version::V2) {
                    fprintf(sfrout, "%.3lf %.3lf %.3lf %.3lf %.3lf %.1lf %.1lf %.1lf ",
                        snr.mean_cnr(), snr.dark_cnr(), snr.bright_cnr(),
                        snr.dark_snr(), snr.bright_snr(), snr.contrast(),
                        snr.oversampling(), edge_length
                    );
                }
                
                const vector<double>& sfr = blocks[i].get_sfr(l);
                for (size_t j=0; j < sfr.size(); j++) {
                    fprintf(sfrout, "%lf ", sfr[j]);
                }
                fprintf(sfrout, "\n");
                
                cv::Point3d deviation = blocks[i].get_line_deviation(l);
                fprintf(devout, "%d %lf %lf %.8lf %lf %lf\n", int(i), ec.x, ec.y, deviation.x, deviation.y, deviation.z);
            }
            
            if (output_version >= Output_version::V2) {
                bool success = blocks[i].serialize(serout);
                serialization_errors |= !success;
            }
        }    
        fclose(fout);
        fclose(sfrout);
        fclose(devout);
        
        if (output_version >= Output_version::V2 && serout) {
            fclose(serout);
        }
        
        if (serialization_errors) {
            logger.error("%s\n", "One or more errors during edge serialization");
        }
    }
    
  private:
    double angle_reduce(double x) {
        double quad1 = fabs(fmod(x, M_PI/2.0));
        if (quad1 > M_PI/4.0) {
            quad1 = M_PI/2.0 - quad1;
        }
        quad1 = quad1 / M_PI * 180;
        return quad1;
    }
    
    string ofname;
    string sfrname;
    string devname;
    string serialname;
    bool    lpmm_mode;
    double  pixel_size;
    Output_version::type output_version;
    Job_metadata metadata;
};

#endif
