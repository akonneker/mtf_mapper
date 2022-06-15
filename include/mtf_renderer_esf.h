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
#ifndef MTF_RENDERER_ESF_H
#define MTF_RENDERER_ESF_H

#include "mtf_renderer.h"
#include "include/output_version.h"
#include "common_types.h"

class Mtf_renderer_esf : public Mtf_renderer {
  public:
    Mtf_renderer_esf(const std::string& fname_esf, const std::string& fname_lsf, Output_version::type output_version)
      :  ofname_esf(fname_esf), ofname_lsf(fname_lsf), output_version(output_version) {
      
    }
    
    void render(const vector<Block>& blocks) {
        FILE* fout_esf = fopen(ofname_esf.c_str(), "wt");
        FILE* fout_lsf = fopen(ofname_lsf.c_str(), "wt");

        printf("output version = %d", (int)output_version);

        if (output_version >= Output_version::V2) {
            fprintf(fout_esf, "# MTF Mapper ESF, output format version %d \n", int(output_version));
            fprintf(fout_esf, "# column  1: block_id\n");
            fprintf(fout_esf, "# column  2: edge centroid x (pixels)\n");
            fprintf(fout_esf, "# column  3: edge centroid y (pixels)\n");
            fprintf(fout_esf, "# column  4: ESF sample spacing (pixels)\n");
            fprintf(fout_esf, "# column  5: number of ESF samples, N\n");
            fprintf(fout_esf, "# next N columns: ESF samples\n");

            fprintf(fout_lsf, "# MTF Mapper LSF, output format version %d \n", int(output_version));
            fprintf(fout_lsf, "# column  1: block_id\n");
            fprintf(fout_lsf, "# column  2: edge centroid x (pixels)\n");
            fprintf(fout_lsf, "# column  3: edge centroid y (pixels)\n");
            fprintf(fout_lsf, "# column  4: LSF sample spacing (pixels)\n");
            fprintf(fout_lsf, "# column  5: number of LSF samples, N\n");
            fprintf(fout_lsf, "# next N columns: LSF samples\n");
        }

        for (size_t i=0; i < blocks.size(); i++) {
            for (size_t k=0; k < 4; k++) {
                const vector<double>& esf = blocks[i].get_esf(k);

                if (output_version >= Output_version::V2) {
                    fprintf(fout_esf, "%d %lf %lf 0.125 %d ", int(i), blocks[i].get_edge_centroid(k).x, blocks[i].get_edge_centroid(k).y, int(esf.size()));
                }

                for (size_t j=0; j < esf.size(); j++) {
                    fprintf(fout_esf, "%lf ", esf[j]);
                }
                fprintf(fout_esf, "\n");
                
                double sum = 0;
                for (size_t j=1; j < esf.size()-1; j++) {
                    sum += (esf[j+1] - esf[j-1])*0.5;
                }

                if (output_version >= Output_version::V2) {
                    fprintf(fout_lsf, "%d %lf %lf 0.125 %d ", int(i), blocks[i].get_edge_centroid(k).x, blocks[i].get_edge_centroid(k).y, int(esf.size()));
                }
                
                double sign = sum < 0 ? -1 : 1;
                fprintf(fout_lsf, "%lf ", 0.0);
                for (size_t j=1; j < esf.size()-1; j++) {
                    fprintf(fout_lsf, "%lf ", sign*(esf[j+1] - esf[j-1])*0.5);
                }
                fprintf(fout_lsf, " %lf\n", 0.0);
            }
        }    
        fclose(fout_esf);
        fclose(fout_lsf);
    }
    
    string ofname_esf;
    string ofname_lsf;
    Output_version::type output_version;

};

#endif
