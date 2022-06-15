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

#include "include/ca_renderer_print.h"

void Ca_renderer_print::render(const vector<Block>& blocks) {
    FILE* fout = fopen(ofname.c_str(), "wt");
    
    // Only Meridional direction, for now
    fprintf(fout, "# MTF Mapper CA, output format version 2\n");
    fprintf(fout, "# column  1: block_id\n");
    fprintf(fout, "# column  2: edge centroid x (pixels)\n");
    fprintf(fout, "# column  3: edge centroid y (pixels)\n");
    fprintf(fout, "# column  4: Red channel shift relative to Green channel (pixels)\n");
    fprintf(fout, "# column  5: Blue channel shift relative to Green channel (pixels)\n");
    fprintf(fout, "# column  6: Red channel shift relative to Green channel (%% of radial distance)\n");
    fprintf(fout, "# column  7: Blue channel shift relative to Green channel (%% of radial distance)\n");
    for (size_t i=0; i < blocks.size(); i++) {
        if (!blocks[i].valid) continue;
        
        for (size_t k=0; k < 4; k++) {
            Point2d cent = blocks[i].get_edge_centroid(k);
            Point2d dir = cent - img_centre;
            double centre_dist = std::max(1.0, norm(dir));
            dir = dir * (1.0/norm(dir));
            double delta = dir.dot(blocks[i].get_normal(k));
            
            if (!blocks[i].get_edge_valid(k)) continue;
            if (fabs(blocks[i].get_ca(k).x - Edge_info::nodata) < 1e-6 || fabs(blocks[i].get_ca(k).y - Edge_info::nodata) < 1e-6) continue;

            if (fabs(delta) >= cos(45.0/180*M_PI) || blocks.size() == 1 || allow_all_edges) {
                fprintf(fout, "%d %lf %lf %lf %lf %lf %lf\n", 
                    int(i),
                    blocks[i].get_edge_centroid(k).x, blocks[i].get_edge_centroid(k).y,
                    blocks[i].get_ca(k).x, blocks[i].get_ca(k).y,
                    100*blocks[i].get_ca(k).x / centre_dist, 100*blocks[i].get_ca(k).y / centre_dist
                );
            }
        }
    }
    
    fclose(fout);
}