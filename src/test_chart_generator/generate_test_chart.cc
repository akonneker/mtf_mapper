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
#include "svg_page.h"
#include "svg_page_grid.h"
#include "svg_page_perspective.h"
#include "svg_page_lensgrid.h"
#include "svg_page_manualperspective.h"
#include "svg_page_focus.h"
#include "config.h"

#include <tclap/CmdLine.h>

#include <string>
using std::string;
using std::stringstream;

int main(int argc, char** argv) {

    stringstream ss;
    ss << mtfmapper_VERSION_MAJOR << "." << mtfmapper_VERSION_MINOR;
    
    vector<string> allowed_types;
    allowed_types.push_back("perspective");
    allowed_types.push_back("grid");
    allowed_types.push_back("halfgrid");
    allowed_types.push_back("thirdgrid");
    allowed_types.push_back("lensgrid");
    allowed_types.push_back("mfperspective");
    allowed_types.push_back("focus");
    TCLAP::ValuesConstraint<string> type_constraints(allowed_types);
    
    vector<string> allowed_sizes;
    allowed_sizes.push_back("a4");
    allowed_sizes.push_back("A4");
    allowed_sizes.push_back("a3");
    allowed_sizes.push_back("A3");
    allowed_sizes.push_back("a3s");
    allowed_sizes.push_back("A3S");
    allowed_sizes.push_back("a2");
    allowed_sizes.push_back("A2");
    allowed_sizes.push_back("a1");
    allowed_sizes.push_back("A1");
    allowed_sizes.push_back("a0");
    allowed_sizes.push_back("A0");
    TCLAP::ValuesConstraint<string> size_constraints(allowed_sizes);
    
    TCLAP::CmdLine cmd("Generate test charts for MTF50 measurements", ' ', ss.str());
    TCLAP::ValueArg<std::string> tc_type("t", "type", "Chart type (currently \"grid\" or \"perspective\")", false, "perspective", &type_constraints );
    cmd.add(tc_type);
    TCLAP::ValueArg<std::string> tc_size("s", "size", "Chart size (currently \"a4\" to \"a0\")", false, "a3", &size_constraints );
    cmd.add(tc_size);
    TCLAP::ValueArg<std::string> tc_ofname("o", "output", "Output file name (default chart.svg)", false, "chart.svg", "filename" );
    cmd.add(tc_ofname);
    TCLAP::ValueArg<double> tc_distance("d", "distance", "Distance from test chart (mm)", false, 1500, "distance", cmd);
    TCLAP::SwitchArg tc_alpa("","alpa-scale","Add Alpa-specific distance scale to focus chart type", cmd, false);
    
    cmd.parse(argc, argv);
    
    printf("Chart type: %s at %s size\n", tc_type.getValue().c_str(), tc_size.getValue().c_str());
    
    if (tc_type.getValue().compare("perspective") == 0) {
        Svg_page_perspective p(tc_size.getValue(), tc_ofname.getValue());
        printf("Generating chart for viewing distance of %.0lf mm\n", tc_distance.getValue());
        p.set_viewing_parameters(tc_distance.getValue(), -45/180.0*M_PI); // must still become paramers
        p.render();
    } else {
        grid_scale scale = INVALIDGRID;
        if (tc_type.getValue().compare("grid") == 0) {
            scale = FULLGRID;
        }
        if (tc_type.getValue().compare("halfgrid") == 0) {
            scale = HALFGRID;
        }
        if (tc_type.getValue().compare("thirdgrid") == 0) {
            scale = THIRDGRID;
        }
        
        if (scale != INVALIDGRID) {
            Svg_page_grid p(tc_size.getValue(), tc_ofname.getValue(), scale);
            p.render();
        } else {
            if (tc_type.getValue().compare("lensgrid") == 0) {
                Svg_page_lensgrid p(tc_size.getValue(), tc_ofname.getValue());
                p.render();
                return 0;
            } 
            
            if (tc_type.getValue().compare("mfperspective") == 0) {
                if (tc_size.getValue() == "A3" || tc_size.getValue() == "a3"){
                    Svg_page_manualperspective p(tc_size.getValue(), tc_ofname.getValue());
                    p.set_viewing_parameters(tc_distance.getValue(), -45/180.0*M_PI);
                    p.render();
                } else {
                    printf("Sorry, only A3 size supported for this chart type at the moment\n");
                    exit(-1);
                }
            } else {
                if (tc_type.getValue().compare("focus") == 0) {
                    Svg_page_focus p(tc_size.getValue(), tc_ofname.getValue(), tc_alpa.getValue());
                    p.set_viewing_parameters(tc_distance.getValue(), -45/180.0*M_PI);      
                    p.render(); 
                } else {
                    printf("Illegal chart type %s\n", tc_type.getValue().c_str());
                }
            }
        }
    } 
    
    return 0;
}
