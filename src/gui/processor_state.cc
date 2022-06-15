/*
Copyright 2021 Frans van den Bergh. All rights reserved.

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

#include "processor_state.h"

QString Processor_state::updated_arguments(void) const {
    QString arguments(base_arguments);
    
    if (single_roi_mode) {
        if (!arguments.contains("--annotate") && !arguments.contains("-a ")) {
            arguments = arguments + QString(" -a");
        }
        if (!arguments.contains("--edges") && !arguments.contains("-q ")) {
            arguments = arguments + QString(" -q");
        }
        if (!arguments.contains("--single-roi")) {
            arguments = arguments + QString(" --single-roi");
        }
        
        // disable all the other mutually exclusive outputs
        arguments.replace("--lensprofile ", " ");
        arguments.replace("-s ", " ");
        arguments.replace("--surface ", " ");
        arguments.replace("-p ", " ");
        arguments.replace("--profile ", " ");
        arguments.replace("--focus ", " ");
        arguments.replace("--chart-orientation ", " ");
    }
    
    if (focus_mode) {
        if (!arguments.contains("--focus")) {
            arguments = arguments + QString(" --focus");
        }
        // disable all the other mutually exclusive outputs
        arguments.replace("--lensprofile ", " ");
        arguments.replace("-q ", " ");
        arguments.replace("--edges ", " ");
        arguments.replace("-a ", " ");
        arguments.replace("--annotate ", " ");
        arguments.replace("-s ", " ");
        arguments.replace("--surface ", " ");
        arguments.replace("-p ", " ");
        arguments.replace("--profile ", " ");
    }
    
    
    if (imatest_mode) {
        if (!arguments.contains("--imatest-chart")) {
            arguments = arguments + QString(" --imatest-chart");
        }
    }
    return arguments;
}

