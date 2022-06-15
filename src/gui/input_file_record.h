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

#ifndef INPUT_FILE_RECORD_H
#define INPUT_FILE_RECORD_H

#include <QString>

#include "raw_developer.h"

class Input_file_record {
  public:
    enum class state_t {
        SUBMITTED,
        COMPLETED,
        FAILED
    };
    
    Input_file_record(void) {}
    
    Input_file_record(const QString& input_fname, 
        const QString& args, const QString& temp_dir,
        std::shared_ptr<Raw_developer> raw_developer) 
      : input_fname(input_fname), arguments(args),
        temp_dir(temp_dir), state(state_t::SUBMITTED),
        raw_developer(raw_developer) {}
    
    QString get_input_name(void) const {
        return input_fname;
    }
    
    QString get_output_name(void) const {
        return output_fname;
    }
    
    void set_output_name(const QString& fname) {
        output_fname = fname;
    }
    
    QString get_arguments(void) const {
        return arguments;
    }
    
    QString get_temp_dir(void) const {
        return temp_dir;
    }
    
    state_t get_state(void) const {
        return state;
    }
    
    void set_state(state_t new_state) {
        state = new_state;
    }
    
    std::shared_ptr<Raw_developer> get_raw_developer(void) const {
        return raw_developer;
    }
    
    bool is_valid(void) const {
        return input_fname.length() > 0;
    }
    
    QString input_fname;
    QString arguments;
    QString temp_dir;
    QString output_fname;
    state_t state;
    std::shared_ptr<Raw_developer> raw_developer;
};

#endif
