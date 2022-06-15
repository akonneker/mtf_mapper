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

#ifndef PROCESSOR_STATE_H
#define PROCESSOR_STATE_H

#include <QString>
#include <QStringList>


#include <memory>

#include "processing_command.h"
#include "raw_developer.h"

class Processor_state {
  public:
    Processor_state(void) {}
    
    Processor_state(const QStringList& input_files,
      const QString& gnuplot_binary, const QString& exiv2_binary,
      const QString& arguments,
      std::shared_ptr<Raw_developer> raw_developer)
    : input_files(input_files), gnuplot_binary(gnuplot_binary), 
      exiv2_binary(exiv2_binary), base_arguments(arguments),
      raw_developer(raw_developer) {}
    
    void set_single_roi_mode(bool roi) {
        single_roi_mode = roi;
    }
    
    void set_focus_mode(bool focus) {
        focus_mode = focus;
    }
    
    void set_imatest_mode(bool imatest) {
        imatest_mode = imatest;
    }
    
    void set_manual_roi_mode(bool manual_roi) {
        manual_roi_mode = manual_roi;
    }
    
    Processing_command::state_t initial_state(void) const {
        return manual_roi_mode ? Processing_command::state_t::AWAIT_ROI : Processing_command::state_t::READY;
    }
    
    bool is_valid(void) const {
        return input_files.size() > 0;
    }
    
    QString updated_arguments(void) const;
    
    bool single_roi_mode = false;
    bool focus_mode = false;
    bool imatest_mode = false;
    bool manual_roi_mode = false;
    
    QStringList input_files;
    QString gnuplot_binary;
    QString exiv2_binary;
    QString base_arguments;
    std::shared_ptr<Raw_developer> raw_developer;
};

#endif
