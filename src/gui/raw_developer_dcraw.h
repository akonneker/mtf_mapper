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
#ifndef RAW_DEVELOPER_DCRAW_H
#define RAW_DEVELOPER_DCRAW_H

#include "raw_developer.h"
#include <vector>
using std::vector;

class Raw_developer_dcraw : public Raw_developer {
  public:
    Raw_developer_dcraw(const QString& dcraw_exec) : dcraw_exec(dcraw_exec) {}
    virtual ~Raw_developer_dcraw(void) {}
    
  private:
    virtual QString program_name([[maybe_unused]] bool bayer_mode) {
        return dcraw_exec;
    }
    
    virtual QStringList arguments(bool bayer_mode, const QString& input, [[maybe_unused]] const QString& output) {
        QStringList args;
        
        if (bayer_mode) {
            args << "-4" << "-T" << "-D" << "-c" << input;
        } else {
            args << "-w" <<  "-4" << "-T" << "-q" << "3" << "-H" << "1" << "-c" << input;
        }
        
        return args;
    }
    
    virtual bool accepts_suffix(const QString& suffix) {
        bool found = false;
        for (size_t i=0; !found && i < suffixes.size(); i++) {
            found = suffix.compare(suffixes[i], Qt::CaseInsensitive) == 0;
        }
        return found;
    };
  
    QString dcraw_exec;
    vector<QString> suffixes = {
      "NEF",  // Nikon
      "ARW",  // Sony
      "PEF",  // Pentax
      "IIQ",  // Phase One
      "MOS",  // Leaf
      "ORF",  // Olympus
      "RW2",  // Panasonic
      "RAF",  // Fujifilm -> bayer mode will probably break horribly
      "DNG",  // Pentax/Ricoh, maybe others
      "CR2"  // Canon
    };
};

#endif

