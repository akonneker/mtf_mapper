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
#ifndef RAW_DEVELOPER_LIBRAW_H
#define RAW_DEVELOPER_LIBRAW_H

#include "raw_developer.h"
#include <QDir>
#include <QFileInfo>

class Raw_developer_libraw : public Raw_developer {
  public:
    Raw_developer_libraw(const QString& dcraw_emu_exec, const QString& unprocessed_raw_exec)
    : dcraw_emu_exec(dcraw_emu_exec), unprocessed_raw_exec(unprocessed_raw_exec) {}
    
    virtual ~Raw_developer_libraw(void) {}
    
  private:
    virtual QString program_name(bool bayer_mode) {
        if (bayer_mode) {
            return unprocessed_raw_exec;
        } else {
            return dcraw_emu_exec;
        }
    }
    
    virtual QStringList arguments(bool bayer_mode, const QString& input, [[maybe_unused]] const QString& output) {
        QStringList args;
        
        if (bayer_mode) {
            // output is hard-coded in unprocessed_raw to be input.tiff, so
            // in this case we are actually working from a copy of input
            args << "-T" << input;
        } else {
            // output goes to stdout
            args << "-w" << "-4" << "-T" << "-q" << "3" << "-o" << "0" 
                << "-s" << "0" << "-H" << "1" << "-Z" << "-" << input;
        }
        
        return args;
    }
    
    virtual QString pre_proc(bool bayer_mode, const QString& input, const QString& output) { 
        QString effective_ifname;
        if (bayer_mode) {
            effective_ifname = QFileInfo(output).absolutePath() + "/" + QFileInfo(input).fileName();
            // make a copy of the input file in the output (temp) dir
            if (QFile::copy(input, effective_ifname)) {
                copied_input = effective_ifname;
                logger.info("copied [%s] to [%s] before raw development\n",
                    input.toLocal8Bit().constData(),
                    effective_ifname.toLocal8Bit().constData()
                );
            }
        } else {
            effective_ifname = input;
        }
        return effective_ifname;
    }
    
    virtual void post_proc(bool bayer_mode, [[maybe_unused]] const QString& input, const QString& output) {
        if (copied_input.length() > 0 && QFile::exists(copied_input)) {
            QFile::remove(copied_input);
        }
        if (bayer_mode) {
            QString raw_ofname = copied_input + ".tiff";
            logger.info("after raw development, move [%s] to [%s] \n",
                raw_ofname.toLocal8Bit().constData(),
                output.toLocal8Bit().constData()
            );
            if (QFile::exists(raw_ofname)) {
                QFile::remove(output); // ensure we can overwrite this file if it exists
                QFile(raw_ofname).rename(output);
            } else {
                QFile::remove(output); // signal that there was no output here either
                logger.debug("Raw developer [%s] did not produce the expected output file [%s]\n",
                    unprocessed_raw_exec.toLocal8Bit().constData(),
                    raw_ofname.toLocal8Bit().constData()
                );
            }
        }
    }
  
    virtual bool accepts_suffix(const QString& suffix) {
        bool found = false;
        for (size_t i=0; !found && i < suffixes.size(); i++) {
            found = suffix.compare(suffixes[i], Qt::CaseInsensitive) == 0;
        }
        return found;
    };
  
    QString dcraw_emu_exec;
    QString unprocessed_raw_exec;
    QString copied_input;
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
      "CR2",  // Canon
      "CR3"  // newer Canon
    };
};

#endif
