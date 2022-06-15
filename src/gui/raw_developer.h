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
#ifndef RAW_DEVELOPER_H
#define RAW_DEVELOPER_H

#include "include/logger.h"

#include <memory>

#include <QString>
#include <QStringList>
#include <QProcess>
#include <QObject>
#include <QSettings>

class Raw_developer {
  public:
    Raw_developer(void) {}
    virtual ~Raw_developer(void) {}
    
    int process(const QString& input, const QString& output, bool bayer_mode) {
          QProcess dcp;
          
          dcp.setProgram(program_name(bayer_mode));
          dcp.setStandardOutputFile(output);
          
          QString effective_input = pre_proc(bayer_mode, input, output);
          
          dcp.setArguments(arguments(bayer_mode, effective_input, output));
          
          logger.debug("arguments to raw devloper [%s]:\n", program_name(bayer_mode).toLocal8Bit().constData());
          for (int kk = 0; kk < dcp.arguments().size(); kk++) {
              logger.debug("[%d]=%s\n", kk, dcp.arguments().at(kk).toLocal8Bit().constData());
          }
          dcp.start();
          dcp.waitForFinished(-1);
          int dc_rval = dcp.exitStatus() == QProcess::NormalExit && dcp.exitCode() == 0;
          if (!dc_rval) {
              logger.error("Error. raw developer call failed on input image %s [exit status=%d, exitcode=%d]\n", 
                  effective_input.toLocal8Bit().constData(), dcp.exitStatus(), dcp.exitCode()
              );
          } else {
              post_proc(bayer_mode, input, output);
          }
          return dc_rval;
    }
    
    bool accepts(const QString& suffix) {
        return accepts_suffix(suffix);
    }
    
  private:
    virtual QString program_name(bool bayer_mode) = 0;
    virtual QStringList arguments(bool bayer_mode, const QString& input, const QString& output) = 0;
    virtual bool accepts_suffix(const QString& suffix) = 0;
    
    // to support making a copy of the input file first
    virtual QString pre_proc([[maybe_unused]] bool bayer_mode, [[maybe_unused]] const QString& input, [[maybe_unused]] const QString& output) { 
        return input;
    }
    
    // optional clean up
    virtual void post_proc([[maybe_unused]] bool bayer_mode, [[maybe_unused]] const QString& input, [[maybe_unused]] const QString& output) {}
};

#include "settings_helpers_tab.h"
#include "raw_developer_dcraw.h"
#include "raw_developer_libraw.h"

class Raw_developer_factory {
  public:
    static std::shared_ptr<Raw_developer> build(const Settings_helpers_tab& settings) {
        Raw_developer* res = nullptr;
        switch (settings.get_raw_developer()) {
        case Settings_helpers_tab::raw_developer_libraw:
            res = new Raw_developer_libraw(settings.get_dcraw_emu_binary(), settings.get_unprocessed_raw_binary());
            break;
        case Settings_helpers_tab::raw_developer_dcraw:
            res = new Raw_developer_dcraw(settings.get_dcraw_binary());
            break;
        default:
            logger.error("Unknown raw developer type %d in Settings_helpers_tab, found in Raw developer factory\n", 
                settings.get_raw_developer()
            );
        }
        return std::shared_ptr<Raw_developer>(res);
    }
};

#endif
