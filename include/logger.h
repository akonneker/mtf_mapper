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
#ifndef LOGGER_H
#define LOGGER_H

#include <stdio.h>
#include <string>
using std::string;

class Logger {
  public: 
      typedef enum {LOGGER_NONE=0, LOGGER_INFO=1, LOGGER_DEBUG=2, LOGGER_ERROR=4} log_level_t;
      Logger(void) : destination(stdout) {}

      void redirect(const string& fname, bool append=false) {
          destination = fopen(fname.c_str(), append ? "a" : "w");
          if (!destination) {
              fprintf(stderr, "Could not create %s for writing\n", fname.c_str());
              destination = stdout;
          }
      }

      ~Logger(void) {
          if (destination != stdout) {
              fclose(destination);
          }
      }
      
      void flush(void) {
          fflush(destination);
      }

      void enable_level(log_level_t ll) {
          log_level |= ll;
      }

      void disable_level(log_level_t ll) {
          log_level &= ~ll;
      }

      template <class... T> void info(T... t) {
          if (log_level & LOGGER_INFO) {
              fprintf(destination, t...);
          }
      }

      template <class... T> void debug(T... t) {
          if (log_level & LOGGER_DEBUG) {
              fprintf(destination, t...);
          }
      }

      template <class... T> void error(T... t) {
          if (log_level & LOGGER_ERROR) {
              fprintf(destination, t...);
              fflush(destination);
          }
      }
  private:
      FILE* destination;
      int log_level = LOGGER_INFO | LOGGER_ERROR;
};

extern Logger logger;

#endif
