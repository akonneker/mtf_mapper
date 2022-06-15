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
#ifndef TCLAP_CONSTRAINTS_H
#define TCLAP_CONSTRAINTS_H

#include "tclap/Constraint.h"
#include <string>
using std::string;

template <class T>
class tclap_zero_one : public TCLAP::Constraint<T> {
  public:
    virtual ~tclap_zero_one(void) {}
    
    string description(void) const {
        return "Real value restricted to range [0.0, 1.0]";
    }
    
    string shortID(void) const {
        return "range [0, 1]";
    }
    
    bool check(const T& v) const {
        return v >= 0 && v <= 1.0;
    }
    
};

template <class T>
class tclap_ranged : public TCLAP::Constraint<T> {
  public:
    tclap_ranged(const T& lower, const T& upper) : lower(lower), upper(upper) {}
    virtual ~tclap_ranged(void) {}
    
    string description(void) const {
        char buf[1024];
        if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>) {
            sprintf(buf, "A real value restricted to range [%g, %g]", lower, upper);
            return string(buf);
        } else {
            sprintf(buf, "An integer value restricted to range [%d, %d]", lower, upper);
            return string(buf);
        }
    }
    
    string shortID(void) const {
        char buf[1024];
        if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>) {
            sprintf(buf, "range [%g, %g]", lower, upper);
            return string(buf);
        } else {
            sprintf(buf, "range [%d, %d]", lower, upper);
            return string(buf);
        }
    }
    
    bool check(const T& v) const {
        return v >= lower && v <= upper;
    }
    
  private:
    T lower;
    T upper;
};

template <class T>
class tclap_nonneg : public TCLAP::Constraint<T> {
  public:
    virtual ~tclap_nonneg(void) {}
    
    string description(void) const {
        return "A non-negative value, i.e., greater or equal to zero";
    }
    
    string shortID(void) const {
        if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>) {
            return "positive number";
        } else {
            return "positive integer";
        }
    }
    
    bool check(const T& v) const {
        return v >= 0;
    }
    
};

#endif
