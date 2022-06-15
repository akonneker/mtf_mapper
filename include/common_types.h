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
#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H

#include <vector>
#include <list>
#include <map>
#include <cmath>
#include <algorithm>

using std::vector;
using std::pair;
using std::list;
using std::map;
using std::max;
using std::min;

#include <iostream>
#include <string>

using std::cout;
using std::endl;

// OpenCV headers
#include <opencv2/core/core.hpp>

using cv::Point2d;

typedef vector<Point2d> Pointlist;
typedef map<int, Pointlist> Boundarylist;

using std::make_pair;

#ifdef _WIN32
	#define EXE_SUFFIX ".exe"
    #define _WIN32_WINNT 0x0502
#else
	#define EXE_SUFFIX ""
#endif

#ifndef M_PI
    #define M_PI 3.14159265358979
#endif

#ifndef uint16_t
    typedef unsigned short int uint16_t;
#endif

#ifndef uint32_t
    typedef unsigned int uint32_t;
#endif

#ifndef uint32_t
    typedef int int32_t;
#endif

#if _MSC_VER == 1600 
	typedef unsigned short int uint16_t;
	typedef unsigned int uint32_t;
	typedef int int32_t;

	#define M_PI 3.14159265358979
	
    __inline long int 
	lrint (double flt) 	{	
		int intgr;

		_asm {
			fld flt
			fistp intgr
		}
			
		return intgr;
	} 

#endif

#include "include/scanline.h"

#define SQR(x) ((x)*(x))

// stat functions
#include <sys/stat.h>
#ifndef _WIN32
#define STAT stat
#else
#define STAT _stat
#endif

#ifdef _WIN32
    #ifndef S_ISDIR
        #define S_ISDIR(mode)  (((mode) & S_IFMT) == S_IFDIR)
    #endif
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#endif

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

#endif //COMMON_TYPES_H


