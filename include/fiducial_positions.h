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
#ifndef FIDUCIAL_POSITIONS_H
#define FIDUCIAL_POSITIONS_H

#include "common_types.h"
#include <array>

class Fiducial {
  public:
    Fiducial(double img_x=0, double img_y=0, double real_x=0, double real_y=0, int code=0) 
    : icoords(img_x, img_y), rcoords(real_x, real_y), code(code) {}

    Point2d icoords;
    Point2d rcoords;
    int code;
};

const int n_fiducials = 22;

const Fiducial main_fiducials[n_fiducials] = {
    {0, 0, 30.000000, 0.000000, 0},
    {0, 0, -30.000000, 0.000000, 0},
    {0, 0, 52.130517, 36.502181, 1},
    {0, 0, 36.502181, -52.130517, 2},
    {0, 0, -52.130517, -36.502181, 3},
    {0, 0, -36.502181, 52.130517, 4}, 
    {0, 0, 81.915204, -57.357644, 5},
    {0, 0, -58.791585, -83.963085, 6},
    {0, 0, -86.010965, 60.225526, 7}, 
    {0, 0, 61.659467, 88.058845, 8},
    {0, 0, 112.763114, 41.042417, 9},
    {0, 0, 41.897468, -115.112346, 10},
    {0, 0, -117.461578, -42.752518, 11},
    {0, 0, -43.607568, 119.810809, 12}, 
    {0, 0, 44.462619, 122.160041, 13},
    {0, 0, 124.509272, -45.317669, 14},
    {0, 0, -46.172719, -126.858504, 15},
    {0, 0, -129.207735, 47.027770, 16}, 
    {0, 0, 98.994949, 98.994949, 17},
    {0, 0, 100.762716, -100.762716, 18},
    {0, 0, -102.530483, -102.530483, 19},
    {0, 0, -104.298250, 104.298250, 20}
};

const vector< vector<int> > fiducial_code_mapping = {
    // for focus chart
    {0, 0,  5,  6,  7,  8, 1, 2, 3, 4, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20}, // A0
    {0, 0,  9, 10, 11, 12, 5, 6, 7, 8, 1, 2, 3, 4, 13, 14, 15, 16, 17, 18, 19, 20},    // A1
    {0, 0, 13, 14, 15, 16, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 17, 18, 19, 20},     // A2
    {0, 0,  1,  2,  3,  4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20}, // A3
    {0, 0, 17, 18, 19, 20, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1, 2, 3, 4},     // A4
    
    // for lensgrid chart, which has slightly different fiducial positions
    {0, 0, 2,  6, 10, 14, 18, 1, 7, 8, 9, 3, 11, 12, 13, 4, 15, 16, 17, 5, 19, 20},    // A0
    {0, 0, 3,  7, 11, 15, 19, 6, 2, 8, 9, 10, 1, 12, 13, 14, 4, 16, 17, 18, 5, 20},    // A1
    {0, 0, 4,  8, 12, 16, 20, 6, 7, 2, 9, 10, 11, 3, 13, 14, 15, 1, 17, 18, 19, 5},    // A2
    {0, 0, 1,  5,  9, 13, 17, 2, 7, 8, 3, 10, 11, 12, 4, 14, 15, 16, 6, 18, 19, 20},   // A3
    {0, 0, 6, 11, 16, 17,  5, 1, 7, 8, 9, 10, 2, 12, 13, 14, 15, 3, 4, 18, 19, 20},    // A4; this code might be weaker than the rest in terms of Hamming distance
};

// scaling factors of fiducials relative to A3 size (main_fiducials are designed for A3)
const std::array<double, 10> fiducial_scale_factor = {
    // scale factors for 45-degree focus chart
    2.83095238095238095238,  // A0
    2.00238095238095238095,  // A1 
    1.41428571428571428571,  // A2 
    1,                       // A3
    0.707142857142857,       // A4
    
    // repeated for lensgrid chart
    2.83095238095238095238,  // A0
    2.00238095238095238095,  // A1 
    1.41428571428571428571,  // A2 
    1,                       // A3
    0.707142857142857        // A4
};

const std::vector< pair<double, double> > fiducial_position_scale_factor = {
  // position scales in main_fiducials are unchanged for focus chart type
  {1.0, 1.0},
  {1.0, 1.0},
  {1.0, 1.0},
  {1.0, 1.0},
  {1.0, 1.0},
  
  // position scales in main_fiducials are stretched to improve sensitivty to rotation
  // in lensgrid chart type
  {1.05, 1.58},
  {1.05, 1.58},
  {1.05, 1.58},
  {1.05, 1.58},
  {1.05, 1.58}
};

// this enum provides an index into fiducial_code_mapping
typedef enum {
  // focus chart type of various sizes
  A0 = 0,
  A1 = 1,
  A2 = 2,
  A3 = 3,
  A4 = 4,
  
  // lensgrid chart type of various sizes
  A0L = 5,
  A1L = 6,
  A2L = 7,
  A3L = 8,
  A4L = 9
} fiducial_mapping_index;

#endif
