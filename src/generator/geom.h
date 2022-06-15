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
#ifndef GENERAL_GEOM_H
#define GENERAL_GEOM_H

#include <cstdio>

#include "aabb.h"

//==============================================================================
class Geometry  {
  public:
      Geometry(double cx=0, double cy=0, double own_area=0, Aa_bb bounds=Aa_bb())
      : cx(cx), cy(cy), own_area(own_area), bounds(bounds) {
      }
      
      Geometry(Aa_bb in_bounds)
      : cx(0.5*(in_bounds.min_x + in_bounds.max_y)), 
        cy(0.5*(in_bounds.min_y + in_bounds.max_y)), 
        own_area(in_bounds.area), bounds(in_bounds) {
      }
      
      virtual ~Geometry(void) {
      }
      
      virtual double intersection_area(const Geometry&, double, double) const {
          printf("\n\nnot defined\n\n");
          return 0;
      }

      virtual bool is_inside(double, double) const {
          printf("\n\nnot defined\n\n");
          return 0;
      }
      
      virtual void print(void) const {
          printf("\n\nnot defined\n");
      }

      double cx;
      double cy;
      double own_area;
      
      Aa_bb bounds;
};

#endif // GENERAL_GEOM_H
