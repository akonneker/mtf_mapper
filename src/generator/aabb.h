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
#ifndef AABB_GEOM_H
#define AABB_GEOM_H

//==============================================================================
class Aa_bb  { // Axis-aligned bounding box
  public:
      Aa_bb(double min_x=0, double min_y=0, double max_x=0, double max_y=0) 
      : min_x(min_x), min_y(min_y), max_x(max_x), max_y(max_y),
        area((max_x - min_x)*(max_y - min_y)) {
      }
      
      inline bool bounds_overlap(const Aa_bb& b, double xoffset=0, double yoffset=0) const {
          
          if ((max_x+xoffset) < b.min_x || // *this is left of b
              (min_x+xoffset) > b.max_x) {  // *this is right of b
              
              return false;
          }
          
          if ((max_y+yoffset) < b.min_y || // *this is above b
              (min_y+yoffset) > b.max_y) {  // *this is below b
              
              return false;
          }
          
          return true;
      }
      
      inline bool is_inside(double x, double y, double xoffset=0, double yoffset=0) const {
          // xoffset and yoffset probably do not apply here ...
          if (x > (max_x+xoffset) ||  
              x < (min_x+xoffset)) {  
              
              return false;
          }
          
          if (y > (max_y+yoffset) || 
              y < (min_y+yoffset)) {  
              
              return false;
          }
          
          return true;
      }
      

      double min_x;
      double min_y;
      double max_x;
      double max_y;
      
      double area;
};

#endif // AABB_GEOM_H
