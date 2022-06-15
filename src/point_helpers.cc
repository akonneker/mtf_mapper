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
#include "include/common_types.h"
#include "include/gradient.h"

Point2d centroid(const Pointlist& p) {
    double mx = 0;
    double my = 0;
    
    for (size_t i=0; i < p.size(); i++) {
        mx += p[i].x;
        my += p[i].y;
    }
    mx /= p.size();
    my /= p.size();
    
    return Point2d(mx, my);
}

Point2d average_dir(const Gradient& g, int x, int y) {
    double mx = 0;
    double my = 0;
    
    double wsum = 0;
    
    int offsets[9][2] = { {0,0}, {-1,0}, {1,0}, {0,-1}, {0,1},
                                 {-1,1}, {1,1}, {1,-1}, {-1,-1} };
    
    for (int k=0; k < 9; k++) {
    
        int lx = x + offsets[k][0];
        int ly = y + offsets[k][1];
        
        
        if (lx >= 0 && lx < g.width() &&
            ly >= 0 && ly < g.height()) {
    
            mx += g.grad_x(lx,ly) * g.grad_magnitude(lx,ly);
            my += g.grad_y(lx,ly) * g.grad_magnitude(lx,ly);
            wsum += g.grad_magnitude(lx,ly);
        }
        
    }
    
    mx /= wsum;
    my /= wsum;
    return Point2d(-my, mx);
}

Point2d normalize(const Point2d& p) {
    Point2d q;
    double norm = sqrt(p.ddot(p));
    q.x = p.x / norm;
    q.y = p.y / norm;
    return q; 
}
