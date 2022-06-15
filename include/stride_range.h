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

#ifndef STRIDE_RANGE_H
#define STRIDE_RANGE_H

#include "threadpool.h"

class Stride_range {
  public:
    Stride_range(size_t start, size_t end, size_t stride)
    : first(start), last(end), stride(stride)
    {}
    
    size_t begin(void) const {
        return first;
    }
    
    size_t end(void) const {
        return last + 2*stride; // excluded end value
    }
    
    size_t& increment(size_t& v) const {
        if (v == last) {
            v = last + 2*stride;
            return v;
        } else {
            if (v + stride > last) {
                v = last;
                return v;
            }
        }
        v += stride;
        return v;
    }
    
    template<class T>
    static void parallel_for(T& ftor, ThreadPool& tp, size_t ceiling) {
        size_t stride = min(ceiling, tp.size());
        size_t rpt = ceiling / stride; // rows per thread
        size_t remainder = ceiling - rpt*stride;
        
        vector< std::future<void> > futures;
        for (size_t b=0; b < stride; b++) {
            size_t lower = b;
            size_t upper = b < remainder ? lower + rpt*stride : lower + (rpt-1)*stride;
            Stride_range sr(lower, upper, stride);
            
            futures.emplace_back( 
                tp.enqueue( [sr,&ftor] {
                    ftor(sr);
                })
            );
            
        }
        for (size_t i=0; i < futures.size(); i++) {
            futures[i].wait();
        }
    }
    
    size_t first;
    size_t last;
    size_t stride;
};

#endif
