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
#ifndef COMPONENT_LABELLING_H
#define COMPONENT_LABELLING_H

#include "common_types.h"
#include "include/logger.h"
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <map>

// A class for extracting the boundaries of a binary image.
//
// It expects black objects on white backgrounds.
//
// Implements the method in (F. Chang, C.-J. Chen, C.-J. Jen, A linear-time
// component labeling algorithm using contour tracing technique, Computer
// Vision and Image Understanding, 93:206-220, 2004)

//==============================================================================
class Component_labeller {
public:
    Component_labeller(void);
    Component_labeller(const cv::Mat& in_img,
        int min_boundary_length = 10, bool snapshot = false, 
        int max_boundary_length = 5000);

    ~Component_labeller(void);
    
    void release(void) {
        _pix_data.clear();
        _pix_data.shrink_to_fit();
        _labels.clear();
        _labels.shrink_to_fit();
        configured = false;
    }

    void configure(const cv::Mat& in_img,
        int min_boundary_length = 10, 
        int max_boundary_length = 5000,
        bool snapshot = false);


    const Boundarylist& get_boundaries(void) const {
        assert(configured);
        return _boundaries;
    }
    
    int largest_hole(int label) const {
        auto it = _holes.find(label);
        if (it != _holes.end()) {
            return it->second;
        }
        return 0;
    }

    inline int operator[](int index) const {
        return _labels[index];
    }

    inline int operator()(int x, int y) const {
        if (!(x >= 0 && x < _width && y >= 0 && y < _height)) {
            logger.debug("trying to access %d, %d (size is %d,%d)\n", x, y, _width, _height);
            return -1;
        }
        return _labels[y * _width + x];
    }

    inline int get_width(void) const {
        return _width;
    }

    inline int get_height(void) const {
        return _height;
    }
    
    // set borders to the background color, since objects that
    // touch the image borders cause havoc with the main algorithm
    static void zap_borders(cv::Mat& masked_img, int fill=255) {
        for (int y=0; y < masked_img.rows; y++) {
            masked_img.at<uchar>(y,0) = fill;
            masked_img.at<uchar>(y,1) = fill;
            masked_img.at<uchar>(y,masked_img.cols-1) = fill;
            masked_img.at<uchar>(y,masked_img.cols-2) = fill;
        }
        for (int x=0; x < masked_img.cols; x++) {
            masked_img.at<uchar>(0,x) = fill;
            masked_img.at<uchar>(1,x) = fill;
            masked_img.at<uchar>(masked_img.rows-1,x) = fill;
            masked_img.at<uchar>(masked_img.rows-2,x) = fill;
        }
    }
    
    void inflate_boundaries(double radius);

private:
    typedef enum {
        INTERNAL = 0,
        EXTERNAL = 1,
        EXTERNAL_FIRST = 2,
        INTERNAL_FIRST = 3
    } mode_type;

    void _find_components(void);

    void _contour_tracing(int x, int y, int label, mode_type mode);

    void _tracer(int x, int y, int& nx, int& ny,
        int& from, mode_type mode, bool mark_white = true);

    void _draw_snapshot(void);

    int _width;
    int _height;

    vector<uint8_t> _pix_data;
    uint8_t* _pix;
    vector<int32_t> _labels;
    Boundarylist _boundaries;
    std::map<int, int> _holes;

    int _min_boundary_length;
    int _max_boundary_length;

    int C;

    bool configured;
};

#endif // COMPONENT_LABELLING_H


