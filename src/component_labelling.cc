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
#include "include/component_labelling.h"
#include <opencv2/imgcodecs/imgcodecs.hpp>

//------------------------------------------------------------------------------
Component_labeller::Component_labeller(void) 
: _width(0), _height(0),
  _min_boundary_length(0), 
  _max_boundary_length(0),
  configured(false) {
    // nothing to do
}

//------------------------------------------------------------------------------
Component_labeller::Component_labeller(const cv::Mat& in_img,
    int min_boundary_length, bool snapshot, int max_boundary_length)
: _width(in_img.cols), _height(in_img.rows),
  _min_boundary_length(min_boundary_length), 
  _max_boundary_length(max_boundary_length), configured(false)
{
    configure(in_img, min_boundary_length, max_boundary_length, snapshot);
}

//------------------------------------------------------------------------------
void Component_labeller::configure(const cv::Mat& in_img,
        int min_boundary_length, int max_boundary_length, bool snapshot) {
    
    if (configured) {
        _boundaries.clear();
    }

    _width  = in_img.cols;
    _height = in_img.rows;
    _min_boundary_length = min_boundary_length;
    _max_boundary_length = max_boundary_length;
    configured = true;

    _labels.resize(in_img.rows*in_img.cols, 0);

    // make a copy of the pixel data, but add an extra
    // white line to the top of the image
    _pix_data.resize(in_img.rows*in_img.cols + _width);
    memset(_pix_data.data(), 255, _width);
    _pix = _pix_data.data() + _width;
    memcpy(_pix, in_img.data, in_img.rows*in_img.cols);
    memset(_pix, 255, _width);

    // set first and last columns to white
    for (int y=0; y < _height; y++) {
        _pix[y*_width] = 255;
        _pix[y*_width + _width-1] = 255;
    }

    C = 1; // current label

    _find_components();

    if (snapshot) {
        _draw_snapshot();
    }

}


//------------------------------------------------------------------------------
Component_labeller::~Component_labeller(void) {
    if (configured) {
        release();
    }
    configured = false;
}


//------------------------------------------------------------------------------
void Component_labeller::_find_components(void) {

    for (int y=0; y < _height-1; y++) {
        for (int x=0; x < _width; x++) {

            int pos = y * _width + x;
            // find next black pixel
            if (_pix[pos] != 0) {
                continue;
            }

            int done = false;
            // pixel at [pos] is black
            if (_labels[pos] == 0 && _pix[pos - _width] != 0) {
                // not labelled, and pixel above it is white : step 1
                _labels[pos] = C;

                // trace external contour, label as 'C'
                _contour_tracing(x, y, C, EXTERNAL);

                C++;
                done = true;
            } else {
                if (_pix[pos + _width] != 0 && _labels[pos + _width] == 0) {
                    // pixel below current is unmarked step 2:
                    if (_labels[pos] == 0) {
                        // internal contour
                        _labels[pos] = _labels[pos-1];

                        // trace internal contour, label as '_labels[pos]'
                        _contour_tracing(x, y, _labels[pos], INTERNAL);

                    } else {
                        // trace internal contour, label as '_labels[pos]'
                        _contour_tracing(x, y, _labels[pos], INTERNAL);

                    }
                    done = true;
                }
            }
            if (!done && _labels[pos] == 0) {
                // step 3
                _labels[pos] = _labels[pos-1];
            }
        }
    }
}

//------------------------------------------------------------------------------
void Component_labeller::_contour_tracing(int x, int y, int label, mode_type mode) {
    int nx;
    int ny;
    int from = 0;

    // call tracer
    _tracer(x, y, nx, ny, from, mode == EXTERNAL ? EXTERNAL_FIRST : INTERNAL_FIRST);

    // tracer returns nx < 0 if initial point was isolated
    if (nx < 0) {
        // isolated point. already labelled
        return;
    }

    // external trace, first pixel of new component.
    // mark the white pixel above, if necessary
    int pos = (y-1) * _width + x;
    if (mode == EXTERNAL && _pix[pos] != 0) {
        _labels[pos] = -1;
    }

    int initial_nx = nx;
    int initial_ny = ny;
    int initial_x = x;
    int initial_y = y;

    bool done = false;

    Pointlist boundary;
    boundary.push_back(Point2d(nx,ny));

    while (!done) {
        x = nx;
        y = ny;
        // call tracer
        _tracer(x, y, nx, ny, from, mode, true);

        boundary.push_back(Point2d(nx,ny));

        if (nx < 0) {
            break;
        }

        if (_labels[x + y*_width] == 0) {
            _labels[x + y*_width] = label;
        }

        // must check if next point is (initial_nx, initial_ny)
        int dummy_from = from;
        int dummy_nx;
        int dummy_ny;
        _tracer(nx, ny, dummy_nx, dummy_ny, dummy_from, mode, false);
        if (dummy_nx == initial_nx && dummy_ny == initial_ny &&
            nx == initial_x && ny == initial_y) {
            done = true;
        }
    }
    if (mode == EXTERNAL && 
        boundary.size() >= (size_t)_min_boundary_length &&
        boundary.size() <= (size_t)_max_boundary_length) {
        _boundaries.insert(make_pair(label, boundary));
    }
    
    if (mode == INTERNAL && boundary.size() >= 2) {
        auto it = _holes.find(label);
        if (it != _holes.end()) {
            _holes[label] = std::max(_holes[label], int(boundary.size()));
        } else {
            _holes.insert(make_pair(label, int(boundary.size())));
        }
    }
    

}

//------------------------------------------------------------------------------
void Component_labeller::_tracer(int x, int y, int& nx, int& ny,
    int& from, mode_type mode, bool mark_white) {

    int neighbours[8][2] = {
        {1, 0}, {1, 1}, {0, 1}, {-1,1},
        {-1, 0}, {-1, -1}, {0, -1}, {1, -1}
    };

    int dir = 0;
    switch (mode) {
    case EXTERNAL_FIRST:
        dir = 7;
        break;
    case INTERNAL_FIRST:
        dir = 3;
        break;
    default:
        dir = (from + 2) % 8;
    }

    for (int n=0; n < 8; n++) {
        nx = (x+neighbours[dir][0]);
        ny = (y+neighbours[dir][1]);

        // only consider pixels inside the image
        if (nx >= _width || nx < 0 || ny >= _height || ny < 0) {
            continue;
        }

        int pos = nx + _width*ny;

        if (_pix[pos] != 0) {
            // white pixel, mark it ... not sure if this is sufficient?
            if (mark_white) {
                _labels[pos] = -1;
            }
        } else {
            // found the next black pixel
            from = (dir + 4) % 8;
            return;
        }

        dir = (dir + 1) % 8;
    }

    // isolated point, return invalid next coordinates
    nx = ny = -1;
}


//------------------------------------------------------------------------------
void Component_labeller::_draw_snapshot(void) {

    cv::Mat img(_height, _width, CV_8UC1);
    for (int i=0; i < img.rows*img.cols; i++) {
        img.data[i] = _labels[i] < 0 ? 255 : _labels[i];
    }

    imwrite(std::string("labelled.png"), img);
}


void Component_labeller::inflate_boundaries(double radius) {

    constexpr size_t min_boundary_length = 4 * 10; // a 10x10 square, maybe?

    for (auto& ble: _boundaries) {
        auto& b = ble.second;
        
        if (b.size() > min_boundary_length) {
            Point2d centroid(0, 0);
            
            for (auto& p: b) {
                centroid += p;
            }
            centroid *= 1.0/b.size();
            
            // move each point outward by radius along the vector
            // connecting the point to the centroid
            Point2d dir;
            for (auto& p: b) {
                dir = p - centroid;
                double l = norm(dir);
                if (l > 0) { // do not try to inflate the centroid itself!
                    dir *= (l+radius)/l;
                    p = dir + centroid;
                }
            }
            
        }
        
    }
}
