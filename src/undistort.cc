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


#include "include/undistort.h"

// find closest point in radmap lookup table, then apply
// quadratic interpolation to refine the value
cv::Point2d Undistort::transform_point(double c, double r) {
    double px = (c + offset.x - centre.x);
    double py = (r + offset.y - centre.y);
    
    double rad = 2*sqrt(px*px + py*py); // convert to half-pixel spacing to match how radmap was built
    int rad_f = (int)rad;
    
    if (rad_f > (int)radmap.size() - 3) {
        if (radmap.size() == 0) { // this should only happen if we call transform_point before unmap, which should not happen
            logger.error("%s\n", "Warning: radmap not initialized in transform_point. Calling build_radmap, but fix the code!");
            build_radmap();
        } else {
            logger.error("Warning: radmap range exceeded in transform_point. Clamping (rad-%lf, rad_f=%d)\n", rad, rad_f);
            rad_f = (int)radmap.size() - 3;
        }
    }
    
    double w = rad - rad_f;
    // quadratic interpolation through half-pixel spaced points
    double fc = radmap[rad_f];
    double fb = -1.5*radmap[rad_f] + 2*radmap[rad_f+1] - 0.5*radmap[rad_f+2];
    double fa =  0.5*radmap[rad_f] - radmap[rad_f+1] + 0.5*radmap[rad_f+2];
    double rad_d = 2*(fa*w*w + fb*w + fc);
    
    // avoid divide-by-zero
    if (rad < 1e-8) {
        return centre;
    }
    return Point2d(px, py) * (rad_d/rad) + centre - offset;
}

void Undistort::build_radmap(void) {
    radmap.clear();
    double maxrad = 1.1*sqrt(centre.x*centre.x + centre.y*centre.y); // TODO: add better support for non-centered CoD
    for (int rad=0; rad <= 2*(maxrad + 2); rad++) { // half-pixel spacing
        Point2d tp = slow_transform_point(centre.x + 0.5*rad, centre.y);
        double rd = norm(tp - Point2d(centre.x, centre.y));
        radmap.push_back(rd);
    }
}

cv::Mat Undistort::unmap_base(const cv::Mat& in_src, cv::Mat& rawimg, int pad_left, int pad_top) {
    // to preserve Bayer CFA alignment
    pad_left += pad_left % 2;
    pad_top += pad_top % 2;
    
    cv::Mat src;
    copyMakeBorder(in_src, src, pad_top, pad_top, pad_left, pad_left, cv::BORDER_CONSTANT, cv::Scalar::all(0));
    centre.x = src.cols / 2;
    centre.y = src.rows / 2;
    
    logger.debug("Padding with %d left/right pixels, %d top/bottom pixels\n", pad_left, pad_top);
    
    if (pad_left > 0 || pad_top > 0) {
        logger.debug("%s\n", "Padding distorted/Bayer image.");
        cv::Mat rcopy = rawimg.clone();
        copyMakeBorder(rcopy, rawimg, pad_top, pad_top, pad_left, pad_left, cv::BORDER_CONSTANT, cv::Scalar::all(0));
        last_padding = cv::Point2i(pad_left, pad_top);
    }
    
    build_radmap();

    cv::Mat map_x(src.rows, src.cols, CV_32FC1);
    cv::Mat map_y(src.rows, src.cols, CV_32FC1);
    
    cv::Mat blurred;
    cv::blur(src, blurred, cv::Size(3, 3));
    
    Point2d prev;
    for (int r=0; r < src.rows; r++) {
        prev = transform_point(-1, r);
        for (int c=0; c < src.cols; c++) {
            Point2d tp = transform_point(c, r);
            map_x.at<float>(r, c) = tp.x;
            map_y.at<float>(r, c) = tp.y;
            
            // add some blurring to suppress noise before magnification
            // this avoids the edges becoming so rough that the rectangle
            // test on the undistorted image fails
            double stretch = norm(prev - tp);
            if (stretch > 0.9) {
                blurred.at<uint16_t>(r, c) = src.at<uint16_t>(r, c);
            } else {
                // 0.9 -> src 
                // 0.4 -> blurred 
                double w = stretch < 0.4 ? 0 : 2*(stretch - 0.4);
                blurred.at<uint16_t>(r, c) = w*src.at<uint16_t>(r, c) + (1-w)*blurred.at<uint16_t>(r, c);
            }
            prev = tp;
        }
    }
    cv::Mat timg;
    cv::remap(blurred, timg, map_x, map_y, cv::INTER_LINEAR, cv::BORDER_CONSTANT, cv::Scalar::all(0));
    
    return timg;
}

void Undistort::estimate_padding(const cv::Mat& src, int& pad_left, int& pad_top) {
    // this method scans along the longest edge of the image to see if an edge of reasonable length,
    // say 80 pixels will map to a marginal-length edge (say 20 pixels)
    // this becomes the cut-off on the longest edge, which also 
    // defines the short edge cut-off

    const double src_w = 80;  // pixels
    const double dest_w = 20; // do not try to expand beyond this

    Point2d extreme = inverse_transform_point(1, 1);

    if (!allow_crop) {
        pad_left = extreme.x < 0 ? -ceil(extreme.x) : 0;
        pad_top = extreme.y < 0 ? -ceil(extreme.y) : 0;
        logger.debug("extreme: %lf, %lf. pad_left=%d, pad_top=%d\n", extreme.x, extreme.y, pad_left, pad_top);
        return;
    }

    Point2d sdir(1, 0);
    if (src.rows > src.cols) {
        sdir = Point2d(0, 1);
    }

    const double rad_min = centre.ddot(sdir);
    double rad_max = centre.ddot(sdir) - extreme.ddot(sdir);
    Point2d dir = centre * (1.0 / norm(centre));
    double rad = rad_min;
    // rather sweep along the top and left edges (but still measure radially?)
    for (; rad <= rad_max; rad += 10.0) {
        Point2d p;
        if (sdir.x == 1.0) {
            p = Point2d(centre.x - rad, 0);
        }
        else {
            p = Point2d(0, centre.y - rad);
        }
        Point2d p2 = p - src_w*dir;

        Point2d tp = slow_transform_point(p.x, p.y);
        Point2d tp2 = slow_transform_point(p2.x, p2.y);
        double dist = norm(tp2 - tp);
        if (dist < dest_w) break;
    }
    if (sdir.x == 1.0) {
        pad_left = std::max(-(centre.x - rad - src_w), 0.0);
        Point2d fwd = slow_transform_point(-pad_left, 0);
        pad_top = fwd.y + src_w;
    }
    else {
        pad_top = std::max(-(centre.y - rad - src_w), 0.0);
        Point2d fwd = slow_transform_point(0, -pad_top);
        pad_left = fwd.x + src_w;
    }
}

void Undistort::apply_padding(vector<cv::Mat>& images) {
    if (last_padding.x > 0 || last_padding.y > 0) {
        for (size_t i = 0; i < images.size(); i++) {
            cv::Mat rcopy = images[i].clone();
            copyMakeBorder(rcopy, images[i], last_padding.y, last_padding.y , last_padding.x, last_padding.x, cv::BORDER_CONSTANT, cv::Scalar::all(0));
        }
    }
}
