/*
Copyright 2020 Frans van den Bergh. All rights reserved.

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
#ifndef GL_IMAGE_PANEL_EDGES_NESTED_H
#define GL_IMAGE_PANEL_EDGES_NESTED_H

#include "gl_image_panel_edges.h"

static double edist(const QPointF& a, const QPointF& b) {
    QPointF d(b-a);
    return sqrt(d.x()*d.x() + d.y()*d.y());
}

GL_image_panel_edges::GL_roi::GL_roi(void) {
    pts[0] = QPointF(-1, -1);
    pts[1] = QPointF(-1, -1);
}

GL_image_panel_edges::GL_roi::GL_roi(QPointF p1, QPointF p2) {
    pts[0] = p1;
    pts[1] = p2;
    n = 2;
}

void GL_image_panel_edges::GL_roi::add_point(QPointF p) {
    if (n < 1) {
        pts[0] = p;
        n++;
    } else {
        pts[1] = p;
    }
}

QPointF& GL_image_panel_edges::GL_roi::get(int i) {
    assert(i <= 1);
    return pts[i];
}

double GL_image_panel_edges::GL_roi::length(void) const {
    if (n >= 1) {
        return edist(pts[0], pts[1]);
    }
    return -1.0;
}

int GL_image_panel_edges::GL_roi::handle_selected(QPointF p, double dist_thresh) {
    for (size_t i=0; i < pts.size(); i++) {
        double d = edist(p, pts[i]);
        
        if (d < dist_thresh) {
            return (int)i;
        }
    }
    return -1;
}

int GL_image_panel_edges::GL_roi::box_selected(QPointF p, double dist_thresh, double width_thresh) {
    int handle_idx = handle_selected(p, dist_thresh);
    
    if (handle_idx < 0) {
        QPointF dir = pts[1] - pts[0];
        double l = sqrt(dir.x()*dir.x() + dir.y()*dir.y());
        dir /= l;
        QPointF norm(-dir.y(), dir.x());
        
        double par = QPointF::dotProduct(p - pts[0], dir);
        double perp = QPointF::dotProduct(p - pts[0], norm);
        
        if (par >= 0 && par < l && fabs(perp) <= width_thresh) {
            return 1;
        }
    }
    
    return 0;
}



GL_image_panel_edges::GL_closebox::GL_closebox(QPointF handle_a, QPointF handle_b,
    [[maybe_unused]] int img_width, [[maybe_unused]] int img_height) {
    
    if (handle_a.x() > handle_b.x()) {
        std::swap(handle_a, handle_b);
    }
    
    dir = handle_b - handle_a;
    double l = sqrt(dir.x()*dir.x() + dir.y()*dir.y());
    dir /= l;
    QPointF norm(-dir.y(), dir.x());
    
    vector<QPointF> candidates = {
        handle_a + l*dir - norm*28,
        handle_a + l*dir + norm*28,
        handle_a - norm*28,
        handle_a + norm*28
    };
    
    bool found = false;
    
    // first go through candidates in standard order,
    // looking for one that falls inside the image
    for (auto& p : candidates) {
        if (p.x() > 0 && p.x() < img_width - 1 &&
            p.y() > 0 && p.y() < img_height - 1) {
          
            found = true;
            pos = p;
            break;
        }          
    }
    
    
    // next, try to choose the one closest to the image centre
    if (!found) {
        const QPointF ic(img_width/2, img_height/2);
        double min_dist = 1e20;
        
        for (const auto& p : candidates) {
            double dist = edist(p, ic);
            if (dist < min_dist) {
                pos = p;
                min_dist = dist;
            }          
        }
    }
}

bool GL_image_panel_edges::GL_closebox::selected(QPointF p, double dist_thresh) {
    QPointF norm(-dir.y(), dir.x());
    
    double par = QPointF::dotProduct(p - pos, dir);
    double perp = QPointF::dotProduct(p - pos, norm);
    
    if (fabs(par) <= dist_thresh && fabs(perp) <= dist_thresh) {
        return true;
    }
    
    return false;
}
    

#endif
