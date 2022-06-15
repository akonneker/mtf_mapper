/*
Copyright 2021 Frans van den Bergh. All rights reserved.

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

#include "histo_widget.h"
#include <QPainter>

#include <cmath>

Histo_widget::Histo_widget(QWidget* parent) 
 : QWidget(parent) {
 
    setMaximumWidth(256);
    setMinimumWidth(256);
    setMinimumHeight(64);
    warning_icon.addFile(":/Icons/Warning");
}

void Histo_widget::paintEvent([[maybe_unused]] QPaintEvent* event) {
    
    QPainter painter(this);
    painter.setCompositionMode(QPainter::CompositionMode_SourceOver);
    painter.fillRect(0, 0, width(), height(), 
        histo_overlap ? QColor::fromRgb(217, 123, 40) : QColor::fromRgb(128+80, 128+80, 128+80)
    );
    
    int data_d_max = 0;
    for (auto& v : data_d) {
        data_d_max = std::max(data_d_max, (int)log(v+1));
    }
    
    int data_l_max = 0;
    for (auto& v : data_l) {
        data_l_max = std::max(data_l_max, (int)log(v+1));
    }
    
    if (data_d.size() > 0 && data_d_max > 0) {
        constexpr int bg_mid = 128;
        painter.setPen(QColor::fromRgb(0, bg_mid - 35, bg_mid - 35));
        for (size_t i=0; i < data_d.size(); i++) {
            double h_d = log(data_d[i]+1) * height() / double(data_d_max);
            painter.drawLine(i, height() - h_d, i, height());
        }
        
        painter.setPen(QColor::fromRgb(0, bg_mid + 35, bg_mid + 35));
        for (size_t i=0; i < data_l.size(); i++) {
            double h_l = log(data_l[i]+1) * height() / double(data_l_max);
            painter.drawLine(i, height() - h_l, i, height());
        }
        
        painter.setPen(QColor::fromRgb(255, 80, 80));
        for (size_t i=0; i < data_d.size(); i++) {
        
            double h_d = log(data_d[i]+1) * height() / double(data_d_max);
            double h_l = log(data_l[i]+1) * height() / double(data_l_max);
            
            if (h_d > 0 && h_l > 0) {
                int h_e = std::min(h_d, h_l);
                painter.drawLine(i, height() - h_e, i, height());
            }
        }
    }
    
    if (data_d.size() > 0 && data_l.size() > 0) {
        if (data_d.front() > 0 || data_l.front() > 0) {
            warning_icon.paint(&painter, 0, 0, 24, 24);
        }
        
        if (data_d.back() > 0 || data_l.back() > 0) {
            warning_icon.paint(&painter, width() - 24, 0, 24, 24);
        }
    }
}

void Histo_widget::set_histogram(histo_t dark, histo_t light) {
    
    // smooth out the histograms a bit for a nicer display
    data_d = histo_t(dark.size(), 0);
    data_l = histo_t(light.size(), 0);
    
    // save the original extrema so we can replace them later
    // these extrema may have a non-linear encoding indicating
    // black crush or saturation
    vector<int> extrema(4, 0);
    if (dark.size() >= 256 && light.size() >= 256) {
        extrema[0] = dark[0];
        extrema[1] = light[0];
        extrema[2] = dark[255];
        extrema[3] = light[255];
    }
    
    // one pass of morphological max
    constexpr int hw = 1;
    for (int i=hw; i < (int)dark.size() - hw; i++) {
        for (int d=-hw; d <= hw; d++) {
            if ((i+d) >= hw && (i+d) < (int)dark.size() - hw) {
                data_d[i] = std::max(data_d[i], dark[i+d]);
                data_l[i] = std::max(data_l[i], light[i+d]);
            }
        }
    }
    // and one pass of box filtering
    dark = data_d;
    light = data_l;
    std::fill(data_d.begin(), data_d.end(), 0);
    std::fill(data_l.begin(), data_l.end(), 0);
    for (int i=hw; i < (int)dark.size() - hw; i++) {
        for (int d=-hw; d <= hw; d++) {
            int j = std::max(hw, std::min(i+d, (int)dark.size() - 1 - hw));
            data_d[i] += dark[j];
            data_l[i] += light[j];
        }
    }
    
    // replace the extrema
    if (dark.size() >= 256 && light.size() >= 256) {
        data_d[0] = extrema[0];
        data_l[0] = extrema[1];
        data_d[255] = extrema[2];
        data_l[255] = extrema[3];
    }
    
    histo_overlap = false;
    if (data_d.size() > 0 && data_l.size() > 0) {
        size_t last_dark = 0;
        for (size_t i=0; i < data_d.size(); i++) {
            if (data_d[i] > 0) {
                last_dark = i;
            }
        }
        size_t first_light = 0;
        for (int i=data_l.size() - 1; i >= 0; i--) {
            if (data_l[i] > 0) {
                first_light = i;
            }
        }
        histo_overlap = first_light <= last_dark;
    }
    
    update();
}


