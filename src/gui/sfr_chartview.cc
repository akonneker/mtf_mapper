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
#include <QtWidgets> 
#include "sfr_chartview.h"
#include "sfr_dialog.h"

Sfr_chartview::Sfr_chartview(QChart* chart, Sfr_dialog* sfr_dialog) : QChartView(chart), chart(chart), sfr_dialog(sfr_dialog) {
}

void Sfr_chartview::mouseMoveEvent(QMouseEvent* event) {
    QPointF pos = chart->mapToValue(event->pos());
    sfr_dialog->notify_mouse_position(pos.x(), false);
    QChartView::mouseMoveEvent(event);
}

void Sfr_chartview::mousePressEvent(QMouseEvent* event) {
    QPointF pos = chart->mapToValue(event->pos());
    sfr_dialog->notify_mouse_position(pos.x(), (event->buttons() & Qt::LeftButton) ==  Qt::LeftButton);
    QChartView::mousePressEvent(event);
}

void Sfr_chartview::enterEvent(QEvent* event) {
    QChartView::enterEvent(event);
    viewport()->setCursor(Qt::CrossCursor);
}

