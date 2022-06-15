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
#ifndef GL_IMAGE_VIEWER_H
#define GL_IMAGE_VIEWER_H

#include <QtWidgets/QWidget>
#include <QtWidgets/QAbstractScrollArea>
#include "gl_image_panel.h"

class GL_image_viewer : public QAbstractScrollArea {
    Q_OBJECT
    
  public:
    explicit GL_image_viewer(QWidget* parent);
    
    bool viewportEvent(QEvent* e) override;
    void scrollContentsBy(int dx, int dy) override;
    void wheelEvent(QWheelEvent* e) override;
    
    void mouseMoveEvent(QMouseEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;
    void mouseReleaseEvent(QMouseEvent* event) override;
    void keyPressEvent(QKeyEvent* event) override;
    
    void set_GL_widget(GL_image_panel* w);
    bool load_image(const QString& fname);
    bool load_image(QImage* qimg);
    void set_clickable(bool b);
    void set_resize_on_load(bool b);

    QSize sizeHint(void) const override;
    
  private:
    void zoom_action(double direction, int zx, int zy);
    
    GL_image_panel* widget;
    
    bool panning = false;
    QPoint pan;
    QPoint click;
    
    bool zooming = false;
    QPoint zoom_pos;
    QPoint zoom_pos_temp;
    
    QPoint last_mouse_pos;
    
    bool must_update_bars = true;
    bool is_clickable = false;
    QSize hinted_size = QSize(0, 0);
    bool resize_on_load = false;
    
  public slots:
    void clear_overlay();
};

#endif
