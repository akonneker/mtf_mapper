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
#ifndef GL_IMAGE_PANEL_DOTS_H
#define GL_IMAGE_PANEL_DOTS_H

#include "gl_image_panel.h"
#include <memory>

class mtfmapper_app;

class GL_image_panel_dots : public GL_image_panel {
  public:
    explicit GL_image_panel_dots(QWidget* parent = 0, mtfmapper_app* app = 0);
    
    virtual void click_marker(QPoint pos, bool add=false) override;
    virtual void clear_overlay(void) override { dot_list.clear(); }
    virtual void initialize_overlay(void) override;
    virtual void paint_overlay(void) override;
    
    void mouseMoveEvent(QMouseEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;
    void mouseReleaseEvent(QMouseEvent* event) override;
    
  private:
    void make_dots();
    void draw_dot(double x, double y, double r, double g, double b);
    
    mtfmapper_app* app = nullptr;
    
    std::shared_ptr<QOpenGLShaderProgram> dots_program;
    QOpenGLBuffer dots_vbo;
    std::map<std::string, vector<Sfr_marker> > dot_list;
};

#endif

