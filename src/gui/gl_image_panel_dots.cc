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
#include "include/logger.h"
#include "gl_image_panel_dots.h"
#include <QOpenGLShaderProgram>
#include <QOpenGLTexture>
#include <QMouseEvent>
#include <opencv2/imgcodecs/imgcodecs.hpp>

#include "mtfmapper_app.h"

#include <cmath>

GL_image_panel_dots::GL_image_panel_dots(QWidget *parent, mtfmapper_app* app)
    : GL_image_panel(parent), app(app) {
      
    logger.debug("GL_image_panel_dots ctor: OpenGL version: %d.%d, samples=%d\n", format().majorVersion(), format().minorVersion(), format().samples());
}


void GL_image_panel_dots::initialize_overlay(void) {
    make_dots();
    

    // shader program for click-dots
    QOpenGLShader *dots_vshader = new QOpenGLShader(QOpenGLShader::Vertex, this);
    const char *dots_vsrc =
        "#version 120\n"
        "attribute highp vec4 vertex;\n"
        "attribute mediump vec4 texCoord;\n"
        "varying mediump vec4 texc;\n"
        "uniform mediump mat4 projectionMatrix;\n"
        "uniform mediump mat4 modelMatrix;\n"
        "uniform mediump mat4 viewMatrix;\n"
        "void main(void)\n"
        "{\n"
        "    gl_Position = projectionMatrix * viewMatrix * modelMatrix * vertex;\n"
        "    texc = texCoord;\n"
        "}\n";
    dots_vshader->compileSourceCode(dots_vsrc);
    
    QOpenGLShader *dots_fshader = new QOpenGLShader(QOpenGLShader::Fragment, this);
    const char *dots_fsrc =
        "#version 120\n"
        "uniform vec4 dcolour;\n"
        "uniform sampler2D texture;\n"
        "varying mediump vec4 texc;\n"
        "void main(void)\n"
        "{\n"
        "    vec2 uv = texc.st - vec2(0.5, 0.5);\n"
        "    float dist = sqrt(dot(uv, uv));\n"
        "    float t = smoothstep(0.3, 0.5, dist);\n"
        "    gl_FragColor = vec4(dcolour.xyz, 1 - t);\n"
        "}\n";
    dots_fshader->compileSourceCode(dots_fsrc);
    
    dots_program = std::shared_ptr<QOpenGLShaderProgram>(new QOpenGLShaderProgram);
    dots_program->addShader(dots_vshader);
    dots_program->addShader(dots_fshader);
    dots_program->bindAttributeLocation("vertex", prog_vert_att);
    dots_program->bindAttributeLocation("texCoord", prog_texcoord_att);
    dots_program->link();

    dots_program->bind();
    dots_program->setUniformValue("texture", 0);
}

void GL_image_panel_dots::paint_overlay(void) {
    dots_program->bind();
    dots_program->setUniformValue("viewMatrix", view);
    dots_program->setUniformValue("projectionMatrix", projection);
    
    dots_vbo.bind();
    dots_program->enableAttributeArray(prog_vert_att);
    dots_program->enableAttributeArray(prog_texcoord_att);
    dots_program->setAttributeBuffer(prog_vert_att, GL_FLOAT, 0, 3, 5 * sizeof(GLfloat));
    dots_program->setAttributeBuffer(prog_texcoord_att, GL_FLOAT, 3 * sizeof(GLfloat), 2, 5 * sizeof(GLfloat));
    
    constexpr double colours[3][3] = {
        {0.129, 0.627, 0.875},
        {0.604, 0.792, 0.329},
        {0.965, 0.651, 0.149}
    };
    
    if (dot_list[current_fname].size() > 0) {
        for (const auto& s: dot_list[current_fname]) {
            if (s.dot_no >= 1 && s.dot_no <= 3) {
                int j = s.dot_no - 1;
                draw_dot(s.p.x(), s.p.y(), colours[j][0], colours[j][1], colours[j][2]);    
            }
        }
    }
}


void GL_image_panel_dots::draw_dot(double x, double y, double r, double g, double b) {
    double lsf = 1.0;
    if (scale_factor < 1) {
        lsf = 0.5/scale_factor + 0.75;
    }
    if (scale_factor == 1.0) {
        lsf = 1.25;
    }
  
    QMatrix4x4 model = QMatrix4x4();
    model = QMatrix4x4();
    model.translate(x, y, 0.1);
    model.scale(1.2*lsf);
    dots_program->setUniformValue("modelMatrix", model);
    dots_program->setUniformValue("dcolour", 1, 1, 1, 1.0);
    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
    
    model = QMatrix4x4();
    model.translate(x, y, 0.2);
    model.scale(1.05*lsf);
    dots_program->setUniformValue("modelMatrix", model);
    dots_program->setUniformValue("dcolour", 0, 0, 0, 1.0);
    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
    
    
    model = QMatrix4x4();
    model.translate(x, y, 0.3);
    model.scale(lsf);
    dots_program->setUniformValue("modelMatrix", model);
    dots_program->setUniformValue("dcolour", r, g, b, 1.0);
    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
}

void GL_image_panel_dots::make_dots(void) {
    QVector<GLfloat> vd;
    
    const double rad = 10;
    vd.append(-rad);
    vd.append(-rad);
    vd.append(0.1);
    vd.append(0);
    vd.append(0);
    
    vd.append(rad);
    vd.append(-rad);
    vd.append(0.1);
    vd.append(1);
    vd.append(0);
    
    vd.append(rad);
    vd.append(rad);
    vd.append(0.1);
    vd.append(1);
    vd.append(1);
    
    vd.append(-rad);
    vd.append(rad);
    vd.append(0.1);
    vd.append(0);
    vd.append(1);
    

    dots_vbo.create();
    dots_vbo.bind();
    dots_vbo.allocate(vd.constData(), vd.count() * sizeof(GLfloat));
}

void GL_image_panel_dots::click_marker(QPoint pos, bool add) {
    
    int ix = pos.x() - imgsize.width()/2;
    int iy = pos.y() - imgsize.height()/2;
    
    int dot_no = 0;
    for (auto it=dot_list.begin(); it != dot_list.end(); it++) {
        dot_no += it->second.size();
    }
    
    if (add) {
        if (dot_no < 3) {
            dot_list[current_fname].push_back(Sfr_marker(ix, iy, dot_no + 1));
        } else {
            // just rebuild the whole list, skipping the one we are replacing
            std::map<std::string, vector<Sfr_marker> > new_cache;
            for (const auto& f: dot_list) {
                for (const auto& s: f.second) {
                    if (s.dot_no != dot_no) {
                        new_cache[f.first].push_back(s);
                    }
                }
            }
            dot_list = new_cache;
            dot_list[current_fname].push_back(Sfr_marker(ix, iy, dot_no));
        }
    } else {
        dot_list.clear();
        dot_list[current_fname].push_back(Sfr_marker(ix, iy, 1));
    }
}

void GL_image_panel_dots::mousePressEvent(QMouseEvent* event) { 
    click = event->pos();
    event->ignore();
}

static double sqr(double x) {
    return x*x;
}

void GL_image_panel_dots::mouseReleaseEvent(QMouseEvent* event) {
    if (event->button() == Qt::LeftButton) {
        double d = sqrt( sqr(event->pos().x() - click.x()) + sqr(event->pos().y() - click.y()) );
        
        if (d < 10) { // mouse "release" close enough to mouse "press" to consider it a click (rather than drag)
            
            QPoint img_coords = locate(click);
            
            bool shift_down = event->modifiers().testFlag(Qt::ShiftModifier);
            bool ctrl_down = event->modifiers().testFlag(Qt::ControlModifier);
            
            bool valid = app->edge_selected(img_coords.x(), img_coords.y(), ctrl_down, shift_down);
                
            if (valid) {
                click_marker(img_coords, shift_down);
                update();
            }
        } 
    }
    event->ignore();
}

void GL_image_panel_dots::mouseMoveEvent(QMouseEvent* event) {
    event->ignore();
}
