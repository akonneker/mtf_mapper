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
#ifndef GL_IMAGE_PANEL_H
#define GL_IMAGE_PANEL_H

#include <vector>
using std::vector;
#include <map>
using std::pair;
using std::make_pair;
#include <string>
using std::string;

#include <cmath>
#include <algorithm>

#include <opencv2/opencv.hpp>

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLBuffer>
#include <QMatrix4x4>

#include "sfr_marker.h"
#include "cache_entry.h"
#include "image_viewport.h"

QT_FORWARD_DECLARE_CLASS(QOpenGLShaderProgram)
QT_FORWARD_DECLARE_CLASS(QOpenGLTexture)

class GL_image_panel : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    explicit GL_image_panel(QWidget *parent = 0);
    ~GL_image_panel();

    QSize minimumSizeHint() const override;
    QSize sizeHint() const override;
    QSize img_size() const { return imgsize; }
    
    void move(int nx, int ny);
    QPoint zoom(int step, int x, int y);
    QPoint locate(QPoint pos);
    bool load_image(const QString& fname);
    bool load_image(QImage& qimg);
    void set_cache_size(uint64_t size);
    void set_cache_state(bool s) { cache_enabled = s; }
    int image_channels(void) const { return img_channels; }
    int image_depth(void) const { return img_depth; }
    
    virtual void click_marker(QPoint pos, bool add=false) = 0;
    virtual void clear_overlay(void) = 0;
    virtual void initialize_overlay(void) = 0;
    virtual void paint_overlay(void) = 0;
    
    void set_default_image(QImage* qimg) { default_image = qimg; }
    Image_viewport get_viewport(void) { return vp; }
    
    void mouseMoveEvent(QMouseEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;
    void mouseReleaseEvent(QMouseEvent* event) override;
    
    cv::Mat get_cv_img(void) const { return current_cv_image; }
    void set_max_scale_factor(double factor) { max_scale_factor = factor; }
    void set_gamma(float gamma) { gamma_value = gamma; }
    int get_img_min(void) { return img_minmax.first; }
    int get_img_max(void) { return img_minmax.second; }
    
    static int program_counter;

protected:
    void initializeGL() override;
    void paintGL() override;
    void resizeGL(int width, int height) override;
    bool load_image(cv::Mat cvimg);
    void reset_scroll_range(void);
    void trim_cache(void);
    
    QMatrix4x4 view;
    QMatrix4x4 projection;
    std::string current_fname;
    
    QSize imgsize = QSize(0,0);
    double scale_factor = 1.0;
    double max_scale_factor = 2.0;
    
    QPoint located_pos;
    QPoint click;
    
    static constexpr int prog_vert_att = 0;
    static constexpr int prog_texcoord_att = 1;
    
private:
    vector<QOpenGLTexture*> textures;
    QOpenGLShaderProgram* program;
    QOpenGLBuffer vbo;
    
    std::map<std::string, Cache_entry> image_cache;
    uint64_t image_cache_size = 0;
    uint64_t max_image_cache_size = uint64_t(1) << 30;
    bool cache_enabled = true;
    
    int img_channels = 1;
    int img_depth = 8;
    bool img_assumed_linear = false;
    std::pair<int, int> img_minmax{-1,-1};
    
    Image_viewport vp;
    
    QImage* default_image = nullptr;
    cv::Mat current_cv_image;
    
    float gamma_value = 1.0;
};

#endif
