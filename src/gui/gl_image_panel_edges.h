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
#ifndef GL_IMAGE_PANEL_EDGES_H
#define GL_IMAGE_PANEL_EDGES_H

#include <limits>

#include "gl_image_panel.h"
#include "histogram_type.h"
#include "include/bayer.h"

#include <QPointF>
#include <array>
#include <assert.h>
#include <memory>

class GL_image_panel_edges : public GL_image_panel {
  Q_OBJECT
  
  public:
    explicit GL_image_panel_edges(QWidget *parent = 0);
    bool save_rois(const QString& fname);
    bool load_rois(const QString& fname);
    
    virtual void click_marker(QPoint pos, bool add=false) override;
    virtual void clear_overlay(void) override;
    virtual void initialize_overlay(void) override;
    virtual void paint_overlay(void) override;
    
    void mouseMoveEvent(QMouseEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;
    void mouseReleaseEvent(QMouseEvent* event) override;
    
    void line_endpoint(QPoint pos);
    void set_cfa_mask(Bayer::cfa_mask_t mask);
    
  signals:
    void update_edge_length(double edge_length);
    void update_histogram(histo_t dark, histo_t light);
    void enable_save_button(void);
    void disable_save_button(void);
    
  private:
    void make_dot_vbo();
    void draw_dot(double x, double y, double r, double g, double b);
    void draw_line(QPointF start, QPointF end, double r, double g, double b, double width = 1.0);
    void draw_roi(QPointF start, QPointF end, double r, double g, double b);
    void draw_box(QPointF start, QPointF end, float box_width, double r, double g, double b, double t = 1.0);
    void draw_close_symbol(QPointF pos, QPointF dir, double r, double g, double b);
    void check_roi_boxes_and_handles(QPointF img_coords);
    void broadcast_histogram(void);
    
    std::shared_ptr<QOpenGLShaderProgram> line_program;
    std::shared_ptr<QOpenGLShaderProgram> box_program;
    std::shared_ptr<QOpenGLShaderProgram> dots_program;
    
    QOpenGLBuffer edges_vbo;
    
    QOpenGLBuffer dot_vbo;
    
    QPointF img_centre;
    Bayer::cfa_mask_t cfa_mask = Bayer::cfa_mask_t::ALL;
    
    typedef enum {
        NONE,
        AFTER_FIRST_CLICK,
        MOVING,
        DRAGGING_HANDLE,
        ROI_SELECTED
    } state_t;
    
    QPoint release_pt;
    
    state_t state = NONE;
    
    class GL_roi {
      public:
        GL_roi(void);
        GL_roi(QPointF p1, QPointF p2);
        
        void add_point(QPointF p);
        QPointF& get(int i);
        double length(void) const;
        
        int handle_selected(QPointF p, double dist_thresh = 10.0);
        int box_selected(QPointF p, double dist_thresh = 10.0, double width_thresh = 28.0);
        
      private:
        size_t n = 0;  
        std::array<QPointF, 2> pts;
    };
    
    GL_roi* current_roi = nullptr;
    int current_roi_handle_idx = -1;
    vector<GL_roi> rois;
    
    class GL_closebox {
      public:
        GL_closebox(
            QPointF handle_a = QPointF(-1, -1), 
            QPointF handle_b = QPointF(-1, -1), 
            int img_width = std::numeric_limits<int>::max(), 
            int img_height = std::numeric_limits<int>::max()
        );
        
        bool selected(QPointF p, double dist_thresh = 10.0);
        
        QPointF get_pos(void) { return pos; }
        QPointF get_dir(void) { return dir; }
        
        bool is_valid(void) const { return valid; }
        void make_valid(void) { valid = true; }
        void make_invalid(void) { valid = false; }
        
      private:
        QPointF pos;
        QPointF dir;
        bool valid = false;
    };
    
    GL_closebox closebox;
};

#endif

