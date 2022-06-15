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
#ifndef SFR_DIALOG_H
#define SFR_DIALOG_H

#include <QDialog>
#include <QtCharts>
using namespace QtCharts;

#include "sfr_entry.h"
#include "sfr_chartview.h"
#include "entry_view.h"

#include <vector>
using std::vector;
#include <memory>

class Sfr_dialog : public QDialog {
  Q_OBJECT
  
  public:
    Sfr_dialog(QWidget *parent, const QRect initial_geom = QRect(-1, -1, 0, 0));
        
    void replace_entry(const Sfr_entry& entry);
    void add_entry(const Sfr_entry& entry);
    void notify_mouse_position(double value, bool click=false);
    void clear(void);
    void keyPressEvent(QKeyEvent* event);
    void keyReleaseEvent(QKeyEvent* event);
    int pip_number(void) const { return entries.size(); }
    
  signals:
    void sfr_dialog_closed();
    void send_geometry(QRect r);
  
  protected:
    void paintEvent(QPaintEvent* event) override;
    void reject(void) override;

  private:
    void set_label_background(QLabel* label, const string& condition);
    int xscale(void);
    void hide_xscale(void);
    void show_xscale(void);
  
    vector<Sfr_entry> entries;
    QChart* chart;
    Sfr_chartview* chart_view;
    QValueAxis* x_axis;
    QValueAxis* y_axis;
    vector<QLineSeries*> series; 
    vector<QGraphicsSimpleTextItem*> mtf50_text;
    double cursor_domain_value;
    vector<vector<QLabel*>> table_labels;
    QLabel* x_label;
    QLabel* x_label_value;
    QGraphicsRectItem* mtf50_rect;
    QGridLayout* label_layout;
    QPushButton* save_img_button;
    QPushButton* save_data_button;
    std::shared_ptr<QIcon> mtfmapper_logo;
    QAtomicInt repainting;
    bool lock_cursor = false;
    Entry_view view;
    QComboBox* box_view;
    Qt::KeyboardModifiers last_modifier;
    QCheckBox* lp_mm_cb;
    
    QLabel* x_scale_label;
    QSlider* x_scale_slider = nullptr;
    vector<int> x_scale_levels;
    
  public slots:
    bool update_lp_mm_mode(void);
    void save_image(void);
    void save_data(void);
    void plot_type_changed(int index);
    void lp_mm_toggled(void);
    void x_scale_changed(int);
};


#endif

