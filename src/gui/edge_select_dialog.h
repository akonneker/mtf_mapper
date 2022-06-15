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
#ifndef EDGE_SELECT_DIALOG_H
#define EDGE_SELECT_DIALOG_H

#include <QDialog>

#include "gl_image_viewer.h"
#include "gl_image_panel_edges.h"
#include "manual_roi_help_dialog.h"
#include "histo_widget.h"
#include "histogram_type.h"
#include "include/bayer.h"
#include <memory>

class QPushButton;
class QTextEdit;
class QGroupBox;


class Edge_select_dialog : public QDialog {
  Q_OBJECT
  
  public:
    Edge_select_dialog(QWidget* parent);
    
    bool load_image(QString img_name, 
        QStringList arguments = QStringList(), 
        QString raw_filename = "");
        
    GL_image_viewer* get_viewer(void) { return img_viewer; }
    GL_image_panel* get_panel(void) { return img_panel; }
    
    void set_roi_file(const QString& fname) { roi_file = fname; }
    bool autoload_roi(void) const { return rb_autoload_roi->isChecked(); }
    bool load_rois(QString fname);
    bool save_rois(QString fname);

  private:
    void extract_bayer_info(const QStringList& arguments);
    GL_image_viewer*  img_viewer;
    GL_image_panel_edges*   img_panel;
    QPushButton* cancel_button;
    QPushButton* accept_button;
    QPushButton* help_button;
    QPushButton* load_button;
    QPushButton* save_button;
    QPushButton* apply_all_button;
    QPushButton* abort_all_button;
    QPushButton* clear_roi_button;
    QWidget* parent = nullptr;
    QLabel* img_filename;
    QLabel* text_img_filename;
    QLabel* img_type;
    QLabel* text_img_type;
    QLabel* img_progress;
    QLabel* text_img_progress;
    QLabel* edge_length;
    QLabel* text_edge_length;
    QLabel* img_max_val;
    QLabel* text_img_max_val;
    QToolButton* gamma_switch;
    QRadioButton* rb_autoload_roi;
    QRadioButton* rb_blank_roi;
    Histo_widget* histogram;
    std::shared_ptr<QImage> icon_image;

    QString roi_file;
    Bayer::bayer_t bayer_channel;
    Bayer::cfa_pattern_t cfa_pattern;
    Bayer::cfa_mask_t cfa_mask;

    Manual_roi_help_dialog* help_dialog;
  
  signals:
    void start_manual_queue_processing(QString filename);
    void manual_queue_skip_all();

  public slots:
    void open();
    void export_roi();
    void show_help();
    void load_roi();
    void save_roi();
    void update_queue_size(int queue_size);
    void update_edge_length(double edge_length);
    void update_histogram(histo_t dark, histo_t light);
    void disable_save_button();
    void enable_save_button();
    void apply_all_button_clicked();
    void abort_all_button_clicked();
    void gamma_state(bool);
    void clear_roi_button_clicked();
};

#endif

