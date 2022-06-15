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
#ifndef MTFMAPPER_APP_H
#define MTFMAPPER_APP_H
 
#include <QMainWindow> 
#include <QModelIndex>
#include <QFileDialog>
#include <QStringList>
#include <QList>
#include <QtConcurrent/QtConcurrent>
#include "settings_dialog.h"
#include "worker_thread.h"
#include "gl_image_viewer.h"
#include "gl_image_panel_dots.h"
#include "exiv2_property.h"
#include "about_dialog.h"
#include "help_dialog.h"
#include "processing_command.h"
#include "processor_state.h"

#include "sfr_entry.h"
#include "sfr_dialog.h"

#include "edge_select_dialog.h"

#include <opencv2/opencv.hpp>

#include <vector>
using std::vector;
#include <queue>
using std::queue;
#include <memory>

class QPushButton;
class QLabel;
class QGroupBox;
class QGraphicsScene;
class QGraphicsView;
class QStringList;
class QSpinBox;
class QGraphicsPixmapItem;
class QTreeView;
class QThread;
class QCheckBox;
class QProgressBar;
class QSplitter;

class mtfmapper_app : public QMainWindow
{
  Q_OBJECT
        
  public:
    mtfmapper_app(QWidget *parent = 0);
    virtual ~mtfmapper_app(void);
    void dragEnterEvent(QDragEnterEvent* evt) override;
    void dropEvent(QDropEvent* evt) override;
    
  protected:
    void closeEvent(QCloseEvent* event) override;
  
  private:
    void create_actions(void);
    void view_image(const QString& fname);
    void display_exif_properties(int index);
    void clear_temp_files(void);
    void check_and_purge_stale_temp_files(void);
    void check_if_helpers_exist(void);
    void save_action(bool subset = false);
    void open_action(bool roi = false, bool focus = false, bool imatest = false, bool manual_roi = false);
    void update_pip_icons(void);
    
    QMenu*          file_menu;
    QMenu*          settings_menu;
    QMenu*          help_menu;
    QAction*        open_act;
    QAction*        open_roi_act;
    QAction*        open_manual_roi_act;
    QAction*        open_focus_act;
    QAction*        open_imatest_act;
    QAction*        exit_act;
    QAction*        prefs_act;
    QAction*        about_act;
    QAction*        help_act;
    
    QTreeView*      datasets;
    
    QLabel*         img_comment_label;
    QLabel*         img_comment_value;

    QLabel*         af_ft_label;
    QLabel*         af_ft_value;

    QLabel*         focal_length_label;
    QLabel*         focal_length_value;

    QLabel*         focus_distance_label;
    QLabel*         focus_distance_value;
    
    QGroupBox*      horizgroup;
    
    Settings_dialog*    settings;
    About_dialog*       about;
    Help_dialog*        help;
    
    QStringList     input_files;
    QStringList     dataset_files;
    QStringList     tempfiles_to_delete;
    QList<std::shared_ptr<Exiv2_property>>  exif_properties;
    
    QThread*        worker_thread;
    Worker_thread   processor;
    
    QStandardItemModel   dataset_contents;
    
    QStandardItem*  current_dataset_item;
    
    QCheckBox*      tb_img_annotated;
    QCheckBox*      tb_img_profile;
    QCheckBox*      tb_img_gridimg;
    QCheckBox*      tb_img_lensprofile;
    QCheckBox*      tb_img_orientation;
    QCheckBox*      tb_img_ca;
    
    QProgressBar*   progress;
    QPushButton*    abort_button;

    QPushButton*    clear_button;
    QPushButton*    save_button;
    QPushButton*    save_subset_button;
    
    QFrame*           img_frame;
    GL_image_viewer*  img_viewer;
    GL_image_panel*   img_panel;
    
    QSplitter*      splitter;

    QIcon* mtfmapper_logo;
    QImage* icon_image;
    
    Sfr_dialog*     sfr_dialog;
    vector<Sfr_entry> sfr_list;
    std::pair<QString, QStandardItem*> active_annotated_filename;
    std::map<QString, std::pair<int, QStandardItem*>> sfr_pip_map;

    QString zoom_scroll_tt;
    QString annotated_tt;
    
    std::mutex mq_mutex;
    bool mq_busy = false;
    Edge_select_dialog* edge_select_dialog;
    QRect esd_geom = QRect(-1, -1, 0, 0);
    queue<Processing_command> manual_roi_commands;
    
    bool manual_queue_active = false;
    QString manual_queue_roi_filename;
    bool manual_queue_skip_all = false;
    
    QRect sfr_dialog_geom = QRect(-1, -1, 0, 0);
    
  signals:
    void submit_batch(Processor_state state);
    void manual_roi_queue_size(int queue_size);

  public slots:
    void open_auto();
    void open_roi();
    void open_manual_roi();
    void open_focus();
    void open_imatest_chart();
    void dataset_selected(const QModelIndex&);
    void dataset_selected_changed(const QModelIndex&, const QModelIndex&);
    void parent_item(QString s, QString f, QString tempdir);
    void child_item(QString s, QString f);
    void close_item(void);
    void item_for_deletion(QString s);
    void populate_exif_info_from_file(QString s, QString tempdir);
    void add_manual_roi_file(const Processing_command& command);
    void update_progress(int val);
    void enable_manual_queue_processing(QString filename);
    void enable_manual_queue_skip();
    void sfr_dialog_closed(QRect r);
  
    bool edge_selected(int px, int py, bool crtl_down, bool shift_down);
    

    void processor_completed(void);
    void enable_clear_button(void);
    void disable_clear_button(void);
    void clear_button_pressed(void);
    void enable_save_button(void);
    void disable_save_button(void);
    void save_button_pressed(void);
    void save_subset_button_pressed(void);
    
    void mtfmapper_call_failed(Worker_thread::failure_t failure, const QString& input_file);

    void enable_file_open(void);
    void disable_file_open(void);
    
    void set_cache_size(int);
    void settings_saved(void);
};
                               
                                
#endif
