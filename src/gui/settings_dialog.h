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
#ifndef SETTINGS_DIALOG_H
#define SETTINGS_DIALOG_H

#include <QtWidgets>
#include <QSettings>


#include "settings_io_tab.h"
#include "settings_helpers_tab.h"
#include "settings_distortion_tab.h"

class Settings_dialog : public QDialog 
{
  Q_OBJECT
  
  public:
    Settings_dialog(QWidget *parent);
    QString get_argument_string(bool focus_mode);
    void check_mtf_lower(void);
    QString peek_argument_line(void) const;
    void reset_argument_line(void);
    
    
    QSettings   settings;
    
    QPushButton* accept_button;
    QPushButton* cancel_button;
    
    QTabWidget* tab_widget;

    Settings_helpers_tab* helpers;
    Settings_distortion_tab* distortion;
    Settings_io_tab* io;
    
    int gnuplot_img_width;
    
  signals:
    void argument_string(QString s);  
    void set_cache_size(int);
    void settings_saved(void);
    
  public slots:
    void open();
    void save_and_close();
    void set_gnuplot_img_width(int w);
    void lpmm_toggled();
};

#endif

