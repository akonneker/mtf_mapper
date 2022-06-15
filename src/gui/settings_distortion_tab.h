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
#ifndef SETTINGS_DISTORTION_TAB_H
#define SETTINGS_DISTORTION_TAB_H

#include <QtWidgets>

class Settings_distortion_tab : public QWidget {
    Q_OBJECT

  public:
    Settings_distortion_tab(QCheckBox* parent_lp_mm, QWidget* parent = nullptr);

    QLabel* ea_f_label;
    QLineEdit* ea_f_line;
    QLabel* sg_f_label;
    QLineEdit* sg_f_line;

    QRadioButton* rb_lens_pw_quad;
    QRadioButton* rb_lens_quad;
    QRadioButton* rb_lens_none;
    QRadioButton* rb_lens_radial;
    QRadioButton* rb_lens_equiangular;
    QRadioButton* rb_lens_stereo;

    QCheckBox* parent_lp_mm;

    static const QString setting_ea_f;
    static const QString setting_ea_f_default;
    static const QString setting_sg_f;
    static const QString setting_sg_f_default;
    static const QString setting_lens;

  private:

  public slots:
    void equiangular_toggled();
    void stereographic_toggled();

};

#endif

