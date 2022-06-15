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
#ifndef SETTINGS_IO_TAB_H
#define SETTINGS_IO_TAB_H

#include <QtWidgets>

class Settings_io_tab : public QWidget {
  Q_OBJECT
  
  public:
    Settings_io_tab(QWidget *parent = nullptr);

    QLabel* arguments_label;
    QLineEdit* arguments_line;
    QLabel* threshold_label;
    QLineEdit* threshold_line;
    QLabel* pixsize_label;
    QLineEdit* pixsize_line;
    QLabel* pixsize_units;

    QLabel* cache_label;
    QLineEdit* cache_line;
    QLabel* contrast_label;
    QLineEdit* contrast_line;
    QLabel* lpmm_label;
    QLineEdit* lp1_line;
    QLineEdit* lp2_line;
    QLineEdit* lp3_line;

    QCheckBox* cb_linear_gamma;
    QCheckBox* cb_annotation;
    QCheckBox* cb_profile;
    QCheckBox* cb_grid;
    QCheckBox* cb_orientation;
    QCheckBox* cb_lensprofile;
    QCheckBox* cb_autocrop;
    QCheckBox* cb_lpmm;
    QCheckBox* cb_gnuplot_scaled;
    QCheckBox* cb_lensprofile_fixed;
    QCheckBox* cb_ca_active;
    QCheckBox* cb_fullsfr;
    QCheckBox* cb_nosmoothing;

    QComboBox* box_colour; // Bayer channel, actually
    QComboBox* box_esf_model;
    QComboBox* box_ca_type;
    QLabel* bayer_label;
    QLabel* esf_model_label;
    QLabel* ca_type_label;

    QLabel* zscale_label;
    QSlider* zscale_slider;

    QCheckBox* cb_surface_max;
    QLineEdit* surface_max_value;
    QLabel* surface_max_units;

    QPalette normal_pal;
    QPalette disabled_pal;

    static const QString setting_threshold;
    static const QString setting_threshold_default;
    static const QString setting_pixsize;
    static const QString setting_pixsize_default;
    static const QString setting_bayer;
    static const QString setting_esf_model;
    static const QString setting_linear_gamma;
    static const Qt::CheckState setting_linear_gamma_default;
    static const QString setting_annotation;
    static const Qt::CheckState setting_annotation_default;
    static const QString setting_profile;
    static const Qt::CheckState setting_profile_default;
    static const QString setting_grid;
    static const Qt::CheckState setting_grid_default;
    static const QString setting_lensprofile;
    static const Qt::CheckState setting_lensprofile_default;
    static const QString setting_orientation;
    static const Qt::CheckState setting_orientation_default;
    static const QString setting_autocrop;
    static const Qt::CheckState setting_autocrop_default;
    static const QString setting_lpmm;
    static const Qt::CheckState setting_lpmm_default;
    static const QString setting_gnuplot_scaled;
    static const Qt::CheckState setting_gnuplot_scaled_default;
    static const QString setting_zscale;
    static const int setting_zscale_default;
    static const QString setting_cache;
    static const double setting_cache_default;
    static const QString setting_mtf_contrast;
    static const double setting_mtf_contrast_default;
    static const QString setting_lp1;
    static const QString setting_lp1_default;
    static const QString setting_lp2;
    static const QString setting_lp2_default;
    static const QString setting_lp3;
    static const QString setting_lp3_default;
    static const QString setting_arguments;
    static const QString setting_lensprofile_fixed;
    static const Qt::CheckState setting_lensprofile_fixed_default;
    static const QString setting_surface_max_value;
    static const QString setting_surface_max_value_default;
    static const QString setting_surface_max_flag;
    static const Qt::CheckState setting_surface_max_flag_default;
    static const QString setting_ca_active;
    static const Qt::CheckState setting_ca_active_default;
    static const QString setting_ca_type;
    static const int setting_ca_type_default;
    static const QString setting_fullsfr;
    static const Qt::CheckState setting_fullsfr_default;
    static const QString setting_nosmoothing;
    static const Qt::CheckState setting_nosmoothing_default;

  private:
  public slots:
    void lpmm_toggled();

};

#endif

