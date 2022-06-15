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
#include <QtWidgets> 
#include "settings_dialog.h"

#include "common.h"

#include <iostream>
using std::cout;
using std::endl;



Settings_dialog::Settings_dialog(QWidget *parent ATTRIBUTE_UNUSED)
 : settings("mtfmapper", "mtfmapper"), gnuplot_img_width(1024)
{

    io = new Settings_io_tab;
    helpers = new Settings_helpers_tab;
    distortion = new Settings_distortion_tab(io->cb_lpmm);
        
    accept_button = new QPushButton(tr("&Accept"), this);
    accept_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    cancel_button = new QPushButton(tr("&Cancel"), this);
    cancel_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
        
    tab_widget = new QTabWidget;
    tab_widget->addTab(io, "Input / Output");
    tab_widget->addTab(distortion, "Distortion");
    tab_widget->addTab(helpers, "Helpers");

    QHBoxLayout* button_layout = new QHBoxLayout;
    button_layout->addWidget(accept_button);
    button_layout->addWidget(cancel_button);
    button_layout->addStretch();
    
    QVBoxLayout* vlayout = new QVBoxLayout;
    vlayout->addWidget(tab_widget);
    vlayout->addLayout(button_layout);
    
    connect(accept_button, SIGNAL(clicked()), this, SLOT( save_and_close() ));
    connect(cancel_button, SIGNAL(clicked()), this, SLOT( close() ));
    connect(io->cb_lpmm, SIGNAL(clicked()), this, SLOT( lpmm_toggled() ));
    
    setLayout(vlayout);
    setWindowTitle("Preferences");
}

void Settings_dialog::open() {

    show();
}

QString Settings_dialog::get_argument_string(bool focus_mode) {
    QString args = QString("-t %1").arg(io->threshold_line->text());
    if (io->cb_linear_gamma->checkState()) {
        args = args + QString(" -l");
    }
    
    if (io->cb_annotation->checkState()) {
        args = args + QString(" -a");
        args = args + QString(" -v 2");
    }
    
    if (io->cb_profile->checkState()) {
        args = args + QString(" -p");
    }
    
    if (io->cb_grid->checkState()) {
        args = args + QString(" -s");
    }
    
    if (io->cb_lensprofile->checkState()) {
        args = args + QString(" --lensprofile");
    }
    
    if (io->cb_orientation->checkState()) {
        args = args + QString(" --chart-orientation");
    }
    
    if (io->cb_ca_active->checkState()) {
        args = args + QString(" --ca");
    }
    
    if (io->cb_autocrop->checkState()) {
        args = args + QString(" --autocrop");
    }

    if (io->cb_lpmm->checkState()) {
        args = args + QString(" --pixelsize " + io->pixsize_line->text());
    }
    
    switch(io->box_colour->currentIndex()) {
    case 0: break;
    case 1: args = args + QString(" --bayer red"); break;
    case 2: args = args + QString(" --bayer green"); break;
    case 3: args = args + QString(" --bayer blue"); break;
    }
    
    switch(io->box_esf_model->currentIndex()) {
    case 0: args = args + QString(" --esf-model kernel"); break;
    case 1: args = args + QString(" --esf-model loess"); break;
    }
    
    if (io->box_ca_type->currentIndex() == 1) {
        args = args + QString(" --ca-fraction");
    }
    
    if (distortion->rb_lens_pw_quad->isChecked()) {
        args = args + QString(" --esf-sampler piecewise-quadratic");
    }
    if (distortion->rb_lens_quad->isChecked()) {
        args = args + QString(" --esf-sampler quadratic");
    }
    if (distortion->rb_lens_none->isChecked()) {
        args = args + QString(" --esf-sampler line");
    }
    if (distortion->rb_lens_radial->isChecked()) {
        args = args + QString(" --optimize-distortion");
    }
    if (distortion->rb_lens_equiangular->isChecked()) {
        args = args + QString(" --equiangular %1").arg(distortion->ea_f_line->text());
    }
    if (distortion->rb_lens_stereo->isChecked()) {
        args = args + QString(" --stereographic %1").arg(distortion->sg_f_line->text());
    }
    
    if (io->cb_gnuplot_scaled->checkState()) {
        args = args + QString(" --gnuplot-width %1").arg(gnuplot_img_width);
    }
    
    if (io->cb_lensprofile_fixed->checkState()) {
        args = args + QString(" --lensprofile-fixed-size");
    }
    
    if (io->cb_surface_max->checkState() && io->surface_max_value->text().length() > 0) {
        double surface_max = io->surface_max_value->text().toDouble();
        if (!io->cb_lpmm->checkState() && surface_max > 1.0) {
            surface_max = 0.0;
        }
        args = args + QString(" --surface-max %1").arg(surface_max);
    }
    
    args = args + QString(" --zscale %1").arg(io->zscale_slider->value()/20.0);
    
    args = args + QString(" --mtf %1").arg(io->contrast_line->text());
    
    // only add --lpX arguments if we are in lp/mm mode
    if (io->cb_lpmm->checkState()) {
        if (io->lp1_line->text().length() > 0) {
            args = args + QString(" --lp1 %1").arg(io->lp1_line->text());
        }
        
        if (io->lp2_line->text().length() > 0) {
            args = args + QString(" --lp2 %1").arg(io->lp2_line->text());
        }
        
        if (io->lp3_line->text().length() > 0) {
            args = args + QString(" --lp3 %1").arg(io->lp3_line->text());
        }
    }

    if (io->cb_fullsfr->checkState()) {
        args = args + QString(" --full-sfr");
    }

    if (io->cb_nosmoothing->checkState()) {
        args = args + QString(" --nosmoothing");
    }
    
    if (!focus_mode) {
        args = args + QString(" -q");
    }

    args = args + QString(" --jpeg");
    
    args = args + QString(" %1").arg(io->arguments_line->text());
    
    return args;
}

void Settings_dialog::save_and_close() {
    helpers->check_gnuplot_binary();
    helpers->check_exiv2_binary();
    helpers->check_dcraw_binary();
    check_mtf_lower();
    settings.setValue(Settings_io_tab::setting_threshold, io->threshold_line->text());
    settings.setValue(Settings_io_tab::setting_pixsize, io->pixsize_line->text());
    settings.setValue(Settings_io_tab::setting_mtf_contrast, io->contrast_line->text());
    settings.setValue(Settings_io_tab::setting_linear_gamma, io->cb_linear_gamma->checkState());
    settings.setValue(Settings_io_tab::setting_annotation, io->cb_annotation->checkState());
    settings.setValue(Settings_io_tab::setting_profile, io->cb_profile->checkState());
    settings.setValue(Settings_io_tab::setting_lpmm, io->cb_lpmm->checkState());
    settings.setValue(Settings_io_tab::setting_grid, io->cb_grid->checkState());
    settings.setValue(Settings_io_tab::setting_lensprofile, io->cb_lensprofile->checkState());
    settings.setValue(Settings_io_tab::setting_orientation, io->cb_orientation->checkState());
    settings.setValue(Settings_io_tab::setting_ca_active, io->cb_ca_active->checkState());
    settings.setValue(Settings_io_tab::setting_autocrop, io->cb_autocrop->checkState());
    settings.setValue(Settings_io_tab::setting_gnuplot_scaled, io->cb_gnuplot_scaled->checkState());
    settings.setValue(Settings_io_tab::setting_lensprofile_fixed, io->cb_lensprofile_fixed->checkState());
    settings.setValue(Settings_io_tab::setting_surface_max_flag, io->cb_surface_max->checkState());
    settings.setValue(Settings_helpers_tab::setting_gnuplot, helpers->gnuplot_line->text());
    settings.setValue(Settings_helpers_tab::setting_exiv, helpers->exiv_line->text());
    settings.setValue(Settings_helpers_tab::setting_dcraw, helpers->dcraw_line->text());
    settings.setValue(Settings_helpers_tab::setting_dcraw_emu, helpers->dcraw_emu_line->text());
    settings.setValue(Settings_helpers_tab::setting_unprocessed_raw, helpers->unproc_raw_line->text());
    settings.setValue(Settings_helpers_tab::setting_raw_developer, helpers->box_raw_developer->currentIndex());
    settings.setValue(Settings_io_tab::setting_zscale, io->zscale_slider->value());
    settings.setValue(Settings_io_tab::setting_cache, io->cache_line->text());
    settings.setValue(Settings_io_tab::setting_lp1, io->lp1_line->text());
    settings.setValue(Settings_io_tab::setting_lp2, io->lp2_line->text());
    settings.setValue(Settings_io_tab::setting_lp3, io->lp3_line->text());
    settings.setValue(Settings_io_tab::setting_arguments, io->arguments_line->text());
    settings.setValue(Settings_io_tab::setting_bayer, io->box_colour->currentIndex());
    settings.setValue(Settings_io_tab::setting_esf_model, io->box_esf_model->currentIndex());
    settings.setValue(Settings_io_tab::setting_ca_type, io->box_ca_type->currentIndex());
    settings.setValue(Settings_io_tab::setting_surface_max_value, io->surface_max_value->text());
    settings.setValue(Settings_io_tab::setting_fullsfr, io->cb_fullsfr->checkState());
    settings.setValue(Settings_io_tab::setting_nosmoothing, io->cb_nosmoothing->checkState());
    
    if (distortion->rb_lens_pw_quad->isChecked()) {
        settings.setValue(Settings_distortion_tab::setting_lens, 0);
    }
    if (distortion->rb_lens_quad->isChecked()) {
        settings.setValue(Settings_distortion_tab::setting_lens, 1);
    }    
    if (distortion->rb_lens_none->isChecked()) {
        settings.setValue(Settings_distortion_tab::setting_lens, 2);
    }
    if (distortion->rb_lens_radial->isChecked()) {
        settings.setValue(Settings_distortion_tab::setting_lens, 3);
    }
    if (distortion->rb_lens_equiangular->isChecked()) {
        settings.setValue(Settings_distortion_tab::setting_lens, 4);
        io->cb_lpmm->setCheckState(Qt::Checked);
        settings.setValue(Settings_io_tab::setting_lpmm, io->cb_lpmm->checkState());
    }
    if (distortion->rb_lens_stereo->isChecked()) {
        settings.setValue(Settings_distortion_tab::setting_lens, 5);
        io->cb_lpmm->setCheckState(Qt::Checked);
        settings.setValue(Settings_io_tab::setting_lpmm, io->cb_lpmm->checkState());
    }
    
    set_cache_size(settings.value(Settings_io_tab::setting_cache, io->setting_cache_default).toInt());
    settings_saved();
    
    close();
}

void Settings_dialog::check_mtf_lower(void) {
    double min_contrast = io->contrast_line->text().toDouble();
    if (min_contrast < 10) {
        QMessageBox::warning(
            this, 
            QString("Low MTF-XX value warning"), 
            QString(
                "Selecting MTF-XX values below 10 may result in unexpected behaviour (some edges may not be"
                " detected, or results may be very sensitive to noise)."
            )
        );
    }
}


void Settings_dialog::set_gnuplot_img_width(int w) {
    gnuplot_img_width = w < 1024 ? 1024 : w;
}


void Settings_dialog::lpmm_toggled() {
    double resolution_scale = 1.0;
    if (io->pixsize_line->text().length() > 0) {
        resolution_scale = 1000.0/io->pixsize_line->text().toDouble();
    }
    if (io->cb_lpmm->checkState() == Qt::Unchecked) {
        io->surface_max_units->setText("c/p");
        char buffer[100];
        sprintf(buffer, "%.3lf", io->surface_max_value->text().toDouble() / resolution_scale);
        io->surface_max_value->setText(buffer);
    } else {
        io->surface_max_units->setText("lp/mm");
        char buffer[100];
        sprintf(buffer, "%.2lf", io->surface_max_value->text().toDouble() * resolution_scale);
        io->surface_max_value->setText(buffer);
    }
}

QString Settings_dialog::peek_argument_line(void) const {
    return io->arguments_line->text();
}

void Settings_dialog::reset_argument_line(void) {
    io->arguments_line->setText(QString(""));
    settings.setValue(io->setting_arguments, io->arguments_line->text());
}
