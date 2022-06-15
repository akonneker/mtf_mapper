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

#include "settings_distortion_tab.h"

#include "common.h"
#include "nonempty_validator.h"

#include <QtGlobal>
#if QT_VERSION < QT_VERSION_CHECK(5,11,0)
#define horizontalAdvance width
#endif

const QString Settings_distortion_tab::setting_ea_f = "setting_equiangular_f";
const QString Settings_distortion_tab::setting_ea_f_default = "16.0";
const QString Settings_distortion_tab::setting_sg_f = "setting_stereographic_f";
const QString Settings_distortion_tab::setting_sg_f_default = "8.0";
const QString Settings_distortion_tab::setting_lens = "setting_lens_correction_type";

Settings_distortion_tab::Settings_distortion_tab(QCheckBox* parent_lp_mm, QWidget *parent ATTRIBUTE_UNUSED)
: parent_lp_mm(parent_lp_mm) {
    QSettings settings("mtfmapper", "mtfmapper");

    QFontMetrics fm(QApplication::font());
    int reasonable_width = fm.horizontalAdvance("1048576000");

    rb_lens_pw_quad = new QRadioButton("piecewise-quadratic");
    rb_lens_quad = new QRadioButton("quadratic");
    rb_lens_none = new QRadioButton("none");
    rb_lens_radial = new QRadioButton("radial");
    rb_lens_equiangular = new QRadioButton("equiangular");
    rb_lens_stereo = new QRadioButton("stereographic");

    QDoubleValidator* dv_f = new Nonempty_DoubleValidator(0.01, 999.0, 2, 16.0, this);
    ea_f_label = new QLabel(tr("focal length"), this);
    ea_f_line = new QLineEdit(this);
    ea_f_line->setMaxLength(5);
    ea_f_line->setText(settings.value(setting_ea_f, setting_ea_f_default).toString());
    ea_f_line->setValidator(dv_f);
    ea_f_line->setMaximumWidth(reasonable_width);
    sg_f_label = new QLabel(tr("focal length"), this);
    sg_f_line = new QLineEdit(this);
    sg_f_line->setMaxLength(5);
    sg_f_line->setText(settings.value(setting_sg_f, setting_sg_f_default).toString());
    sg_f_line->setValidator(dv_f);
    sg_f_line->setMaximumWidth(reasonable_width);

    rb_lens_pw_quad->setChecked(false);
    rb_lens_quad->setChecked(false);
    rb_lens_none->setChecked(false);
    rb_lens_radial->setChecked(false);
    rb_lens_equiangular->setChecked(false);
    rb_lens_stereo->setChecked(false);

    switch (settings.value(setting_lens, 0).toInt()) {
    case 0: rb_lens_pw_quad->setChecked(true); break;
    case 1: rb_lens_quad->setChecked(true); break;
    case 2: rb_lens_none->setChecked(true); break;
    case 3: rb_lens_radial->setChecked(true); break;
    case 4: rb_lens_equiangular->setChecked(true); break;
    case 5: rb_lens_stereo->setChecked(true); break;
    }

    
    QGridLayout* lens_layout = new QGridLayout;
    lens_layout->addWidget(rb_lens_pw_quad, 0, 0);
    lens_layout->addWidget(rb_lens_quad, 1, 0);
    lens_layout->addWidget(rb_lens_none, 2, 0);
    lens_layout->addWidget(rb_lens_radial, 3, 0);
    lens_layout->addWidget(rb_lens_equiangular, 4, 0);
    lens_layout->addWidget(rb_lens_stereo, 5, 0);
    lens_layout->addWidget(ea_f_label, 4, 1);
    lens_layout->addWidget(ea_f_line, 4, 2);
    lens_layout->addWidget(sg_f_label, 5, 1);
    lens_layout->addWidget(sg_f_line, 5, 2);

    QGroupBox* gb = new QGroupBox("Lens distortion correction");
    gb->setLayout(lens_layout);

    QVBoxLayout* tab_v_layout = new QVBoxLayout;
    tab_v_layout->addWidget(gb);
    tab_v_layout->addStretch(1);

    QHBoxLayout* tab_h_layout = new QHBoxLayout;
    tab_h_layout->addLayout(tab_v_layout);
    tab_h_layout->addStretch(1);

    setLayout(tab_h_layout);

    connect(rb_lens_equiangular, SIGNAL(clicked()), this, SLOT(equiangular_toggled()));
    connect(rb_lens_stereo, SIGNAL(clicked()), this, SLOT(stereographic_toggled()));

}

void Settings_distortion_tab::equiangular_toggled() {
    if (rb_lens_equiangular->isChecked() && parent_lp_mm->checkState() == Qt::Unchecked) {
        parent_lp_mm->setCheckState(Qt::Checked);
    }
}

void Settings_distortion_tab::stereographic_toggled() {
    if (rb_lens_stereo->isChecked() && parent_lp_mm->checkState() == Qt::Unchecked) {
        parent_lp_mm->setCheckState(Qt::Checked);
    }
}


