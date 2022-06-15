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

#include "settings_io_tab.h"

#include "common.h"
#include "nonempty_validator.h"

#include <QtGlobal>
#if QT_VERSION < QT_VERSION_CHECK(5,11,0)
#define horizontalAdvance width
#endif

const QString Settings_io_tab::setting_threshold = "setting_threshold_sauvola";
const QString Settings_io_tab::setting_threshold_default = "0.55";
const QString Settings_io_tab::setting_pixsize = "setting_pixelsize";
const QString Settings_io_tab::setting_pixsize_default = "4.78";
const QString Settings_io_tab::setting_bayer = "setting_bayer";
const QString Settings_io_tab::setting_esf_model = "setting_esf_model";
const QString Settings_io_tab::setting_linear_gamma = "setting_gamma";
const Qt::CheckState Settings_io_tab::setting_linear_gamma_default = Qt::Unchecked;
const QString Settings_io_tab::setting_annotation = "setting_annotation";
const Qt::CheckState Settings_io_tab::setting_annotation_default = Qt::Checked;
const QString Settings_io_tab::setting_profile = "setting_profile";
const Qt::CheckState Settings_io_tab::setting_profile_default = Qt::Checked;
const QString Settings_io_tab::setting_grid = "setting_grid";
const Qt::CheckState Settings_io_tab::setting_grid_default = Qt::Checked;
const QString Settings_io_tab::setting_lensprofile = "setting_lensprofile";
const Qt::CheckState Settings_io_tab::setting_lensprofile_default = Qt::Unchecked;
const QString Settings_io_tab::setting_orientation = "setting_orientation";
const Qt::CheckState Settings_io_tab::setting_orientation_default = Qt::Unchecked;
const QString Settings_io_tab::setting_autocrop = "setting_autocrop";
const Qt::CheckState Settings_io_tab::setting_autocrop_default = Qt::Unchecked;
const QString Settings_io_tab::setting_lpmm = "setting_lpmm";
const Qt::CheckState Settings_io_tab::setting_lpmm_default = Qt::Unchecked;
const QString Settings_io_tab::setting_gnuplot_scaled = "setting_gnuplot_scaled";
const Qt::CheckState Settings_io_tab::setting_gnuplot_scaled_default = Qt::Checked;
const QString Settings_io_tab::setting_zscale = "setting_zscale";
const int Settings_io_tab::setting_zscale_default = 0;
const QString Settings_io_tab::setting_cache = "image_cache_size";
const double Settings_io_tab::setting_cache_default = 1024;
const QString Settings_io_tab::setting_mtf_contrast = "mtf_contrast";
const double Settings_io_tab::setting_mtf_contrast_default = 50.0;
const QString Settings_io_tab::setting_lp1 = "lp1";
const QString Settings_io_tab::setting_lp1_default = "10";
const QString Settings_io_tab::setting_lp2 = "lp2";
const QString Settings_io_tab::setting_lp2_default = "30";
const QString Settings_io_tab::setting_lp3 = "lp3";
const QString Settings_io_tab::setting_lp3_default = "";
const QString Settings_io_tab::setting_arguments = "arguments";
const QString Settings_io_tab::setting_lensprofile_fixed = "lens_profile_fixed_scale";
const Qt::CheckState Settings_io_tab::setting_lensprofile_fixed_default = Qt::Unchecked;
const QString Settings_io_tab::setting_surface_max_value = "surface_max_value";
const QString Settings_io_tab::setting_surface_max_value_default = "0";
const QString Settings_io_tab::setting_surface_max_flag = "surface_max_flag";
const Qt::CheckState Settings_io_tab::setting_surface_max_flag_default = Qt::Unchecked;
const QString Settings_io_tab::setting_ca_active = "setting_ca";
const Qt::CheckState Settings_io_tab::setting_ca_active_default = Qt::Unchecked;
const QString Settings_io_tab::setting_ca_type = "setting_ca_type";
const int Settings_io_tab::setting_ca_type_default = 0;
const QString Settings_io_tab::setting_fullsfr = "setting_fullsfr";
const Qt::CheckState Settings_io_tab::setting_fullsfr_default = Qt::Unchecked;
const QString Settings_io_tab::setting_nosmoothing = "setting_nosmoothing";
const Qt::CheckState Settings_io_tab::setting_nosmoothing_default = Qt::Unchecked;


Settings_io_tab::Settings_io_tab(QWidget *parent ATTRIBUTE_UNUSED) {
    QSettings settings("mtfmapper", "mtfmapper");

    arguments_label = new QLabel(tr("Arguments:"), this);
    arguments_label->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    arguments_line = new QLineEdit(this);
    arguments_line->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);

    QFontMetrics fm(QApplication::font());
    int reasonable_width = fm.horizontalAdvance("1048576000");
    int adv_width = fm.horizontalAdvance("Threshold: ");
    int surf_max_reasonable_width = fm.horizontalAdvance("999.50");

    QString threshold_tt(
        "Relative intensity threshold used when automatically identifying dark targets on the\n"
        "lighter test chart background.\n"
        "If you find that MTF Mapper is not detecting all of your trapezoids then you should \n"
        "try to lower this threshold, maybe try 0.4, 0.3, 0.1. This is often necessary when your\n"
        "test chart does not fill the image, or if you have a broad dark border around the test\n"
        "chart, or if the contrast between the dark targets and the chart background is low.\n\n"
        "The default value of 0.55 was chosen to minimise the number of potential objects that\n"
        "have to be tested to see if they are trapezoids; a lower value might help to find\n"
        "missed trapezoids, but it could slow processing down because more objects have to be\n"
        "screened. So leaving this value at 0.1 is not recommended for general use."
    );
    QDoubleValidator* dv_thresh = new Nonempty_DoubleValidator(0.001, 1.0, 3, 0.55, this);
    threshold_label = new QLabel(tr("Threshold:"), this);
    threshold_label->setMinimumWidth(adv_width);
    threshold_label->setToolTip(threshold_tt);
    threshold_line = new QLineEdit(this);
    threshold_line->setValidator(dv_thresh);
    threshold_line->setMaximumWidth(reasonable_width);
    threshold_line->setToolTip(threshold_tt);

    QString pixsize_tt(
        "Pixel pitch (or size), in microns. Note that this value is only used when the \n"
        "line pairs/mm units flag (above) is enabled. \n"
        "Also note that this value is baked into the output at the time that you \n"
        "process an image using the File/Open menu, but you can change it between \n"
        "processing images from different cameras."
    );
    QDoubleValidator* dv_pixsize = new Nonempty_DoubleValidator(0.001, 999.0, 3, 4.0, this);
    pixsize_label = new QLabel(tr("Pixel size:"), this);
    pixsize_label->setMinimumWidth(adv_width);
    pixsize_label->setToolTip(pixsize_tt);
    pixsize_line = new QLineEdit(this);
    pixsize_line->setValidator(dv_pixsize);
    pixsize_line->setMaximumWidth(reasonable_width);
    pixsize_line->setToolTip(pixsize_tt);
    pixsize_units = new QLabel("\xc2\xb5m", this);

    normal_pal = pixsize_label->palette();
    disabled_pal = pixsize_label->palette();
    disabled_pal.setCurrentColorGroup(QPalette::Disabled);
    disabled_pal.setColorGroup(
        QPalette::Normal, 
        disabled_pal.windowText(), disabled_pal.button(), disabled_pal.light(), disabled_pal.dark(),
        disabled_pal.mid(), disabled_pal.text(), disabled_pal.brightText(), disabled_pal.base(),
        disabled_pal.window()
    );

    QDoubleValidator* dv_cachesize = new Nonempty_DoubleValidator(1, 1024 * 1024, 1, 1024, this);
    cache_label = new QLabel(tr("Cache size:"), this);
    cache_label->setMinimumWidth(adv_width);
    cache_line = new QLineEdit(this);
    cache_line->setValidator(dv_cachesize);
    cache_line->setMaximumWidth(reasonable_width);

    QIntValidator* dv_contrast = new Nonempty_IntValidator(1, 90, 50, this);
    contrast_label = new QLabel(tr("MTF-XX:"), this);
    contrast_label->setMinimumWidth(adv_width);
    contrast_line = new QLineEdit(this);
    contrast_line->setValidator(dv_contrast);
    contrast_line->setMaximumWidth(reasonable_width);


    QIntValidator* dv_lpmm = new QIntValidator(1, 500, this);
    lpmm_label = new QLabel(tr("Lensprofile lp/mm:"), this);
    lp1_line = new QLineEdit(this);
    lp1_line->setValidator(dv_lpmm);
    lp1_line->setMaximumWidth(fm.horizontalAdvance("50000"));
    lp2_line = new QLineEdit(this);
    lp2_line->setValidator(dv_lpmm);
    lp2_line->setMaximumWidth(fm.horizontalAdvance("50000"));
    lp3_line = new QLineEdit(this);
    lp3_line->setValidator(dv_lpmm);
    lp3_line->setMaximumWidth(fm.horizontalAdvance("50000"));

    cb_linear_gamma = new QCheckBox("Linear gamma (8 bit)", this);
    cb_linear_gamma->setToolTip(
        "Force MTF Mapper to treat the input image as if it has a linear intensity encoding.\n"
        "MTF Mapper will use the ICC profile TRC to linearize the input image if one is present.\n"
        "For example, untagged 8-bit images (meaning no ICC profile) are assumed to be encoded in\n"
        "the sRGB TRC, roughly gamma=2.2, and linearized accordingly. If you know that your 8-bit\n"
        "image is in fact linear, such as when using a raw image from a machine vision camera, \n"
        "then you should enable this option.\n\n"
        "In contrast, untagged 16-bit images are assumed to have a linear intensity encoding.\n\n"
        "Enabling this option will override any ICC profile, regardless of whether the image is \n"
        "8-bit ot 16-bit."
    );

    cb_annotation = new QCheckBox("Annotation", this);
    cb_annotation->setToolTip(
        "Produce an output image that is annotated with the MTF50 value of each edge.\n"
        "Once you are viewing the annotated image in the MTF Mapper GUI, you can left-click\n"
        "on the text annotations (usually cyan in colour) to view the SFR curve, or the ESF \n"
        "and LSF plots.\n"
        "You can <shift>-left-click on up to three edges to compare their SFR/ESF/LSF plots.\n"
        "Note that when using <shift>-left-click to compare edges, you may select the edges\n"
        "from any of the annotated images you have produced in the current session.\n\n"
        "This output type should work with practically any test chart containing isolated \n"
        "black trapezoidal shapes on a white background.\n"
        "It also works with the File/Open single edge image method of opening images."
    );

    cb_profile = new QCheckBox("Profile", this);
    cb_grid = new QCheckBox("Grid", this);
    cb_lensprofile = new QCheckBox("Lens profile", this);
    cb_orientation = new QCheckBox("Chart orientation", this);

    cb_autocrop = new QCheckBox("Autocrop", this);
    cb_autocrop->setToolTip(
        "Try to automatically remove any black background around the test chart.\n"
        "Mostly useful when the image circle is much smaller than sensor.\n"
        "It may even work when the chart is framed to fill only part of the \n"
        "image, but no guarantees."
    );
    cb_lpmm = new QCheckBox("Line pairs/mm units", this);
    cb_lpmm->setToolTip(
        "When enabled, the output units used in the following are affected: \n"
        "a) MTF50 values are expressed in lp/mm (cycles/pixel, or c/p otherwise)\n"
        "b) SFR frequencies in SFR plot are expressed in lp/mm (c/p otherwise)\n"
        "c) Lens profile frequencies are interpreted as lp/mm (c/p otherwise)\n"
        "d) CA shifts are measured in microns (pixels, otherwise)\n"
        "\n*Note that the lp/mm unit choice is baked into the outputs at the time when \n"
        "you processed the file via the File/Open menu, except for the SFR plot which is dynamic."
    );

    cb_gnuplot_scaled = new QCheckBox("Scale plots to window", this);
    cb_gnuplot_scaled->setToolTip("Output plot images are scaled to match the main window size.\nUseful on very high resolution displays.");

    cb_lensprofile_fixed = new QCheckBox("Lens profile fixed scale", this);
    cb_lensprofile_fixed->setToolTip(
        "The x-axis scale of the lens profile plot is fixed to the image size.\n"
        "Otherwise the x-axis will scale according to the fraction of the image filled by the chart."
    );

    cb_surface_max = new QCheckBox("3D plot z-axis max value", this);
    cb_surface_max->setToolTip(
        "Specify a fixed upper limit to the z-axis (height) of the 3D surface \n"
        "plots (generated by 'Grid' output option).\n"
        "Useful when generating a stack of plots at different lens apertures, for example."
    );

    cb_ca_active = new QCheckBox("Chromatic Aberration", this);
    cb_ca_active->setToolTip(
        "Measure lateral Chromatic Aberration. This requires a colour input image to make sense.\n"
        "Works with both raw Bayer mosaiced images (but you have to select a 'Bayer channel' other \n"
        "than 'none' in the 'ESF construction' panel), or you can calculate lateral CA from RGB images (RGB per pixel).\n"
        "If you provide a raw camera image and leave the 'Bayer channel' setting to 'none', then dcraw will be used to \n"
        "demosaic the image first, which means you are including the demosaicing approximations in your CA measurement."
    );

    cb_fullsfr = new QCheckBox("Extended SFR domain", this);
    cb_fullsfr->setToolTip(
        "Computes SFR up to 2.0 cycles/pixel, rather than default of 1 cycle/pixel.\n"
        "Effect is only really visible in SFR curve plotting (and only after re-\n"
        "processing an image via File->Open).\n"
    );

    cb_nosmoothing = new QCheckBox("Reduced SFR smoothing", this);
    cb_nosmoothing->setToolTip(
        "Disables Savitzky-Golay SFR smoothing.\n"
    );

    box_colour = new QComboBox;
    box_colour->addItem("none");
    box_colour->addItem("red");
    box_colour->addItem("green");
    box_colour->addItem("blue");

    box_esf_model = new QComboBox;
    box_esf_model->addItem("kernel");
    box_esf_model->addItem("loess");

    box_ca_type = new QComboBox;
    box_ca_type->addItem("direct (pixels)");
    box_ca_type->addItem("radial distance %");

    QDoubleValidator* dv_sm = new Nonempty_DoubleValidator(0, 999.0, 2, 80.0, this);
    surface_max_value = new QLineEdit(this);
    surface_max_value->setMaxLength(5);
    surface_max_value->setText(settings.value(setting_surface_max_value, setting_surface_max_value_default).toString());
    surface_max_value->setValidator(dv_sm);
    surface_max_value->setMaximumWidth(surf_max_reasonable_width);

    threshold_line->setText(settings.value(setting_threshold, setting_threshold_default).toString());
    pixsize_line->setText(settings.value(setting_pixsize, setting_pixsize_default).toString());
    contrast_line->setText(settings.value(setting_mtf_contrast, setting_mtf_contrast_default).toString());
    cache_line->setText(settings.value(setting_cache, setting_cache_default).toString());
    lp1_line->setText(settings.value(setting_lp1, setting_lp1_default).toString());
    lp2_line->setText(settings.value(setting_lp2, setting_lp2_default).toString());
    lp3_line->setText(settings.value(setting_lp3, setting_lp3_default).toString());
    arguments_line->setText(settings.value(setting_arguments, "").toString());
    cb_linear_gamma->setCheckState(
        (Qt::CheckState)settings.value(setting_linear_gamma, setting_linear_gamma_default).toInt()
    );
    cb_annotation->setCheckState(
        (Qt::CheckState)settings.value(setting_annotation, setting_annotation_default).toInt()
    );
    cb_profile->setCheckState(
        (Qt::CheckState)settings.value(setting_profile, setting_profile_default).toInt()
    );
    cb_grid->setCheckState(
        (Qt::CheckState)settings.value(setting_grid, setting_grid_default).toInt()
    );
    cb_lensprofile->setCheckState(
        (Qt::CheckState)settings.value(setting_lensprofile, setting_lensprofile_default).toInt()
    );
    cb_orientation->setCheckState(
        (Qt::CheckState)settings.value(setting_orientation, setting_orientation_default).toInt()
    );
    cb_autocrop->setCheckState(
        (Qt::CheckState)settings.value(setting_autocrop, setting_autocrop_default).toInt()
    );
    cb_lpmm->setCheckState(
        (Qt::CheckState)settings.value(setting_lpmm, setting_lpmm_default).toInt()
    );
    cb_gnuplot_scaled->setCheckState(
        (Qt::CheckState)settings.value(setting_gnuplot_scaled, setting_gnuplot_scaled_default).toInt()
    );
    cb_lensprofile_fixed->setCheckState(
        (Qt::CheckState)settings.value(setting_lensprofile_fixed, setting_lensprofile_fixed_default).toInt()
    );
    cb_surface_max->setCheckState(
        (Qt::CheckState)settings.value(setting_surface_max_flag, setting_surface_max_flag_default).toInt()
    );
    cb_ca_active->setCheckState(
        (Qt::CheckState)settings.value(setting_ca_active, setting_ca_active_default).toInt()
    );
    cb_fullsfr->setCheckState(
        (Qt::CheckState)settings.value(setting_fullsfr, setting_fullsfr_default).toInt()
    );
    cb_nosmoothing->setCheckState(
        (Qt::CheckState)settings.value(setting_nosmoothing, setting_nosmoothing_default).toInt()
    );

    box_colour->setCurrentIndex(settings.value(setting_bayer, 0).toInt());
    box_esf_model->setCurrentIndex(settings.value(setting_esf_model, 0).toInt());
    box_ca_type->setCurrentIndex(settings.value(setting_ca_type, setting_ca_type_default).toInt());


    surface_max_value->setText(settings.value(setting_surface_max_value, setting_surface_max_value_default).toString());

    zscale_label = new QLabel("3D plot z-axis relative scale factor", this);
    zscale_slider = new QSlider(Qt::Horizontal, this);
    zscale_slider->setFocusPolicy(Qt::StrongFocus);
    zscale_slider->setTickPosition(QSlider::TicksAbove);
    zscale_slider->setMinimum(0);
    zscale_slider->setMaximum(20);
    zscale_slider->setTickInterval(5);
    zscale_slider->setSingleStep(1);
    zscale_slider->setValue(settings.value(setting_zscale, setting_zscale_default).toInt());


    QGroupBox* voGroupBox = new QGroupBox(tr("Output types"), this);
    QVBoxLayout* vo_layout = new QVBoxLayout;
    vo_layout->addWidget(cb_annotation);
    vo_layout->addWidget(cb_profile);
    vo_layout->addWidget(cb_grid);
    vo_layout->addWidget(cb_lensprofile);
    vo_layout->addWidget(cb_orientation);
    vo_layout->addWidget(cb_ca_active);
    voGroupBox->setLayout(vo_layout);

    QGroupBox* inflags_gb = new QGroupBox(tr("Input flags"), this);
    QVBoxLayout* cb_layout = new QVBoxLayout;
    cb_layout->addWidget(cb_linear_gamma);
    cb_layout->addWidget(cb_lpmm);

    QHBoxLayout* pixsize_layout = new QHBoxLayout;
    pixsize_layout->addWidget(pixsize_label);
    pixsize_layout->addWidget(pixsize_line);
    pixsize_layout->addWidget(pixsize_units);
    pixsize_layout->addStretch(2);
    cb_layout->addLayout(pixsize_layout);

    QHBoxLayout* thresh_layout = new QHBoxLayout;
    thresh_layout->addWidget(threshold_label);
    thresh_layout->addWidget(threshold_line);
    thresh_layout->addStretch(2);
    cb_layout->addLayout(thresh_layout);

    cb_layout->addWidget(cb_autocrop);

    inflags_gb->setLayout(cb_layout);

    QGroupBox* bayer_GroupBox = new QGroupBox(tr("ESF construction"), this);
    QGridLayout* rb_layout = new QGridLayout;
    bayer_label = new QLabel("Bayer channel:", this);
    rb_layout->addWidget(bayer_label, 0, 0);
    rb_layout->addWidget(box_colour, 0, 1);
    esf_model_label = new QLabel("ESF model:", this);
    rb_layout->addWidget(esf_model_label, 1, 0);
    rb_layout->addWidget(box_esf_model, 1, 1);
    bayer_GroupBox->setLayout(rb_layout);

    QGroupBox* ca_GroupBox = new QGroupBox(tr("Chromatic aberration"), this);
    QGridLayout* ca_layout = new QGridLayout;
    ca_type_label = new QLabel("CA display type:", this);
    ca_layout->addWidget(ca_type_label, 0, 0);
    ca_layout->addWidget(box_ca_type, 0, 1);
    ca_GroupBox->setLayout(ca_layout);

    QGroupBox* advanced = new QGroupBox(tr("Advanced"), this);
    QVBoxLayout* adv_layout = new QVBoxLayout;

    QHBoxLayout* r3_layout = new QHBoxLayout;
    r3_layout->addWidget(contrast_label);
    r3_layout->addWidget(contrast_line);
    r3_layout->addWidget(new QLabel("% contrast", this));
    r3_layout->addStretch(2);
    adv_layout->addLayout(r3_layout);

    QHBoxLayout* r4_layout = new QHBoxLayout;
    r4_layout->addWidget(lpmm_label);
    r4_layout->addWidget(lp1_line);
    r4_layout->addWidget(lp2_line);
    r4_layout->addWidget(lp3_line);
    r4_layout->addStretch(2);
    adv_layout->addLayout(r4_layout);

    QHBoxLayout* r5_layout = new QHBoxLayout;
    r5_layout->addWidget(cb_lensprofile_fixed);
    adv_layout->addLayout(r5_layout);

    QHBoxLayout* r5_5_layout = new QHBoxLayout;
    r5_5_layout->addWidget(cb_gnuplot_scaled);
    adv_layout->addLayout(r5_5_layout);

    QHBoxLayout* r6_layout = new QHBoxLayout;
    r6_layout->addWidget(zscale_label);
    r6_layout->addStretch(2);
    adv_layout->addLayout(r6_layout);

    QHBoxLayout* r7_layout = new QHBoxLayout;
    r7_layout->addWidget(zscale_slider);
    adv_layout->addLayout(r7_layout);

    surface_max_units = new QLabel(cb_lpmm->checkState() == Qt::Unchecked ? "c/p" : "lp/mm", this);
    surface_max_units->setFixedWidth(fm.horizontalAdvance("lp/mm"));
    QHBoxLayout* r8_layout = new QHBoxLayout;
    r8_layout->addWidget(cb_surface_max);
    r8_layout->addWidget(surface_max_value);
    r8_layout->addWidget(surface_max_units);
    r8_layout->addStretch(2);
    adv_layout->addLayout(r8_layout);

    QHBoxLayout* r9_layout = new QHBoxLayout;
    r9_layout->addWidget(cache_label);
    r9_layout->addWidget(cache_line);
    r9_layout->addWidget(new QLabel("MB", this));
    r9_layout->addStretch(2);
    adv_layout->addLayout(r9_layout);

    adv_layout->addWidget(cb_fullsfr);
    adv_layout->addWidget(cb_nosmoothing);

    QHBoxLayout* r10_layout = new QHBoxLayout;
    r10_layout->addWidget(arguments_label);
    r10_layout->addWidget(arguments_line);
    adv_layout->addLayout(r10_layout);

    advanced->setLayout(adv_layout);

    QGridLayout* vlayout = new QGridLayout;
    vlayout->addWidget(voGroupBox, 0, 0, 1, 2);
    vlayout->addWidget(inflags_gb, 1, 0, 1, 2);
    vlayout->addWidget(bayer_GroupBox, 2, 0, 1, 2);
    vlayout->addWidget(ca_GroupBox, 2, 3, 1, 2);
    vlayout->addWidget(advanced, 0, 3, 2, 2);

    connect(cb_lpmm, SIGNAL(clicked()), this, SLOT(lpmm_toggled()));
    lpmm_toggled(); // to ensure the correct enabled/disabled appearance is used

    setLayout(vlayout);
}

void Settings_io_tab::lpmm_toggled(void) {
    if (cb_lpmm->checkState() == Qt::Unchecked) {
        pixsize_label->setPalette(disabled_pal);
        pixsize_line->setPalette(disabled_pal);
        pixsize_units->setPalette(disabled_pal);
    } else {
        pixsize_label->setPalette(normal_pal);
        pixsize_line->setPalette(normal_pal);
        pixsize_units->setPalette(normal_pal);
    }
}

