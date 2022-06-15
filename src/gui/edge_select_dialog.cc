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
#include "include/logger.h"
#include "edge_select_dialog.h"

Edge_select_dialog::Edge_select_dialog(QWidget* parent)
 : QDialog(parent, Qt::WindowMaximizeButtonHint | Qt::WindowCloseButtonHint), parent(parent) {
    
    QString cancel_tt(
        "Cancel the processing of the file currently shown.\n"
        "Skips on to loading next queued file, if any."
    );
    cancel_button = new QPushButton("Cancel");
    cancel_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    cancel_button->setToolTip(cancel_tt);
    
    QString accept_tt(
        "Accepts currently shown ROIs, and submits image and\n"
        "ROIs for further processing.\n"
        "Proceeds to loading of the next queued file, if any."
    );
    accept_button = new QPushButton("Accept");
    accept_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    accept_button->setToolTip(accept_tt);
    
    QString apply_all_tt(
        "Apply the current ROIs to all file(s) that are\n" 
        "currently queued for manual ROI processing."
    );
    apply_all_button = new QPushButton("Accept\n queued");
    apply_all_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    apply_all_button->setToolTip(apply_all_tt);
    
    QString abort_all_tt(
        "Cancel the processing of all files that are currently\n"
        "queued for manual ROI processing. This is equivalent to\n"
        "manually clicking 'Cancel' for each queued file."
    );
    abort_all_button = new QPushButton("Abort\n queued");
    abort_all_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    abort_all_button->setToolTip(abort_all_tt);

    help_button = new QPushButton("Help");
    help_button->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Minimum);
    
    QString load_tt("Load a previously saved ROI.");
    load_button = new QPushButton("Load");
    load_button->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    load_button->setToolTip(load_tt);
    
    QString save_tt("Save the current ROIs to a file.");
    save_button = new QPushButton("Save");
    save_button->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    save_button->setEnabled(false);
    save_button->setToolTip(save_tt);
    
    QString clear_roi_tt("Remove the current ROIs.");
    clear_roi_button = new QPushButton("Clear ROIs");
    clear_roi_button->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    clear_roi_button->setEnabled(false);
    clear_roi_button->setToolTip(clear_roi_tt);
    
    QString gamma_tt(
        "Toggle between linear (gamma = 1.0) and non-\n"
        "linear (gamma = 2.2) display styles.\n"
        "This does not affect the image data used\n"
        "in subsequent analyses."
    );
    gamma_switch = new QToolButton;
    QIcon gamma_icon;
    gamma_icon.addFile(":/Icons/Gamma_linear", QSize(), QIcon::Normal, QIcon::Off);
    gamma_icon.addFile(":/Icons/Gamma_two", QSize(), QIcon::Normal, QIcon::On);
    gamma_switch->setIcon(gamma_icon);
    gamma_switch->setToolButtonStyle(Qt::ToolButtonIconOnly);
    gamma_switch->setCheckable(true);
    gamma_switch->setFixedSize(66, 36);
    gamma_switch->setIconSize(QSize(64, 34));
    gamma_switch->setToolTip(gamma_tt);
    
    rb_autoload_roi = new QRadioButton("reuse ROIs");
    rb_blank_roi = new QRadioButton("clear ROIs");
    rb_blank_roi->setChecked(true);

    img_viewer = new GL_image_viewer(parent);
    img_viewer->set_clickable(true);
    img_viewer->set_resize_on_load(true);
    
    img_panel = new GL_image_panel_edges(img_viewer);
    img_panel->setMouseTracking(true);
    img_panel->set_max_scale_factor(4.0); // allow more zoom in this dialog
    
    img_viewer->setViewport(img_panel);   // TODO: could combine these
    img_viewer->set_GL_widget(img_panel);
    
    icon_image = std::shared_ptr<QImage>(new QImage(QSize(256, 256), QImage::Format_RGB888));
    icon_image->fill(Qt::red);
    img_panel->set_default_image(icon_image.get());

    help_dialog = new Manual_roi_help_dialog(this);
    help_dialog->setModal(true);
    
    img_filename = new QLabel("filename");
    img_progress = new QLabel("n remaining");
    text_img_filename = new QLabel("Name:");
    text_img_progress = new QLabel("Queued:");
    
    QString edge_length_tt(
        "Edge length of active ROI, in units\n"
        "of pixels at native resolution"
    );
    edge_length = new QLabel("");
    edge_length->setToolTip(edge_length_tt);
    text_edge_length = new QLabel("Edge length:");
    text_edge_length->setToolTip(edge_length_tt);
    
    QString img_max_val_tt(
        "Minimum and maximum pixel values (DNs)\n"
        "found in source image, across all\n"
        "channels (R,G,B) if applicable\n\n"
        "Min/Max do not represent split-edge\n"
        "histogram extremes on 16-bit images"
    );
    img_max_val = new QLabel;
    img_max_val->setToolTip(img_max_val_tt);
    text_img_max_val = new QLabel("DN range:");
    text_img_max_val->setToolTip(img_max_val_tt);
    
    QString type_tt(
        "File type, i.e., Gray, RGB, or Bayer channel, as\n"
        "well as bit depth.\n"
        "Note that the type reflects the choice in the\n"
        "Preferences at the time the File/Open action\n"
        "was performed, especially with regards to the\n"
        "Bayer channel setting. It is also indicative\n"
        "of the type of MTF that will be extracted."
    );
    img_type = new QLabel("type");
    img_type->setToolTip(type_tt);
    text_img_type = new QLabel("Type:");
    text_img_type->setToolTip(type_tt);
    
    QString histogram_tt(
        "The split-edge histogram displays two histograms, one for each\n"
        "half of the ROI, each half assumed to correspond to either a\n"
        "uniformly dark or uniformly light part of the image.\n" 
        "If these histograms overlap, the ROI selection is likely invalid,\n"
        "which is indicated by an amber background with the overlap-\n"
        "ping section displayed in red. Ideally, you want two well-\n"
        "separated histograms.\n\n"
        "Only displays when an ROI is selected or actively modified."
    );
    histogram = new Histo_widget(this);
    histogram->setToolTip(histogram_tt);
    
    QGridLayout* prop_layout = new QGridLayout;
    for (auto w: {text_img_filename, img_filename, text_img_type, text_img_progress, img_progress}) {
        w->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    }
    prop_layout->addWidget(text_img_filename, 0, 0, Qt::AlignTop);
    prop_layout->addWidget(img_filename, 0, 1, 1, 2, Qt::AlignTop);
    prop_layout->addWidget(text_img_type, 1, 0, Qt::AlignTop);
    prop_layout->addWidget(img_type, 1, 1, 1, 2, Qt::AlignTop);
    prop_layout->addWidget(text_img_progress, 2, 0, Qt::AlignTop);
    prop_layout->addWidget(img_progress, 2, 1, Qt::AlignTop);
    prop_layout->addWidget(text_img_max_val, 3, 0);
    prop_layout->addWidget(img_max_val, 3, 1);
    int marg = img_max_val->margin();
    img_max_val->setContentsMargins(marg, marg, 15, marg); // a little bit hacky, but we cannot set the margin of gamma_switch
    prop_layout->addWidget(gamma_switch, 2, 2, 2, 1, Qt::AlignRight);
    
    QGroupBox* prop_box = new QGroupBox("File properties");
    prop_box->setLayout(prop_layout);
    
    QGridLayout* histo_layout = new QGridLayout;
    histo_layout->addWidget(histogram, 0, 0, 1, 4);
    histo_layout->addWidget(text_edge_length, 1, 0);
    histo_layout->addWidget(edge_length, 1, 1, Qt::AlignLeft);
    
    QGroupBox* histo_box = new QGroupBox("Split-edge histogram");
    histo_box->setLayout(histo_layout);
    
    QGridLayout* button_layout = new QGridLayout;
    
    button_layout->addWidget(help_button, 0, 0);
    button_layout->addWidget(apply_all_button, 0, 1);
    button_layout->addWidget(abort_all_button, 0, 2);
    
    button_layout->addWidget(load_button, 2, 0);
    button_layout->addWidget(save_button, 2, 1);
    button_layout->addWidget(clear_roi_button, 2, 2);
    
    button_layout->addWidget(accept_button, 3, 1);
    button_layout->addWidget(cancel_button, 3, 2);
    
    QVBoxLayout* roi_layout = new QVBoxLayout;
    roi_layout->addWidget(rb_autoload_roi);
    roi_layout->addWidget(rb_blank_roi);
    roi_layout->addStretch();
    QGroupBox* roi_box = new QGroupBox("On next image");
    roi_box->setLayout(roi_layout);
    roi_box->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    
    QHBoxLayout* hlayout = new QHBoxLayout;
    hlayout->addWidget(histo_box);
    hlayout->addWidget(prop_box);
    hlayout->addStretch();
    hlayout->addWidget(roi_box, 0, Qt::AlignTop);
    hlayout->addLayout(button_layout);
    
    QGridLayout* vlayout = new QGridLayout;
    vlayout->addWidget(img_viewer, 0, 0);
    vlayout->addLayout(hlayout, 1, 0);
    
    connect(cancel_button, SIGNAL(clicked()), this, SLOT( close() ));
    connect(accept_button, SIGNAL(clicked()), this, SLOT( export_roi() ));
    connect(help_button, SIGNAL(clicked()), this, SLOT( show_help() ));
    connect(load_button, SIGNAL(clicked()), this, SLOT( load_roi() ));
    connect(save_button, SIGNAL(clicked()), this, SLOT( save_roi() ));
    connect(img_panel, SIGNAL( update_edge_length(double) ), this, SLOT( update_edge_length(double) ));
    connect(img_panel, SIGNAL( update_histogram(histo_t, histo_t) ), this, SLOT( update_histogram(histo_t, histo_t) ));
    connect(img_panel, SIGNAL( enable_save_button() ), this, SLOT( enable_save_button() ));
    connect(img_panel, SIGNAL( disable_save_button() ), this, SLOT( disable_save_button() ));
    connect(apply_all_button, SIGNAL( clicked() ), this, SLOT( apply_all_button_clicked() ));
    connect(abort_all_button, SIGNAL( clicked() ), this, SLOT( abort_all_button_clicked() ));
    connect(gamma_switch, SIGNAL( toggled(bool) ), this, SLOT( gamma_state(bool) ));
    connect(clear_roi_button, SIGNAL( clicked() ), this, SLOT( clear_roi_button_clicked() ));
    
    accept_button->setDefault(true);
    setModal(false);
    setWindowModality(Qt::NonModal);
    setLayout(vlayout);
    setWindowTitle("Select one or more edge ROIs");
}

bool Edge_select_dialog::load_image(QString img_name, 
    QStringList arguments, QString raw_filename) {
    
    bool load_success = img_viewer->load_image(img_name);
    
    if (!load_success) {
        return false;
    }
    
    // reset Bayer parameters (for histogram extraction later)
    bayer_channel = Bayer::bayer_t::NONE;
    cfa_pattern = Bayer::cfa_pattern_t::RGGB;
    cfa_mask = Bayer::to_cfa_mask(bayer_channel, cfa_pattern);
    
    QString type_string;
    if (get_panel()->image_channels() == 1) {
        if (arguments.contains("--bayer")) {
            type_string = "Bayer";
            extract_bayer_info(arguments);
            if (bayer_channel != Bayer::bayer_t::NONE) {
                type_string += " " + QString(Bayer::to_string(bayer_channel).c_str());
            }
        } else {
            type_string = "Gray";
        }
    } else {
        type_string = "RGB L";
    }
    if (get_panel()->image_depth() == 8) {
        type_string += ", 8-bit";
        gamma_switch->setChecked(false);
    } else {
        type_string += ", 16-bit";
        gamma_switch->setChecked(true);
    }
    gamma_state(gamma_switch->isChecked());
    img_type->setText(type_string);
    img_filename->setText(QFileInfo(raw_filename).fileName());
    img_max_val->setText(QString("[%0, %1]").arg(img_panel->get_img_min()).arg(img_panel->get_img_max()));
    
    img_panel->set_cfa_mask(cfa_mask);
    img_panel->clear_overlay();
    update();
    // force a repaint on the GL panel once the event loop runs
    QTimer::singleShot(0, [this]() { 
        img_panel->update();
        QCoreApplication::processEvents(); 
    });
    return true;
}

void Edge_select_dialog::open(void) {
    show();
}

void Edge_select_dialog::export_roi(void) {
    bool roi_success = img_panel->save_rois(roi_file);
    if (roi_success) {
        accept();
    } else {
        close();
    }
}

bool Edge_select_dialog::load_rois(QString fname) {
    return img_panel->load_rois(fname);
}

bool Edge_select_dialog::save_rois(QString fname) {
    return img_panel->save_rois(fname);
}

void Edge_select_dialog::load_roi(void) {
    QFileDialog dialog(this, "Select ROI file", QString(), QString());
    dialog.setFileMode(QFileDialog::ExistingFile);
    dialog.setNameFilter(tr("ROI files (*.roi)"));
    
    QStringList filenames;
    if (dialog.exec()) {
        filenames = dialog.selectedFiles();
        if (filenames.size() > 0 && filenames[0].length() > 0) {
            bool load_success = load_rois(filenames[0]);
            if (!load_success) {
                logger.error("Could not load ROI filename %s\n", 
                    filenames[0].toLocal8Bit().constData()
                );
            }
        }
    }
}

void Edge_select_dialog::save_roi(void) {
    QFileDialog dialog(this, "Save ROI file", QString(), QString("*.roi"));
    dialog.setFileMode(QFileDialog::AnyFile);
    dialog.setAcceptMode(QFileDialog::AcceptSave);
    dialog.setNameFilter(tr("ROI files (*.roi)"));
    
    QStringList filenames;
    if (dialog.exec()) {
        filenames = dialog.selectedFiles();
        if (filenames.size() > 0 && filenames[0].length() > 0) {
            if (QFileInfo(filenames[0]).suffix().length() == 0) {
                filenames[0] += ".roi";
            }
            save_rois(filenames[0]);
        }
    }
}

void Edge_select_dialog::show_help(void) {
    help_dialog->setModal(false);
    help_dialog->show();
}

void Edge_select_dialog::update_queue_size(int queue_size) {
    if (queue_size > 0) {
        img_progress->setText(QString("%1 files").arg(queue_size+1));
    } else {
        img_progress->setText(QString("%1 file").arg(queue_size+1));
    }
}

void Edge_select_dialog::update_edge_length(double length) {
    if (length > 0) {
        edge_length->setText(QString("%1").arg(length, 0, 'f', 1));
    } else {
        QTimer::singleShot(200, [this]() { 
            edge_length->setText(QString(""));
        });
        
    }
}

void Edge_select_dialog::update_histogram(histo_t dark, histo_t light) {
    histogram->set_histogram(dark, light);
}

void Edge_select_dialog::extract_bayer_info(const QStringList& arguments) {
    int bayer_idx = arguments.indexOf("--bayer");
    if (bayer_idx >= 0) {
        std::string bayer_channel_str = arguments[std::min(bayer_idx + 1, arguments.size())].toStdString();
        bayer_channel = Bayer::from_string(bayer_channel_str);
        cfa_pattern = Bayer::cfa_pattern_t::RGGB;
        int cfa_idx = arguments.indexOf("--cfa-pattern");
        if (cfa_idx >= 0) {
            std::string cfa_pattern_str = arguments[std::min(cfa_idx + 1, arguments.size())].toStdString();
            cfa_pattern = Bayer::from_cfa_string(cfa_pattern_str);
        }
        cfa_mask = Bayer::to_cfa_mask(bayer_channel, cfa_pattern);
    }
}

void Edge_select_dialog::disable_save_button(void) {
    save_button->setEnabled(false);
    clear_roi_button->setEnabled(false);
}

void Edge_select_dialog::enable_save_button(void) {
    save_button->setEnabled(true);
    clear_roi_button->setEnabled(true);
}

void Edge_select_dialog::apply_all_button_clicked(void) {
    static int manual_seq_number = 0;
    QString manual_roi_file = QString("%1/mtfmapper_batch_%2.roi").arg(QDir::tempPath()).arg(manual_seq_number++);
    bool save_success = img_panel->save_rois(manual_roi_file);
    if (save_success) {
        emit start_manual_queue_processing(manual_roi_file);
        export_roi();
    } else {
        logger.error("Could not save batch ROI filename %s\n", 
            manual_roi_file.toLocal8Bit().constData()
        );
    }
}

void Edge_select_dialog::abort_all_button_clicked(void) {
    emit manual_queue_skip_all();
    close();
}

void Edge_select_dialog::gamma_state(bool enabled) {
    if (enabled) {
        get_panel()->set_gamma(2.2);
    } else {
        get_panel()->set_gamma(1.0);
    }
    get_panel()->update();
}

void Edge_select_dialog::clear_roi_button_clicked(void) {
    img_panel->clear_overlay();
    img_panel->update();
}

