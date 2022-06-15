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
#include "include/logger.h"
#include <QtWidgets>
#include "mtfmapper_app.h"

#include "worker_thread.h"
#include "common.h"

#include <string>
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;
using std::string; 
 
mtfmapper_app::mtfmapper_app(QWidget *parent ATTRIBUTE_UNUSED)
  : processor(this), sfr_dialog(nullptr) 
{
    img_comment_label = new QLabel("Comment: ");
    img_comment_label->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    img_comment_value = new QLabel("N/A");
    img_comment_value->setAlignment(Qt::AlignLeft);
    img_comment_value->setStyleSheet("QLabel { color : SteelBlue; }");
    
    af_ft_label = new QLabel("AF Tune Value: ");
    af_ft_label->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    af_ft_value = new QLabel("N/A");
    af_ft_value->setAlignment(Qt::AlignLeft);
    af_ft_value->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    af_ft_value->setStyleSheet("QLabel { color : SteelBlue; }");

    focal_length_label = new QLabel("Focal Length: ");
    focal_length_label->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    focal_length_value = new QLabel("N/A");
    focal_length_value->setAlignment(Qt::AlignLeft);
    focal_length_value->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    focal_length_value->setStyleSheet("QLabel { color : SteelBlue; }");

    focus_distance_label = new QLabel("Focus Distance: ");
    focus_distance_label->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    focus_distance_value = new QLabel("N/A");
    focus_distance_value->setAlignment(Qt::AlignLeft);
    focus_distance_value->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    focus_distance_value->setStyleSheet("QLabel { color : SteelBlue; }");

    img_frame = new QFrame(this);
    
    tb_img_annotated = new QCheckBox("Annotated image");
    tb_img_annotated->setChecked(true);
    tb_img_profile = new QCheckBox("Profile");
    tb_img_profile->setChecked(true);
    tb_img_gridimg = new QCheckBox("Grid");
    tb_img_gridimg->setChecked(true);
    tb_img_lensprofile = new QCheckBox("Lens profile");
    tb_img_lensprofile->setChecked(true);
    tb_img_orientation = new QCheckBox("Chart orientation");
    tb_img_orientation->setChecked(true);
    tb_img_ca = new QCheckBox("Chromatic aberration");
    tb_img_ca->setChecked(true);
    
    /*
    // TODO: set this up for now, later add it to File/Open dialog
    pixel_pitch = new QComboBox;
    pixel_pitch->setEditable(true);
    pixel_pitch->setPlaceholderText("pitch");
    pixel_pitch->setDuplicatesEnabled(false);
    prop_layout->addWidget(pixel_pitch, 3, 0);
    */
    
    img_viewer = new GL_image_viewer(this);
    img_panel = new GL_image_panel_dots(img_viewer, this);
    img_viewer->setViewport(img_panel);   // TODO: could combine these
    img_viewer->set_GL_widget(img_panel);
    
    
    QStringList labels;
    labels.push_back(QString("Data set"));
    dataset_contents.setHorizontalHeaderLabels(labels);
    
    datasets = new QTreeView;
    datasets->resize(300,400);
    datasets->move(610,0);
    datasets->setModel(&dataset_contents);
    datasets->setRootIsDecorated(true);
    
    progress = new QProgressBar;

    clear_button = new QPushButton("Clear results");
    clear_button->setEnabled(false);

    save_button = new QPushButton("Save all results");
    save_button->setEnabled(false);

    save_subset_button = new QPushButton("Save result subset");
    save_subset_button->setEnabled(false);

    
    QGridLayout* tb_layout = new QGridLayout;
    tb_layout->addWidget(datasets, 0, 0);
    tb_layout->addWidget(clear_button, 1, 0);
    tb_layout->addWidget(save_button, 2, 0);
    tb_layout->addWidget(save_subset_button, 3, 0);
    QGroupBox* vbox2 = new QGroupBox(tr("Data set control"));
    vbox2->setLayout(tb_layout);
    vbox2->setMinimumWidth(200);
    
    QGroupBox* v3GroupBox = new QGroupBox(tr("Image properties"));
    QGridLayout* hlayout = new QGridLayout;
    hlayout->addWidget(img_comment_label, 0, 0);
    hlayout->addWidget(img_comment_value, 0, 1);
    hlayout->addWidget(af_ft_label, 0, 2);
    hlayout->addWidget(af_ft_value, 0, 3);
    hlayout->addWidget(focal_length_label, 1, 0);
    hlayout->addWidget(focal_length_value, 1, 1);
    hlayout->addWidget(focus_distance_label, 1, 2);
    hlayout->addWidget(focus_distance_value, 1, 3);
    v3GroupBox->setLayout(hlayout);

    
    QGroupBox* vGroupBox = new QGroupBox(tr("Current output"));
    QGridLayout* vlayout = new QGridLayout;
    vlayout->addWidget(img_viewer, 0, 0);
    vlayout->addWidget(vbox2, 0, 1);
    vlayout->addWidget(v3GroupBox);
    vGroupBox->setLayout(vlayout);
    
    splitter = new QSplitter(Qt::Horizontal);
    splitter->addWidget(vGroupBox);
    splitter->addWidget(vbox2);
    splitter->setStretchFactor(0, 1);
    splitter->setCollapsible(0, false);
    splitter->setCollapsible(1, false);

    
    abort_button = new QPushButton("Abort");
    abort_button->hide();
    
    QGridLayout *mainLayout = new QGridLayout;
    mainLayout->addWidget(splitter, 0, 0, 1, 2);
    mainLayout->addWidget(progress, 1, 0);
    mainLayout->addWidget(abort_button, 1, 1);
    img_frame->setLayout(mainLayout);
    
    setCentralWidget(img_frame);
    
    worker_thread = new QThread(this);
    
    settings = new Settings_dialog(this);
    about    = new About_dialog(this);
    help     = new Help_dialog(this);
    
    create_actions();
    
    file_menu = new QMenu(tr("&File"), this);
    file_menu->addAction(open_act);
    file_menu->addAction(open_manual_roi_act);
    file_menu->addAction(open_roi_act);
    file_menu->addAction(open_focus_act);
    file_menu->addAction(open_imatest_act);
    file_menu->addSeparator();
    file_menu->addAction(exit_act);
    
    settings_menu = new QMenu(tr("&Settings"), this);
    settings_menu->addAction(prefs_act);

    help_menu = new QMenu(tr("&Help"), this);
    help_menu->addAction(help_act);
    help_menu->addAction(about_act);
    
    menuBar()->addMenu(file_menu);
    menuBar()->addMenu(settings_menu);
    menuBar()->addMenu(help_menu);
    
    connect(datasets->selectionModel(), SIGNAL(currentChanged(const QModelIndex&, const QModelIndex&)), this, SLOT(dataset_selected_changed(const QModelIndex&, const QModelIndex&)));
    
    connect(&processor, SIGNAL(send_parent_item(QString, QString, QString)), this, SLOT(parent_item(QString, QString, QString)));
    connect(&processor, SIGNAL(send_child_item(QString, QString)), this, SLOT(child_item(QString, QString)));
    connect(&processor, SIGNAL(send_close_item()), this, SLOT(close_item()));
    connect(&processor, SIGNAL(send_delete_item(QString)), this, SLOT(item_for_deletion(QString)));
    connect(&processor, SIGNAL(send_exif_filename(QString, QString)), this, SLOT(populate_exif_info_from_file(QString, QString)));
    
    connect(&processor, SIGNAL(send_progress_indicator(int)), this, SLOT(update_progress(int)));
    connect(&processor, SIGNAL(send_progress_indicator_max(int, int)), progress, SLOT(setRange(int, int)));
    connect(this, SIGNAL(submit_batch(Processor_state)), &processor, SLOT(receive_batch(Processor_state)));
    connect(settings, SIGNAL(set_cache_size(int)), this, SLOT(set_cache_size(int)));
    connect(settings, SIGNAL(settings_saved()), this, SLOT(settings_saved()));

    connect(&processor, SIGNAL(send_all_done()), this, SLOT(processor_completed()));
    connect(abort_button, SIGNAL(clicked()), &processor, SLOT(receive_abort()));
    
    connect(&processor, SIGNAL(send_all_done()), this, SLOT(enable_clear_button()));
    connect(clear_button, SIGNAL(clicked()), this, SLOT(clear_button_pressed()));
    connect(clear_button, SIGNAL(clicked()), this, SLOT(disable_save_button()));

    connect(&processor, SIGNAL(send_all_done()), this, SLOT(enable_save_button()));
    connect(save_button, SIGNAL(clicked()), this, SLOT(save_button_pressed()));

    connect(save_subset_button, SIGNAL(clicked()), this, SLOT(save_subset_button_pressed()));

    connect(&processor, SIGNAL(send_all_done()), this, SLOT(enable_file_open()));
    connect(&processor, SIGNAL(mtfmapper_call_failed(Worker_thread::failure_t, const QString&)), this, SLOT(mtfmapper_call_failed(Worker_thread::failure_t, const QString&)));
    connect(&processor, SIGNAL(send_processing_command(const Processing_command&)), this, SLOT(add_manual_roi_file(const Processing_command&)));

    zoom_scroll_tt = QString(
        "Use <ctrl>+scroll-wheel to zoom in/out, or\n"
        "Use right-click & drag mouse up/down to zoom in/out (touchpad-friendly), or\n"
        "Use <plus>/<minus> ('+'/'-') keys to zoom in/out; note that you may have to left-click \n"
        "in the 'Current output' window once before using the +/- keys.\n\n"
        "Use left-click & drag to pan around the image when zoomed in, or\n"
        "Use the mouse wheel to scroll vertically, or <shift>+scroll-wheel to scroll horizontally.\n"
    );

    img_viewer->setToolTip(zoom_scroll_tt);

    annotated_tt = QString(
        "Plot the SFR/ESF/LSF curves of an edge by left-clicking on the numbers over an edge (usually\n"
        "cyan in colour, but they could be yellow or red). Note that you will only see these numbers on\n"
        "edges that were processed successfully.\n\n"
        "You can plot the curves of up to three edges concurrently by holding down <shift> while\n"
        "left-clicking on the cyan numbers; a regular left-click without <shift> will revert to plotting\n"
        "only a single edge.\n\n"
        "Note that you can left-click followed by <shift>-left-click on edges belonging to different\n"
        "annotated output images, allowing you to compare the curves across different settings / cameras."
    );

    // rather a lot of code, but extract the icon and store it in a QImage
    mtfmapper_logo = new QIcon;
    mtfmapper_logo->addFile(":/Icons/AppIcon256");
    
    QGraphicsScene qgs;
    QGraphicsPixmapItem qgpi;
    qgpi.setTransformationMode(Qt::SmoothTransformation);
    qgs.addItem(&qgpi);
    qgpi.setPixmap(mtfmapper_logo->pixmap(256));
    qgs.setSceneRect(QRectF(0, 0, 255, 255));
    icon_image = new QImage(QSize(256, 256), QImage::Format_RGB888);
    icon_image->fill(Qt::white);
    QPainter painter(icon_image);
    qgs.render(&painter);
    img_panel->set_default_image(icon_image);
    
    setWindowTitle(tr("MTF Mapper"));
    resize(920,600);
    setAcceptDrops(true);
    
    edge_select_dialog = new Edge_select_dialog(this);
    connect(this, SIGNAL( manual_roi_queue_size(int) ), edge_select_dialog, SLOT( update_queue_size(int) ));
    connect(edge_select_dialog, SIGNAL( start_manual_queue_processing(QString) ), this, SLOT( enable_manual_queue_processing(QString) ));
    connect(edge_select_dialog, SIGNAL( manual_queue_skip_all() ), this, SLOT( enable_manual_queue_skip() ));
    connect(abort_button, SIGNAL( clicked() ), edge_select_dialog, SLOT( abort_all_button_clicked() ));
    
    sfr_dialog = new Sfr_dialog(this, sfr_dialog_geom);
    connect(sfr_dialog, SIGNAL(sfr_dialog_closed()), img_viewer, SLOT(clear_overlay()));
    connect(sfr_dialog, SIGNAL(send_geometry(QRect)), this, SLOT(sfr_dialog_closed(QRect)));
    connect(settings->accept_button, SIGNAL(clicked()), sfr_dialog, SLOT(update_lp_mm_mode()));
    
    check_if_helpers_exist();
    check_and_purge_stale_temp_files();
    
    // start up processor in a different thread
    QThread* thread = new QThread();
    processor.moveToThread(thread);
    connect(&processor, SIGNAL(work_finished()), thread, SLOT(quit()));
    connect(thread, SIGNAL(finished()), thread, SLOT(deleteLater()));
    thread->start();
}

mtfmapper_app::~mtfmapper_app(void) {
    // first, attempt to just remove the files we created, and know of, since
    // this is safer than deleting the temp directories recursively
    clear_temp_files();

    // now check if anything was left over in the temp directories, and
    // offer to remove them, which may potentially remove user-created files
    // but we do warn them first ...
    check_and_purge_stale_temp_files();

    delete icon_image;
    delete mtfmapper_logo;
}

void mtfmapper_app::clear_temp_files(void) {
    QStringList dirnames;
    for (int i=0; i < tempfiles_to_delete.size(); i++) {
        QString fn(tempfiles_to_delete.at(i));
        QString dn(QFileInfo(fn).absolutePath());
        
        QFile().remove(fn);
        QDir().rmdir(dn);
            
    }
    clear_button->setEnabled(false);
}

void mtfmapper_app::check_and_purge_stale_temp_files(void) {
    // check for, and list, any remaining non-empty temp dirs
    QDir tempdir = QDir::temp();
    tempdir.setNameFilters(QStringList() << "mtfmappertemp_*");
    QStringList dirlist = tempdir.entryList(QDir::Dirs);
    if (dirlist.size() > 0) {
        QListWidget* subset_list = new QListWidget();
        for (int i = 0; i < dirlist.size(); i++) {
            subset_list->addItem(dirlist[i]);
        }
        QDialog* subset_mb = new QDialog(this);
        subset_mb->setWindowTitle("Delete stale temp directories");
        QGridLayout* od_gridbox = new QGridLayout();
        QPlainTextEdit* description = new QPlainTextEdit;
        description->appendPlainText("A previous run of MTF Mapper appears to have left behind some stale temporary files and directories.\n");
        description->appendPlainText(QString("The directories are located under\n%1\n\nThe following directories and all their contents, will be deleted:").arg(QDir::toNativeSeparators(tempdir.absolutePath())));
        description->setReadOnly(true);
        description->setFrameStyle(QFrame::NoFrame);
        QColor current_bg_colour = palette().color(QPalette::Window);
        description->setStyleSheet(QString("QPlainTextEdit[readOnly=\"true\"] { background-color: %0 }").arg(current_bg_colour.name(QColor::HexRgb)));
        QFontMetrics d_fm(description->font());
        int row_height = d_fm.lineSpacing();
        description->setFixedHeight(row_height * 11);
        od_gridbox->addWidget(description, 0, 0, 1, 2);
        int rowheight = subset_list->model()->rowCount() * subset_list->sizeHintForRow(0);
        rowheight = std::min(std::max(rowheight, row_height * (dirlist.size() + 1)), row_height * 10);
        subset_list->setMaximumHeight(rowheight + subset_list->frameWidth() * 2);
        subset_list->setMinimumHeight(rowheight + subset_list->frameWidth() * 2);
        subset_list->setSelectionMode(QAbstractItemView::NoSelection);
        od_gridbox->addWidget(subset_list, 1, 0, 1, 2);
        QPushButton* deletebutton = new QPushButton("Delete", subset_mb);
        QPushButton* cancelbutton = new QPushButton("Cancel", subset_mb);
        od_gridbox->addWidget(deletebutton, 2, 0);
        od_gridbox->addWidget(cancelbutton, 2, 1);
        subset_mb->setLayout(od_gridbox);
        subset_mb->setModal(true);

        bool cancelled = false;
        connect(cancelbutton, &QPushButton::clicked, [&cancelled, &subset_mb]() { cancelled = true;  subset_mb->close(); });
        connect(deletebutton, &QPushButton::clicked, subset_mb, &QDialog::close);
        subset_mb->exec();
        if (!cancelled) {
            QDir starting_dir(tempdir);
            for (int i = 0; i < dirlist.size(); i++) {
                QDir del_dir(starting_dir.absolutePath() + QDir::separator() + dirlist[i]);
                logger.info("\t%s\n", del_dir.absolutePath().toLocal8Bit().constData());
                del_dir.removeRecursively();
            }
        }
    }
}

void mtfmapper_app::create_actions(void) {
    open_act = new QAction(tr("&Open..."), this);
    open_act->setShortcut(tr("Ctrl+O"));
    connect(open_act, SIGNAL(triggered()), this, SLOT(open_auto()));
    
    open_roi_act = new QAction(tr("Op&en single edge image(s)..."), this);
    open_roi_act->setShortcut(tr("Ctrl+E"));
    connect(open_roi_act, SIGNAL(triggered()), this, SLOT(open_roi()));
    
    open_manual_roi_act = new QAction(tr("Open with &manual edge selection ..."), this);
    open_manual_roi_act->setShortcut(tr("Ctrl+M"));
    connect(open_manual_roi_act, SIGNAL(triggered()), this, SLOT(open_manual_roi()));
    
    open_focus_act = new QAction(tr("Open &Focus Position image(s)..."), this);
    open_focus_act->setShortcut(tr("Ctrl+F"));
    connect(open_focus_act, SIGNAL(triggered()), this, SLOT(open_focus()));
    
    open_imatest_act = new QAction(tr("Open &Imatest image(s)..."), this);
    open_imatest_act->setShortcut(tr("Ctrl+I"));
    connect(open_imatest_act, SIGNAL(triggered()), this, SLOT(open_imatest_chart()));
    
    exit_act = new QAction(tr("E&xit"), this);
    exit_act->setShortcut(tr("Ctrl+Q"));
    connect(exit_act, SIGNAL(triggered()), this, SLOT(close()));
    
    prefs_act = new QAction(tr("&Preferences"), this);
    prefs_act->setShortcut(tr("Ctrl-P"));
    connect(prefs_act, SIGNAL(triggered()), settings, SLOT( open() ));

    help_act = new QAction(tr("&Help"), this);
    help_act->setShortcut(tr("Ctrl-H"));
    connect(help_act, SIGNAL(triggered()), help, SLOT( open() ));

    about_act = new QAction(tr("&About"), this);
    connect(about_act, SIGNAL(triggered()), about, SLOT( open() ));
}

void mtfmapper_app::view_image(const QString& fname) {
    bool opencv_succeeded = true;
    try {
        img_viewer->load_image(fname);
    }
    catch (const cv::Exception&) {
        // too bad, OpenCV cannot handle it either ...
        opencv_succeeded = false;
    }
    if (!opencv_succeeded) {
        QMessageBox::information(
            this, tr("Image Viewer"),
            tr("Cannot load %1.").arg(fname)
        );
        return;
    }
} 

void mtfmapper_app::open_auto() {
    open_action(false);
}

void mtfmapper_app::open_roi() {
    open_action(true, false);
}

void mtfmapper_app::open_manual_roi() {
    manual_queue_active = false;
    manual_queue_skip_all = false;
    open_action(false, false, false, true);
}

void mtfmapper_app::open_focus() {
    open_action(false, true);
}

void mtfmapper_app::open_imatest_chart() {
    open_action(false, false, true);
}
 
void mtfmapper_app::open_action(bool roi, bool focus, bool imatest, bool manual_roi) {

    QFileDialog* open_dialog = new QFileDialog(this, tr("Select input files"), QString(), QString());
    open_dialog->setOption(QFileDialog::DontUseNativeDialog);
    
    if (!(roi || focus)) {
        QGroupBox* v4GroupBox = new QGroupBox(tr("Select desired MTF Mapper outputs to produce:"));
        QGridLayout* ft_gridbox = new QGridLayout();
        if (ft_gridbox) {
            // add pixel_pitch QComboBox
            if (manual_roi) {
                ft_gridbox->addWidget(tb_img_annotated, 0, 0);
                ft_gridbox->addWidget(tb_img_ca, 0, 1);
            } else {
                ft_gridbox->addWidget(tb_img_annotated, 0, 0);
                ft_gridbox->addWidget(tb_img_profile, 0, 1);
                ft_gridbox->addWidget(tb_img_gridimg, 0, 2);
                ft_gridbox->addWidget(tb_img_lensprofile, 1, 0);
                ft_gridbox->addWidget(tb_img_orientation, 1, 1);
                ft_gridbox->addWidget(tb_img_ca, 1, 2);
            }
        }
        v4GroupBox->setLayout(ft_gridbox);

        QGridLayout* od_gridbox = qobject_cast<QGridLayout*>(open_dialog->layout());
        od_gridbox->addWidget(v4GroupBox);
    }
    
    open_dialog->setFileMode(QFileDialog::FileMode::ExistingFiles);
    
    // use the state from the settings menu as a starting point
    tb_img_annotated->setCheckState(settings->io->cb_annotation->checkState());
    tb_img_profile->setCheckState(settings->io->cb_profile->checkState());
    tb_img_gridimg->setCheckState(settings->io->cb_grid->checkState());
    tb_img_lensprofile->setCheckState(settings->io->cb_lensprofile->checkState());
    tb_img_orientation->setCheckState(settings->io->cb_orientation->checkState());
    tb_img_ca->setCheckState(settings->io->cb_ca_active->checkState());

    if (open_dialog->exec()) {
        // write state back to settings menu
        settings->io->cb_annotation->setCheckState(tb_img_annotated->checkState());
        settings->io->cb_profile->setCheckState(tb_img_profile->checkState());
        settings->io->cb_grid->setCheckState(tb_img_gridimg->checkState());
        settings->io->cb_lensprofile->setCheckState(tb_img_lensprofile->checkState());
        settings->io->cb_orientation->setCheckState(tb_img_orientation->checkState());
        settings->io->cb_ca_active->setCheckState(tb_img_ca->checkState());
        settings->set_gnuplot_img_width(int(img_viewer->size().height()*1.3));

        input_files = open_dialog->selectedFiles();
        if (input_files.size() > 0) {

            QStringList labels;
            labels.push_back(QString("Data set"));
            dataset_contents.setHorizontalHeaderLabels(labels);
            abort_button->show();
            
            Processor_state ps(
                input_files, 
                settings->helpers->get_gnuplot_binary(),
                settings->helpers->get_exiv2_binary(),
                settings->get_argument_string(focus),
                Raw_developer_factory::build(*settings->helpers)
            );
            
            ps.set_single_roi_mode(roi);
            ps.set_focus_mode(focus);
            ps.set_imatest_mode(imatest);
            ps.set_manual_roi_mode(manual_roi);
            emit submit_batch(ps);
        }
    }
    
}
 
void mtfmapper_app::dataset_selected(const QModelIndex& index) {
    // Horrible hack to determine the index of the actual 
    // filename associated with this entry. There must be a 
    // better way ...
    int count_before = 0;
    int parent_row = 0;
    if (index.parent() != QModelIndex()) {
        parent_row = index.parent().row();
        count_before = index.row() + 1;
    } else {
        parent_row = index.row();
    }
    
    for (int row=parent_row-1; row >= 0; row--) {
        QStandardItem* current_dataset_item = dataset_contents.item(row);
        count_before += current_dataset_item->rowCount() + 1;
    }
    active_annotated_filename = std::pair<QString, QStandardItem*>("", nullptr);
    if (dataset_contents.itemFromIndex(index)->isEnabled()) {
        view_image(dataset_files.at(count_before));
        display_exif_properties(count_before);
        if (dataset_contents.itemFromIndex(index)->text().compare(QString("annotated")) == 0) {
            img_viewer->setToolTip(zoom_scroll_tt + "\n\n" + annotated_tt);
            img_viewer->set_clickable(true);
            sfr_list.clear();
            
            QString sfr_source = QFileInfo(dataset_files.at(count_before)).dir().path() + QString("/serialized_edges.bin");
            // go and fetch the corresponding "edge_sfr" entries
            FILE* fin = fopen(sfr_source.toLocal8Bit().constData(), "rb");
            if (!fin) {
                logger.error("Could not open serialized edge info file [%s], SFR profiles not available\n", sfr_source.toLocal8Bit().constData());
                logger.flush();
            } else {
                // read header
                size_t edge_count = 0;
                Job_metadata job_metadata;
                bool success = Edge_info::deserialize_header(fin, edge_count, job_metadata);
                
                if (!success) {
                    logger.error("%s\n", "Could not read serialized edge info header.");
                    logger.flush();
                } else {
                    size_t edges_processed = 0;
                    while (success && edges_processed < edge_count) {
                        Edge_info b = Edge_info::deserialize(fin, success);
                        if (success) {
                            b.set_metadata(job_metadata);
                            sfr_list.push_back(Sfr_entry(b.centroid.x, b.centroid.y, b));
                            edges_processed++;
                        }
                    }
                    fclose(fin);
                    active_annotated_filename = std::pair<QString, QStandardItem*>(dataset_files.at(count_before), dataset_contents.itemFromIndex(index));
                }
            }
        } else {
            img_viewer->setToolTip(zoom_scroll_tt);
            img_viewer->set_clickable(false);
            sfr_list.clear();
        }
    }
}

void mtfmapper_app::dataset_selected_changed(const QModelIndex& i1, const QModelIndex& i2 ATTRIBUTE_UNUSED) {
    dataset_selected(i1);
}
 
void mtfmapper_app::parent_item(QString s, QString f, QString tempdir) {
    
    // take a peek at the edge info header to obtain the job metadata, which
    // will reveal the channel type
    QIcon channel_icon;
    FILE* fin = fopen((tempdir + QString("/serialized_edges.bin")).toLocal8Bit().constData(), "rb");
    if (fin) {
        Job_metadata metadata;
        size_t dummy_count = 0;
        if (Edge_info::deserialize_header(fin, dummy_count, metadata)) {
            switch (metadata.bayer) {
            case Bayer::bayer_t::NONE: 
                if (metadata.channels > 1) {
                    channel_icon.addFile(":/Icons/Channel_RGB_to_L");
                } else {
                    channel_icon.addFile(":/Icons/Channel_Gray");
                }
                break;
            case Bayer::bayer_t::RED: 
                channel_icon.addFile(":/Icons/Channel_Red");
                break;
            case Bayer::bayer_t::GREEN: 
                channel_icon.addFile(":/Icons/Channel_Green");
                break;
            case Bayer::bayer_t::BLUE: 
                channel_icon.addFile(":/Icons/Channel_Blue");
                break;
            }
        }
        fclose(fin);
    }

    current_dataset_item = new QStandardItem(channel_icon, s);
    dataset_files.push_back(f);
}

void mtfmapper_app::child_item(QString s, QString f) {
    QStandardItem* child = new QStandardItem(s);
    child->setEditable(false);
    current_dataset_item->appendRow(child);
    dataset_files.push_back(f);
    exif_properties.push_back(exif_properties.back());
}

void mtfmapper_app::close_item(void) {
    dataset_contents.appendRow(current_dataset_item);
    datasets->setModel(&dataset_contents);
}

void mtfmapper_app::item_for_deletion(QString s) {
    tempfiles_to_delete.push_back(s);
}


void mtfmapper_app::processor_completed(void) {
    abort_button->hide();
}

void mtfmapper_app::enable_clear_button(void) {
    clear_button->setEnabled(true);
}

void mtfmapper_app::disable_clear_button(void) {
    clear_button->setEnabled(false);
}

void mtfmapper_app::clear_button_pressed(void) {
    clear_temp_files();
    dataset_contents.clear(); 
    dataset_files.clear();
    exif_properties.clear();
    sfr_pip_map.clear();
    sfr_dialog->clear();
    sfr_dialog->close();

    img_comment_value->setText("N/A");
    af_ft_value->setText("N/A");
    focal_length_value->setText("N/A");
    focus_distance_value->setText("N/A");

    img_viewer->load_image(icon_image);

    clear_button->setEnabled(false);
}

void mtfmapper_app::enable_save_button(void) {
    save_button->setEnabled(true);
    save_subset_button->setEnabled(true);
}

void mtfmapper_app::disable_save_button(void) {
    save_button->setEnabled(false);
    save_subset_button->setEnabled(false);
}


void mtfmapper_app::enable_file_open(void) {
    open_act->setEnabled(true);
    open_manual_roi_act->setEnabled(true);
    open_roi_act->setEnabled(true);
    open_focus_act->setEnabled(true);
    open_imatest_act->setEnabled(true);
    exit_act->setEnabled(true);
}

void mtfmapper_app::disable_file_open(void) {
    open_act->setEnabled(false);
    open_manual_roi_act->setEnabled(false);
    open_roi_act->setEnabled(false);
    open_focus_act->setEnabled(false);
    open_imatest_act->setEnabled(false);
    exit_act->setEnabled(false);
}

void mtfmapper_app::save_button_pressed(void) {
    save_action(false);
}

void mtfmapper_app::save_subset_button_pressed(void) {
    save_action(true);
}

void mtfmapper_app::save_action(bool subset) {
    // lock file->open to prevent messing with the dataset list
    bool open_was_enabled = open_act->isEnabled();
    disable_file_open();

    std::map<std::string, int> keepers;
    bool cancelled = false;
    if (subset) {
        QListWidget* subset_list = new QListWidget();
        for (int i = 0; i < dataset_contents.rowCount(); i++) {
            QStandardItem* current_dataset_item = dataset_contents.item(i);
            subset_list->addItem(current_dataset_item->text());
        }
        QDialog* subset_mb = new QDialog(this);
        subset_mb->setWindowTitle("Save result subset");
        QGridLayout* od_gridbox = new QGridLayout();
        od_gridbox->addWidget(new QLabel("Select result subsets to save:"), 0, 0);
        QSize vp_size = img_viewer->maximumViewportSize();
        int rowheight = subset_list->model()->rowCount() * subset_list->sizeHintForRow(0);
        rowheight = std::min(rowheight, vp_size.rheight());
        subset_list->setMaximumHeight(rowheight + subset_list->frameWidth()*2);
        subset_list->setMinimumHeight(rowheight + subset_list->frameWidth()*2);
        subset_list->setSelectionMode(QAbstractItemView::ExtendedSelection);
        od_gridbox->addWidget(subset_list, 1, 0, 1, 2);
        QPushButton* savebutton = new QPushButton("Save", subset_mb);
        QPushButton* cancelbutton = new QPushButton("Cancel", subset_mb);
        od_gridbox->addWidget(savebutton, 2, 0);
        od_gridbox->addWidget(cancelbutton, 2, 1);
        subset_mb->setLayout(od_gridbox);
        subset_mb->setModal(true);
        
        connect(cancelbutton, &QPushButton::clicked, [&cancelled, &subset_mb](){ cancelled = true;  subset_mb->close(); });
        connect(savebutton, &QPushButton::clicked, subset_mb, &QDialog::close);
        subset_mb->exec();

        if (!cancelled) {
            QList<QListWidgetItem *> selection = subset_list->selectedItems();
            for (int i = 0; i < selection.size(); i++) {
                keepers[selection[i]->text().toStdString()] = 1;
            }
        }
        delete subset_mb;
    }
    if (!cancelled) {
        QString save_path = QFileDialog::getExistingDirectory(
            this,
            tr("Choose directory to save results in"),
            QDir::homePath()
        );
        int overwrite_count = 0;
        int idx = 0;
        std::vector< std::pair<QString, QString> > copy_list;
        for (int i = 0; i < dataset_contents.rowCount(); i++) {
            QStandardItem* current_dataset_item = dataset_contents.item(i);
            QString prefix = current_dataset_item->text();
            idx++; // skip the actual image file
            for (int j = 0; j < current_dataset_item->rowCount(); j++) {
                QStandardItem* current_child = current_dataset_item->child(j);
                QString dest_fname = save_path + "/" + prefix;
                QString src_fname = dataset_files.at(idx++);
                if (src_fname.endsWith(".jpg")) {
                    dest_fname += "_" + current_child->text() + ".jpg";
                } else {
                    dest_fname += "_" + current_child->text() + ".png";
                }
                dest_fname.replace(' ', '_');
                bool copy_allowed = false;
                if (subset) {
                    if (keepers.find(prefix.toStdString()) != keepers.end()) {
                        copy_allowed = true;
                    }
                } else {
                    copy_allowed = true;
                }
                if (copy_allowed) {
                    copy_list.push_back(std::make_pair(src_fname, dest_fname));
                    logger.info("cp [%s] [%s]\n", src_fname.toLocal8Bit().constData(), dest_fname.toLocal8Bit().constData());
                    if (QFile::exists(dest_fname)) {
                        overwrite_count++;
                    }
                    vector<QString> raw_names = {"edge_mtf_values.txt", "edge_sfr_values.txt", "raw_esf_values.txt", "raw_psf_values.txt", 
                        "chromatic_aberration.txt", "fiducial_correspondence.txt"};
                    if (j == 0) {
                        for (QString name: raw_names) {
                            QString raw_fname = QFileInfo(src_fname).path() + "/" + name;
                            QString dest_fname = save_path + "/" + prefix + "_" + name;
                            dest_fname.replace(' ', '_');
                            if (QFile::exists(raw_fname)) {
                                copy_list.push_back(std::make_pair(raw_fname, dest_fname));
                                logger.info("cp [%s] [%s]\n", raw_fname.toLocal8Bit().constData(), dest_fname.toLocal8Bit().constData());
                            }
                        }
                    }
                }
            }
        }


        bool overwrite_ok = true;
        if (overwrite_count > 0) {
            QMessageBox::StandardButton response = QMessageBox::question(
                this,
                QString("Saving results"),
                overwrite_count == 1 ?
                tr("One of the output files already exists. Do you want to overwrite it?") :
                tr("Some output files already appear to exist. Do you want to overwrite %1 files?").arg(overwrite_count),
                QMessageBox::Yes | QMessageBox::No,
                QMessageBox::Yes
                );
            if (response == QMessageBox::No) {
                overwrite_ok = false;
            }
        }
        for (auto& cp : copy_list) {
            if (QFile::exists(cp.second)) {
                if (overwrite_ok) {
                    QFile::remove(cp.second);
                    QFile::copy(cp.first, cp.second);
                }
            }
            else {
                QFile::copy(cp.first, cp.second);
            }
        }
    }

    if (open_was_enabled) {
        enable_file_open();
    } // otherwise the worker thread will have to re-enable the open action
}

void mtfmapper_app::display_exif_properties(int index) {
    std::shared_ptr<Exiv2_property> props = exif_properties.at(index);
    img_comment_value->setText(props->get_comment());
    af_ft_value->setText(props->get_af_tune());
    focus_distance_value->setText(props->get_focus_distance());

    focal_length_value->setText(props->get_focal_length() + " / " + props->get_aperture());
}

void mtfmapper_app::populate_exif_info_from_file(QString s, QString tempdir) {

    std::shared_ptr<Exiv2_property> props = std::shared_ptr<Exiv2_property>(
        new Exiv2_property(settings->helpers->get_exiv2_binary(), s, tempdir + "/exifinfo.txt")    
    );
    exif_properties.push_back(props);

    // actually, we could delete it right away ...
    item_for_deletion(tempdir + QString("/exifinfo.txt"));
}

void mtfmapper_app::check_if_helpers_exist(void) {
    bool gnuplot_exists = QFile::exists(settings->helpers->get_gnuplot_binary());
    bool exiv_exists = QFile::exists(settings->helpers->get_exiv2_binary());
    bool dcraw_exists = QFile::exists(settings->helpers->get_dcraw_binary());
    // TODO: check for libraw

    if (!gnuplot_exists) {
        QMessageBox::warning(
            this, 
            QString("gnuplot helper"), 
            QString("gnuplot helper executable not found. Please configure this in the settings.")
        );
    }

    if (!exiv_exists) {
        QMessageBox::warning(
            this, 
            QString("Exiv2 helper"), 
            QString("Exiv2 helper executable not found. Please configure this in the settings.")
        );
    }

    if (!dcraw_exists) {
        QMessageBox::warning(
            this, 
            QString("dcraw helper"), 
            QString("dcraw helper executable not found. Please configure this in the settings.")
        );
    }
}

bool mtfmapper_app::edge_selected(int px, int py, bool /*ctrl_down*/, bool shift_down) {
    bool found = false;
    if (sfr_list.size() > 0) {
        
        size_t close_idx = 0;
        double close_dist = 1e50;
        for (size_t i=0; i < sfr_list.size(); i++) {
            double d = sfr_list[i].distance(px, py);
            if (d < close_dist) {
                close_dist = d;
                close_idx = i;
            }
        }
        
        if (close_dist < 50) {
            found = true;
            if (sfr_dialog->isVisible()) {
                sfr_dialog_geom = sfr_dialog->geometry();
            }
            if (shift_down) {
                // if we are moving the pip to a new image, clear the old one
                if (sfr_dialog->pip_number() == 3) {
                    for (auto& p: sfr_pip_map) {
                        p.second.first &= ~(1 << 2);
                    }
                }
            } else {
                sfr_dialog->clear();
                sfr_pip_map.clear();
            }
            sfr_dialog->add_entry(sfr_list[close_idx]);
            if (sfr_dialog_geom != QRect(-1, -1, 0, 0)) {
                sfr_dialog->setGeometry(sfr_dialog_geom);
            }
            sfr_dialog->show();
        }
    }
    
    if (found) {
        // add pip to the correct file
        if (active_annotated_filename.first.length() > 0 &&
            sfr_dialog->pip_number() > 0) {
            
            sfr_pip_map[active_annotated_filename.first].first |= 1 << (sfr_dialog->pip_number() - 1);
            sfr_pip_map[active_annotated_filename.first].second = active_annotated_filename.second;
        }
        update_pip_icons();
    }
    return found;
}

void mtfmapper_app::closeEvent(QCloseEvent* event) {
    if (sfr_dialog && sfr_dialog->isVisible()) {
        sfr_dialog->close();
    } 
    if (about && about->isVisible()) {
        about->close();
    }
    if (help && help->isVisible()) {
        help->close();
    }
    QMainWindow::closeEvent(event);
}

void mtfmapper_app::set_cache_size(int size) {
    img_panel->set_cache_size(uint64_t(size)*uint64_t(1024*1024));
}

// this slot is used to update the SFR dialog in case someone toggled
// the lp/mm checkbox, which effectively changes the x axis of the SFR plots
void mtfmapper_app::settings_saved(void) {
    if (sfr_dialog) {
        sfr_dialog->repaint();
    }
}

void mtfmapper_app::mtfmapper_call_failed(Worker_thread::failure_t failure, const QString& input_file) {
    if (failure == Worker_thread::failure_t::UNSPECIFIED) {
        if (settings->peek_argument_line().trimmed().length() > 0) {
            QMessageBox call_failed;
            call_failed.setIcon(QMessageBox::Critical);
            call_failed.setWindowTitle("MTF Mapper call failed");
            call_failed.setText(
                "MTF Mapper call failed.\n\n"
                "It is likely that this error has been caused by an invalid option in the Preferences/Arguments field.\n\n"
                "Reset Arguments field in the Preferences dialogue?"
            );
            call_failed.setStandardButtons(QMessageBox::Reset | QMessageBox::Ignore);
            int ret = call_failed.exec();
            switch (ret) {
            case QMessageBox::Reset:
                settings->reset_argument_line();
                break;
            case QMessageBox::Ignore:
            default:
                break;
            }
        }
        else {
            QString logname = QDir::tempPath() + QDir::separator() + QString("mtfmapperlog.txt");
            logger.flush();
            QMessageBox call_failed;
            call_failed.setIcon(QMessageBox::Critical);
            call_failed.setWindowTitle("MTF Mapper call failed");
            call_failed.setText(
                "MTF Mapper call failed.\n\n"
                "Please contact fvdbergh@gmail.com for assistance.\n"
                "Attach the file '" + logname + "' if possible.\n"
            );
            call_failed.setStandardButtons(QMessageBox::Ok);
            call_failed.exec();
        }
    } else {
        QMessageBox call_failed;
        call_failed.setIcon(QMessageBox::Critical);
        call_failed.setWindowTitle("MTF Mapper call failed");
        switch (failure) {
        case Worker_thread::failure_t::IMAGE_OPEN_FAILURE:
            call_failed.setText(
                QString(
                    "MTF Mapper call failed, indicating that the input image [%1] could not be opened.\n\n"
                    "Maybe that was not a valid image file.\n"
                ).arg(input_file)
            );
            break;
        case Worker_thread::failure_t::UNSUPPORTED_IMAGE_ENCODING:
            call_failed.setText(
                QString(
                    "MTF Mapper call failed, indicating that the input image [%1] could not be opened.\n\n"
                    "Only 8-bit unsigned and 16-bit unsigned image encodings are supported.\n"
                ).arg(input_file)
            );
            break;
        case Worker_thread::failure_t::NO_TARGETS_FOUND:
            call_failed.setText(
                QString(
                    "MTF Mapper call failed, indicating that no valid target objects were found.\n\n"
                    "Maybe the input image [%1] did not contain valid targets?\n"
                    "Otherwise you could try lowering the Threshold value in the Preferences dialog.\n"
                ).arg(input_file)
            );
            break;
        case Worker_thread::failure_t::RAW_DEVELOPER_FAILURE:
            call_failed.setText(
                QString(
                    "Raw developer call failed.\n\n"
                    "Maybe the input image [%1] was corrupted?\n"
                ).arg(input_file)
            );
            break;
        default:
            call_failed.setText(
                "MTF Mapper call failed.\n\n"
                "No further information available at this time.\n"
            );
            break;
        }
        call_failed.setStandardButtons(QMessageBox::Ok);
        call_failed.exec();
    }
}

void mtfmapper_app::add_manual_roi_file(const Processing_command& new_command) {
    {
        std::lock_guard<std::mutex> lock(mq_mutex);
        manual_roi_commands.push(new_command);
        emit manual_roi_queue_size((int)manual_roi_commands.size());
        if (mq_busy) {
            return;
        }
        mq_busy = true;
    }
    
    while (true) {
        bool queue_emptied = true;
        Processing_command command;
        {
            std::lock_guard<std::mutex> lock(mq_mutex);
            if (!manual_roi_commands.empty()) {
                command = manual_roi_commands.front();
                manual_roi_commands.pop();
                queue_emptied = false;
            }
        }
        emit manual_roi_queue_size((int)manual_roi_commands.size());
        if (queue_emptied) {
            break;
        }
        
        if (manual_queue_skip_all) {
            processor.remove_file_in_flight();
            item_for_deletion(command.tmp_dirname + QString("/exifinfo.txt"));
            continue;
        }
        
        QString autoload_roi_filename = QString("%1/mtfmapper_autoload.roi").arg(QDir::tempPath());
        QString roi_filename = command.tmp_dirname + "/manual.roi";
        
        // ensure the exposed main window repaints
        update();
        QCoreApplication::processEvents();
        
        if (manual_queue_active) {
            Processing_command modified_command(command);
            modified_command.arguments << "--roi-file" << manual_queue_roi_filename;
            modified_command.set_state(Processing_command::state_t::READY);
            processor.submit_processing_command(modified_command);
        } else {
            disable_clear_button();
            
            bool load_success = edge_select_dialog->load_image(
                command.img_filename, 
                command.arguments,
                command.exif_filename // original raw filename, if applicable
            );
            edge_select_dialog->set_roi_file(roi_filename);
            
            bool dialog_success = false;
            if (load_success) {
                if (edge_select_dialog->autoload_roi()) {
                    edge_select_dialog->load_rois(autoload_roi_filename);
                }
                
                if (esd_geom != QRect(-1, -1, 0, 0)) {
                    edge_select_dialog->setGeometry(esd_geom);
                }
                #ifdef _WIN32
                edge_select_dialog->show();
                #else
                edge_select_dialog->showNormal(); // prevents dialog from remaining maximized on X11
                #endif
                edge_select_dialog->raise();
                edge_select_dialog->activateWindow();
                
                // run a local event loop to simulate the behaviour of exec(), but without
                // blocking the main window
                QEventLoop q;
                connect(edge_select_dialog, SIGNAL(finished(int)), &q, SLOT(quit()));
                connect(edge_select_dialog, SIGNAL(finished(int)), this, SLOT(enable_clear_button()));
                q.exec();
                
                esd_geom = edge_select_dialog->geometry();
                if (edge_select_dialog->result() == QDialog::Accepted) {
                    dialog_success = true;
                    Processing_command modified_command(command);
                    modified_command.arguments << "--roi-file" << roi_filename;
                    modified_command.set_state(Processing_command::state_t::READY);
                    processor.submit_processing_command(modified_command);
                    if (edge_select_dialog->autoload_roi()) {
                        edge_select_dialog->save_rois(autoload_roi_filename);
                        item_for_deletion(autoload_roi_filename);
                    } 
                }
            } 
            if (!load_success || !dialog_success) {
                processor.remove_file_in_flight();
                item_for_deletion(command.tmp_dirname + QString("/exifinfo.txt"));
            }
            
            if (!load_success) {
                emit mtfmapper_call_failed(Worker_thread::failure_t::IMAGE_OPEN_FAILURE, command.img_filename);
            }
        }
    }
    {
        std::lock_guard<std::mutex> lock(mq_mutex);
        mq_busy = false;
    }
}

void mtfmapper_app::dragEnterEvent(QDragEnterEvent* evt) {
    bool has_local_files = false;
    if (evt->mimeData()->hasUrls()) {
        auto urls = evt->mimeData()->urls();
        for (auto url : urls) {
            has_local_files |= url.isLocalFile();
        }
        if (has_local_files && open_act->isEnabled()) {
            evt->accept();
        }
    }
}

void mtfmapper_app::dropEvent(QDropEvent* evt) {
    QStringList dropped_files;
    auto urls = evt->mimeData()->urls();
    for (auto url : urls) {
        if (url.isLocalFile()) {
            dropped_files.append(url.toLocalFile());
        }
    }

    if (open_act->isEnabled() && dropped_files.size() > 0) { // only process the drop event if we are not currently processing
        // TODO: not ideal to have a copy of the same code from File->Open here, really has to be refactored

        QStringList labels;
        labels.push_back(QString("Data set"));
        dataset_contents.setHorizontalHeaderLabels(labels);
        abort_button->show();
        
        Processor_state ps(
            dropped_files, 
            settings->helpers->get_gnuplot_binary(),
            settings->helpers->get_exiv2_binary(),
            settings->get_argument_string(false),
            Raw_developer_factory::build(*settings->helpers)
        );
        
        emit submit_batch(ps);
    }
}

void mtfmapper_app::update_progress(int val) {
    int prog = progress->maximum() - val;
    progress->setValue(prog);
}

void mtfmapper_app::enable_manual_queue_processing(QString filename) {
    manual_queue_roi_filename = filename;
    manual_queue_active = true;
    emit item_for_deletion(filename);
}

void mtfmapper_app::enable_manual_queue_skip(void) {
    manual_queue_skip_all = true;
}

void mtfmapper_app::update_pip_icons(void) {
    // first, clear all old icons
    std::function<void(QStandardItemModel*, QModelIndex)> loop_each;
    loop_each = [this, &loop_each](QStandardItemModel* model, QModelIndex parent) {
        for (int r=0; r < model->rowCount(parent); r++) {
            QModelIndex index = model->index(r, 0, parent);
            QString name = model->data(index).toString();
            
            if (name.compare(QString("annotated")) == 0) {
                model->itemFromIndex(index)->setIcon(QIcon());
            }
            
            if (model->hasChildren(index)) {
                loop_each(model, index);
            }
        }
    };
    loop_each(&dataset_contents, QModelIndex());
    
    // now draw the new ones
    for (auto& w: sfr_pip_map) {
        if (w.second.first > 0 && w.second.first < 8) {
            QIcon pip_icon;
            QString icon_name = QString(":/Icons/SFR_pip_%0").arg(w.second.first);
            pip_icon.addFile(icon_name);
            w.second.second->setIcon(pip_icon);
        }
    }
}

void mtfmapper_app::sfr_dialog_closed(QRect r) {
    sfr_dialog_geom = r;
    sfr_pip_map.clear();
    update_pip_icons();
}

