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
#include "worker_thread.h"

#include <string>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

#include <QFileInfo>
#include <QSharedPointer>
#include <QCoreApplication>
#include <QProcess>
#include "mtfmapper_app.h"

using std::cout;
using std::endl;
using std::string;

#if QT_VERSION >= QT_VERSION_CHECK(5,14,0)
    #define EMPTY_STRING_PARTS Qt::SkipEmptyParts
#else
    #define EMPTY_STRING_PARTS QString::SkipEmptyParts
#endif

Worker_thread::Worker_thread(QWidget* parent) 
: parent(dynamic_cast<mtfmapper_app*>(parent)), tempdir_number(0), 
  dev_done(false), pc_done(false), batch_done(false), abort(false) {
  
  // initialize the threads a bit later when we are sure the mutexes have been initialized
  dev_thread = std::thread(&Worker_thread::developer_run, this);
  pc_thread = std::thread(&Worker_thread::processor_run, this);
  batch_thread = std::thread(&Worker_thread::batch_run, this);
}

Worker_thread::~Worker_thread(void) {
    abort = true;
    batch_done = true;
    batch_cv.notify_one();
    batch_thread.join();
    dev_done = true;
    dev_cv.notify_one();
    dev_thread.join();
    pc_done = true;
    pc_cv.notify_one();
    pc_thread.join();
}

void Worker_thread::receive_batch(Processor_state state) {
    abort = false;
    QString arguments = state.updated_arguments();
    
    // start by assuming all files will pass raw development
    fif_add(state.input_files.size());
    {
        std::lock_guard<std::mutex> lock(batch_mutex);
        batch_queue.push(state);
    }
    
    for (int i=0; i < state.input_files.size(); i++) {
        std::lock_guard<std::mutex> lock(dev_mutex);
        QString tempdir = tr("%1/mtfmappertemp_%2").arg(QDir::tempPath()).arg(tempdir_number++);
        dev_queue.push(Input_file_record(state.input_files.at(i), arguments, tempdir, state.raw_developer));
    }
    // this should be enough to trigger processing of 'state'
    dev_cv.notify_one();
    batch_cv.notify_one();
}

void Worker_thread::batch_run(void) {
    while (!batch_done) {
        {
            std::unique_lock<std::mutex> lock(batch_mutex);
            batch_cv.wait(lock, [this]{ return batch_done || !batch_queue.empty(); });
        
            if (batch_done) { // for shutting down
                break;
            }
        }
        
        bool batch_queue_done = false;
        while (!batch_queue_done) {
        
            Processor_state state;
            {
                std::lock_guard<std::mutex> lock(batch_mutex);
                state = batch_queue.front();
                batch_queue.pop();
            }
            
            if (!abort) {    
                if (state.is_valid()) {
                    int remaining_files = state.input_files.size();
                    
                    while (!abort && remaining_files > 0) { 
                        Input_file_record input;
                        {
                            std::unique_lock<std::mutex> file_lock(file_mutex);
                            // while blocking in wait, file_mutex is not locked, so developer_run can obtain a lock
                            // to add new entries without blocking indefinitely
                            file_cv.wait(file_lock, [this]{ return !file_queue.empty(); });
                            
                            input = file_queue.front();
                            file_queue.pop();
                            remaining_files--;
                        }
                        
                        if (input.is_valid()) {
                            if (QFileInfo(input.get_output_name()).size() == 0 || input.get_state() == Input_file_record::state_t::FAILED) {
                                add_failure(failure_t::RAW_DEVELOPER_FAILURE, input.get_input_name());
                                fif_add(-1); // one less file to worry about
                                continue;
                            }
                        
                            QStringList mapper_args;
                            
                            // if a "--focal-ratio" setting is already present, then assume this
                            // was a user-provided override, otherwise try to calculate it from the EXIF data
                            if (!input.get_arguments().contains("--focal-ratio")) {
                                Exiv2_property props(state.exiv2_binary, input.get_input_name(), input.get_temp_dir() + "/exifinfo.txt");
                                mapper_args << "--focal-ratio" << props.get_focal_ratio();
                            }

                            mapper_args << "--gnuplot-executable " + state.gnuplot_binary << input.get_output_name() << input.get_temp_dir() 
                                << "--logfile " + input.get_temp_dir() + "/log.txt" 
                                << input.get_arguments().split(QRegExp("\\s+"), EMPTY_STRING_PARTS);
                            
                            Processing_command pc(
                                QCoreApplication::applicationDirPath() + "/mtf_mapper",
                                mapper_args,
                                input.get_output_name(), // output name of raw developer is input to mapper
                                input.get_temp_dir(),
                                input.get_input_name(), // original input file
                                state.initial_state()
                            );
                            
                            {
                                std::lock_guard<std::mutex> pc_lock(pc_mutex);
                                pc_queue.push(pc);
                            }
                            pc_cv.notify_one();
                        }
                    }
                    
                    {
                        // pass on all failures so far (usually raw developer failures before processing completes)
                        std::lock_guard<std::mutex> lock(failure_mutex);
                        for (const auto& failure : failure_list) {
                            emit mtfmapper_call_failed(failure.first, failure.second);
                        }
                        failure_list.clear();
                    }
                
                }
            } 
            
            {
                std::lock_guard<std::mutex> lock(batch_mutex);
                batch_queue_done = batch_queue.empty();
            }
        }
    }
    emit work_finished();
}

void Worker_thread::receive_abort() {
    abort = true;
}

void Worker_thread::process_command(const Processing_command& command) {
    
    QProcess process;
    process.setProgram(command.program);
    process.setArguments(command.arguments);
    const QString& tempdir = command.tmp_dirname;
    
    bool q_output_requested = process.arguments().contains("-q");
    bool e_output_requested = process.arguments().contains("-e") || process.arguments().contains("--esf");

    logger.debug("%s\n", "arguments to mtf mapper:");
    for (int kk = 0; kk < command.arguments.size(); kk++) {
        logger.debug("[%d]=%s\n", kk, command.arguments.at(kk).toLocal8Bit().constData());
    }
    process.start();
    process.waitForFinished(-1);
    int rval = process.exitStatus() == QProcess::NormalExit && process.exitCode() == 0;
    if (!rval) {
        failure_t failure = UNSPECIFIED;
        switch (process.exitCode()) {
        case 2: failure = IMAGE_OPEN_FAILURE; break;
        case 4: failure = NO_TARGETS_FOUND; break;
        case 5: failure = UNSUPPORTED_IMAGE_ENCODING; break;
        default: failure = UNSPECIFIED; break;
        }
        add_failure(failure, command.img_filename);
    } else {
        emit send_delete_item(tempdir + "/log.txt");

        // this call must come from within the worker thread, since we
        // may have to perform a raw conversion in the worker thread
        // which would cause the display image filename (root of each data set object)
        // to differ from the file containing the exif info
        
        emit send_exif_filename(command.exif_filename, tempdir);
        QString fname(QFileInfo(command.img_filename).completeBaseName());
        emit send_parent_item(fname, command.img_filename, tempdir);
        
        if (q_output_requested) {
            emit send_delete_item(tempdir + QString("/edge_sfr_values.txt"));
            emit send_delete_item(tempdir + QString("/edge_mtf_values.txt"));
            emit send_delete_item(tempdir + QString("/edge_line_deviation.txt"));
            if (QFile().exists(tempdir + QString("/serialized_edges.bin"))) {
                emit send_delete_item(tempdir + QString("/serialized_edges.bin"));
            }
        }
        
        if (e_output_requested) {
            if (QFile().exists(tempdir + QString("/raw_esf_values.txt"))) {
                emit send_delete_item(tempdir + QString("/raw_esf_values.txt"));
            }
            if (QFile().exists(tempdir + QString("/raw_psf_values.txt"))) {
                emit send_delete_item(tempdir + QString("/raw_psf_values.txt"));
            }
        }
        
        QString an_file = QString("%1/annotated.png").arg(tempdir);
        if (QFile().exists(an_file)) {
            emit send_child_item(QString("annotated"), an_file);
            emit send_delete_item(an_file);
        }
        QString an_file_jpeg = QString("%1/annotated.jpg").arg(tempdir);
        if (QFile().exists(an_file_jpeg)) {
            emit send_child_item(QString("annotated"), an_file_jpeg);
            emit send_delete_item(an_file_jpeg);
        }
        QString pr_file = QString("%1/profile_image.png").arg(tempdir);
        if (QFile().exists(pr_file)) {
            emit send_child_item(QString("profile"), pr_file);
            emit send_delete_item(pr_file);
            emit send_delete_item(tempdir + QString("/profile.gnuplot"));
            emit send_delete_item(tempdir + QString("/profile.txt"));
            emit send_delete_item(tempdir + QString("/profile_peak.txt"));
        }
        QString gi_file = QString("%1/grid_image.png").arg(tempdir);
        if (QFile().exists(gi_file)) {
            emit send_child_item(QString("grid2d"), gi_file);
            emit send_delete_item(gi_file);
            emit send_delete_item(tempdir + QString("/grid.gnuplot"));
            emit send_delete_item(tempdir + QString("/grid.txt"));
        }
        QString gs_file = QString("%1/grid_surface.png").arg(tempdir);
        if (QFile().exists(gs_file)) {
            emit send_child_item(QString("grid3d"), gs_file);
            emit send_delete_item(gs_file);
        }
        QString fp_file = QString("%1/focus_peak.png").arg(tempdir);
        if (QFile().exists(fp_file)) {
            emit send_child_item(QString("focus"), fp_file);
            emit send_delete_item(fp_file);
            emit send_delete_item(tempdir + QString("/profile_curve.txt"));
            emit send_delete_item(tempdir + QString("/profile_points.txt"));
            // we might as well try to delete all the bayer-channel specific outputs, whether they are generated, or not
            emit send_delete_item(tempdir + QString("/green_profile_curve.txt"));
            emit send_delete_item(tempdir + QString("/green_profile_points.txt"));
            emit send_delete_item(tempdir + QString("/blue_profile_curve.txt"));
            emit send_delete_item(tempdir + QString("/blue_profile_points.txt"));
            emit send_delete_item(tempdir + QString("/red_profile_curve.txt"));
            emit send_delete_item(tempdir + QString("/red_profile_points.txt"));
        }
        QString lp_file = QString("%1/lensprofile.png").arg(tempdir);
        if (QFile().exists(lp_file)) {
            emit send_child_item(QString("lensprofile"), lp_file);
            emit send_delete_item(lp_file);
            emit send_delete_item(tempdir + QString("/lensprofile.txt"));
            emit send_delete_item(tempdir + QString("/lensprofile.gnuplot"));
        }
        QString co_file = QString("%1/chart_orientation.png").arg(tempdir);
        if (QFile().exists(co_file)) {
            emit send_child_item(QString("chart orientation"), co_file);
            emit send_delete_item(co_file);
        }
        QString ca_file = QString("%1/ca_image.png").arg(tempdir);
        if (QFile().exists(ca_file)) {
            emit send_child_item(QString("chromatic aberration"), ca_file);
            emit send_delete_item(ca_file);
            emit send_delete_item(tempdir + QString("/ca_grid.gnuplot"));
            emit send_delete_item(tempdir + QString("/ca_grid.txt"));
        }
        if (QFile().exists(tempdir + QString("/chromatic_aberration.txt"))) {
            emit send_delete_item(tempdir + QString("/chromatic_aberration.txt"));
        }
        QString fids_file = QString("%1/fiducial_correspondence.txt").arg(tempdir);
        if (QFile().exists(fids_file)) {
            emit send_delete_item(fids_file);
        }
        QString roi_file = QString("%1/manual.roi").arg(tempdir);
        if (QFile().exists(roi_file)) {
            emit send_delete_item(roi_file);
        }
        
        emit send_close_item();
    }
}

void Worker_thread::developer_run(void) {
    while (!dev_done) {
        {
            std::unique_lock<std::mutex> lock(dev_mutex);
            dev_cv.wait(lock, [this]{ return dev_done || !dev_queue.empty(); });
            
            if (dev_done) { // for shutting down
                break;
            }
        }
        
        bool dev_queue_done = false;
        while (!dev_queue_done) {
            
            Input_file_record state;
            {
                std::lock_guard<std::mutex> lock(dev_mutex);
                state = dev_queue.front();
                dev_queue.pop();
            }
            
            if (!abort) {
            
                QDir().mkdir(state.get_temp_dir());

                QFileInfo fi(state.get_input_name());
                if (state.get_raw_developer()) {
                    if (state.get_raw_developer()->accepts(fi.suffix())) {
                        state.set_output_name(
                            QString(state.get_temp_dir() + "/" + fi.completeBaseName() + QString(".tiff"))
                        );
                        
                        // TODO: we can add a cache here that we can pass to process; if the input file + args match
                        // just copy the previous raw developer output if it still exists. store cache in worker_thread
                        // for persistance
                        int dev_success = state.get_raw_developer()->process(
                            state.get_input_name(), // input file to raw developer
                            state.get_output_name(), // output file of raw developer
                            state.get_arguments().contains(QString("--bayer")) // bayer_mode
                        );
                        
                        if (dev_success) {
                            state.set_state(Input_file_record::state_t::COMPLETED);
                        } else {
                            state.set_state(Input_file_record::state_t::FAILED);
                        }
                        
                        emit send_delete_item(state.get_output_name()); // queued, will only be deleted on program exit, or clear() command
                    } else {
                        state.set_output_name(state.get_input_name());
                        state.set_state(Input_file_record::state_t::COMPLETED);
                    }
                    
                    {
                        std::lock_guard<std::mutex> file_lock(file_mutex);
                        file_queue.push(state);
                    }
                    file_cv.notify_one();
                    
                } else {
                    // this should not happen at runtime
                    logger.error("Raw developer not set in Worker_thread::run()\n");
                    break;
                }
            } else {
                fif_add(-1);
            }
            
            {
                std::lock_guard<std::mutex> lock(dev_mutex);
                dev_queue_done = dev_queue.empty();
            }
        }
    }
}

void Worker_thread::processor_run(void) {
    while (!pc_done) {
        {
            std::unique_lock<std::mutex> lock(pc_mutex);
            pc_cv.wait(lock, [this]{ return pc_done || !pc_queue.empty(); });
        
            if (pc_done) { // for shutting down
                break;
            }
        }
        
        bool pc_queue_done = false;
        while (!pc_queue_done) {
        
            Processing_command cmd;
            {
                std::lock_guard<std::mutex> lock(pc_mutex);
                cmd = pc_queue.front();
                pc_queue.pop();
            }
            
            if (!abort) {    
                if (cmd.is_valid()) {
                    if (cmd.get_state() == Processing_command::state_t::AWAIT_ROI) {
                        send_processing_command(cmd);
                    } else {
                        process_command(cmd);
                        fif_add(-1);
                    }
                }
            } else {
                fif_add(-1);
            }
            
            {
                std::lock_guard<std::mutex> lock(pc_mutex);
                pc_queue_done = pc_queue.empty();
            }
        }
        
        {
            std::lock_guard<std::mutex> lock(failure_mutex);
            for (const auto& failure : failure_list) {
                emit mtfmapper_call_failed(failure.first, failure.second);
            }
            failure_list.clear();
        }
    }
}

void Worker_thread::submit_processing_command(const Processing_command& command) {
    {
        std::lock_guard<std::mutex> pc_lock(pc_mutex);
        pc_queue.push(command);
    }
    pc_cv.notify_one();    
}

void Worker_thread::add_failure(failure_t fail, const QString& fname) {
    std::lock_guard<std::mutex> lock(failure_mutex);
    failure_list.push_back(std::make_pair(fail, fname));
}

void Worker_thread::fif_add(int delta) {
    std::lock_guard<std::mutex> lock(fif_mutex);
    fif_count += delta;
    
    if (delta > 0) {
        send_progress_indicator_max(0, fif_count + 1);
    }
    
    emit send_progress_indicator(fif_count);
    
    if (fif_count <= 0) {
        emit send_all_done();
    }
}

