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
#ifndef WORKER_THREAD_H
#define WORKER_THREAD_H

#include "processing_command.h"
#include "raw_developer.h"
#include "input_file_record.h"
#include "processor_state.h"

#include <vector>
using std::vector;
#include <chrono>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <atomic>
#include <queue>
using std::queue;

#include <QThread>
#include <QStringList>
#include <QStandardItemModel>
#include <QDir>
#include <QRunnable>

class mtfmapper_app;

class Worker_thread  : public QObject
{
  Q_OBJECT
  
  public:
    typedef enum {
        UNSPECIFIED=0,
        IMAGE_OPEN_FAILURE=1,
        NO_TARGETS_FOUND=2,
        UNSUPPORTED_IMAGE_ENCODING=3,
        RAW_DEVELOPER_FAILURE=4
    } failure_t;

    Worker_thread(QWidget *parent);
    ~Worker_thread(void);
    
    void process_command(const Processing_command& command);
    void submit_processing_command(const Processing_command& command);
    
    void remove_file_in_flight(void) {
        fif_add(-1);
    }
    
  signals:
    void send_parent_item(QString s, QString f, QString tempdir);
    void send_child_item(QString s, QString f);
    void send_close_item(void); 
    void send_delete_item(QString s); 
    void send_exif_filename(QString s, QString tempdir);
    
    void send_progress_indicator(int p);
    void send_progress_indicator_max(int lb, int ub);
    void send_all_done(void);
    void mtfmapper_call_failed(Worker_thread::failure_t failure, const QString& input_file);
    
    void send_processing_command(const Processing_command& command);
    void work_finished(void);
    
  public slots:
    void receive_abort(void);
    void receive_batch(Processor_state state);
    
  private:
    void developer_run(void);
    void processor_run(void);
    void batch_run(void);
    void add_failure(failure_t fail, const QString& fname);
    void fif_add(int delta);
  
    mtfmapper_app* parent;
    
    std::atomic_int tempdir_number;
    
    bool dev_done = false;
    std::thread dev_thread;
    std::mutex dev_mutex;
    std::condition_variable dev_cv;
    queue<Input_file_record> dev_queue;
    
    std::mutex file_mutex;
    std::condition_variable file_cv;
    queue<Input_file_record> file_queue;
    
    bool pc_done = false;
    std::thread pc_thread;
    std::mutex pc_mutex;
    std::condition_variable pc_cv;
    queue<Processing_command> pc_queue;
    
    bool batch_done = false;
    std::thread batch_thread;
    std::mutex batch_mutex;
    std::condition_variable batch_cv;
    queue<Processor_state> batch_queue;
    
    std::mutex failure_mutex;
    vector<std::pair<failure_t, QString>> failure_list;
    
    // files-in-flight
    std::mutex fif_mutex;
    int fif_count = 0;
    
    bool abort;
};

#endif

