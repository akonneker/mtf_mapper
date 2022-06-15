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

#include "manual_roi_help_dialog.h"


Manual_roi_help_dialog::Manual_roi_help_dialog(QWidget* parent) 
    : QDialog(parent, Qt::WindowSystemMenuHint | Qt::WindowTitleHint | Qt::WindowCloseButtonHint) {

    dismiss_button = new QPushButton("Dismiss");
    dismiss_button->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    
    prev_button = new QPushButton("prev");
    prev_button->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    prev_button->setEnabled(false);
    
    next_button = new QPushButton("next");
    next_button->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    
    page_counter = new QLabel("1");
    page_counter->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    
    pages.push_back(std::shared_ptr<QPixmap>(new QPixmap(":Images/ManualROIHelp")));
    pages.push_back(std::shared_ptr<QPixmap>(new QPixmap(":Images/ManualROIHelp_p2")));
    pages.push_back(std::shared_ptr<QPixmap>(new QPixmap(":Images/ManualROIHelp_p3")));
    
    page_counter->setText(QString("%1/%2").arg(page_number + 1).arg(pages.size()));

    label = new QLabel(this);
    label->setAlignment(Qt::AlignCenter);
    label->setPixmap(*pages[page_number]);

    scroller = new QScrollArea;
    scroller->setWidget(label);
    
    update_page();    
    scroller->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    scroller->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    scroller->setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
    
    QGridLayout* button_layout = new QGridLayout;
    button_layout->setColumnStretch(0, 10);
    button_layout->addWidget(prev_button, 0, 1);
    button_layout->setColumnStretch(1, 0);
    button_layout->addWidget(page_counter, 0, 2);
    button_layout->setColumnStretch(2, 0);
    button_layout->addWidget(next_button, 0, 3);
    button_layout->setColumnStretch(3, 0);
    button_layout->setColumnStretch(4, 10);
    button_layout->addWidget(dismiss_button, 0, 5);
    button_layout->setColumnStretch(5, 0);

    QGridLayout* vlayout = new QGridLayout(this);
    vlayout->addWidget(scroller, 0, 0);
    vlayout->addLayout(button_layout, 1, 0);

    connect(dismiss_button, SIGNAL(clicked()), this, SLOT(close()));
    connect(prev_button, SIGNAL(clicked()), this, SLOT(page_down()));
    connect(next_button, SIGNAL(clicked()), this, SLOT(page_up()));

    setLayout(vlayout);
    setWindowTitle("Manual ROI selection help");
}

void Manual_roi_help_dialog::open(void) {
    show();
}

void Manual_roi_help_dialog::update_page(void) {
    QScreen* screen = QGuiApplication::primaryScreen();
    QRect geom = screen->geometry();

    QPixmap& image = *pages[page_number];
    label->setPixmap(image);
    int swidth = image.width() > geom.width() ? geom.width() - 40 : image.width() + 16;
    int sheight = image.height() > geom.height() ? geom.height() - 80 : image.height() + 16;

    scroller->setMinimumSize(swidth, sheight);
    update();
}

void Manual_roi_help_dialog::page_up(void) {
    page_number = std::min((int)pages.size() - 1, page_number + 1);
    page_counter->setText(QString("%1/%2").arg(page_number + 1).arg(pages.size()));
    if (page_number == (int)pages.size() - 1) {
        next_button->setEnabled(false);
    }
    prev_button->setEnabled(true);
    update_page();
}

void Manual_roi_help_dialog::page_down(void) {
    page_number = std::max(0, page_number - 1);
    page_counter->setText(QString("%1/%2").arg(page_number + 1).arg(pages.size()));
    if (page_number == 0) {
        prev_button->setEnabled(false);
    }
    next_button->setEnabled(true);
    update_page();
}
