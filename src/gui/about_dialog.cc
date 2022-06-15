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
#include "about_dialog.h"

#include "common.h"
#include "config.h"

#define SFI(x) #x
#define VER(x) SFI(x)

About_dialog::About_dialog(QWidget* parent ATTRIBUTE_UNUSED) {

    dismiss_button = new QPushButton("Dismiss");
    dismiss_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );

    body = new QTextEdit(parent);
    body->setHtml(
        "<html>"
        "<body bgcolor=\"#efefef\">"
        "<br><br>"
        "<center><big>MTF Mapper</big></center>"
        "<center>version " VER(mtfmapper_VERSION_MAJOR) "." VER(mtfmapper_VERSION_MINOR) "." VER(mtfmapper_VERSION_SUB) "</center>"
        "<center><small>by</small></center>"
        "<center>Frans van den Bergh</center>"
        "<br><br>"
        "<center>visit MTF Mapper online at </center>"
        "<center><a href=\"http://sourceforge.net/projects/mtfmapper/\">"
          "http://sourceforge.net/projects/mtfmapper/</a>"
        "</center>"
        "</body>"
        "</html>"
    );
    body->setFixedWidth(400);
    body->setReadOnly(true);

    QGridLayout* vlayout = new QGridLayout(this);
    vlayout->addWidget(body, 0, 0, 1, 3);
    vlayout->addWidget(dismiss_button, 1, 1);
    
    connect(dismiss_button, SIGNAL(clicked()), this, SLOT( close() ));
    
    setLayout(vlayout);
    setWindowTitle("About");
}


void About_dialog::open(void) {
    body->show();
    show();
}


