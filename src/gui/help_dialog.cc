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
#include "help_dialog.h"

#include "common.h"
#include "config.h"

#define HELP_DIR_FORMAT(x) "##x"

Help_dialog::Help_dialog(QWidget* parent ATTRIBUTE_UNUSED) : QDialog(parent) {

    dismiss_button = new QPushButton("Dismiss");
    dismiss_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );

	body = new QTextBrowser(parent);
    #if _WIN32
    if (mtfmapper_HTML_HELP_DIR.compare(".//html/mtf_mapper_user_guide.html") == 0) {
        body_text = QDir(QDir::toNativeSeparators(QCoreApplication::applicationDirPath() + QString("/../html/mtf_mapper_user_guide.html"))).canonicalPath();
    } else {
        body_text = QDir(QString(mtfmapper_HTML_HELP_DIR.c_str())).canonicalPath();
    }
    #else
    body_text = QDir(QString(mtfmapper_HTML_HELP_DIR.c_str())).canonicalPath();
    #endif
	body->setMinimumWidth(600);
    body->setMinimumHeight(200);

    QGridLayout* vlayout = new QGridLayout(this);

    vlayout->addWidget(body, 0, 0, 1, 3);
    vlayout->addWidget(dismiss_button, 1, 1);
    
    connect(dismiss_button, SIGNAL(clicked()), this, SLOT( close() ));
    
    
    setLayout(vlayout);
    setWindowTitle("Help");
}


void Help_dialog::open(void) {
	bool result = QDesktopServices::openUrl(QUrl::fromLocalFile(body_text));
	if (!result) {
		body->setHtml(
			QString(
				"<HTML>"
				"<BODY>"
				"<p style=\"text - align:center\">"
				"Could not launch the browser to open help file.<br>"
				"Try copying the link below, and paste it in your browser.<br>"
				"<a href=\"file://"
			) + body_text +
			QString("\">file://") + body_text +
			QString(
				"</a>"
				"</p>"
				"</BODY>"
				"</HTML>"
			)
		);
		body->show();
		show();
	}
}


