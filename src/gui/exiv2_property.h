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
#ifndef EXIV2_PROPERTY_H
#define EXIV2_PROPERTY_H

#include <QString>

class Exiv2_property {
  public:

    typedef enum {
        NIKON=1,
        CANON,
        OLYMPUS,
        OTHER,
        NONE
    } Exiv_mode;

    Exiv_mode mode;

    Exiv2_property(QString bin_name, QString ifname, QString tfname);

    Exiv_mode get_mode(void) const {
        return mode;
    }

    QString get_af_tune(void) const {
        return p_af_tune;
    }

    QString get_comment(void) const {
        return p_comment;
    }

    QString get_focus_distance(void) {
        return p_focus_distance;
    }

    QString get_focal_length(void) {
        return p_focal_length;
    }

    QString get_aperture(void) {
        return p_aperture;
    }
    
    QString get_focal_ratio(void) {
        return p_focal_ratio;
    }

	void set_exiv2_binary(const QString& s) {
		exiv2_binary = s;
	}

  private:
    QString extract_property(QString propname);
    QString extract_af_tune(void);
    QString extract_comment(void);
    QString extract_focus_distance(void);
    QString extract_focal_length(void);
    QString extract_aperture(void);
    QString extract_focal_ratio(void);


    char*   eat_whitespace(char* p);
    char*   eat_non_whitespace(char* p);

    QString exiv2_binary;

    QString ifname;
    QString tfname;

    QString p_af_tune;
    QString p_comment;
    QString p_focus_distance;
    QString p_focal_length;
    QString p_aperture;
    QString p_focal_ratio;
};

#endif

