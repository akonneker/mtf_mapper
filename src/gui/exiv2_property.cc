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
#include "exiv2_property.h"
#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

using std::vector;
using std::string;
using std::ifstream;
using std::cout;
using std::endl;

#include <QProcess>

#if QT_VERSION >= QT_VERSION_CHECK(5,14,0)
    #define EMPTY_STRING_PARTS Qt::SkipEmptyParts
#else
    #define EMPTY_STRING_PARTS QString::SkipEmptyParts
#endif

Exiv2_property::Exiv2_property(QString bin_name, QString ifname, QString tfname)
: exiv2_binary(bin_name),ifname(ifname), tfname(tfname)
{

    // ping exiv2 to get camera make
    QString make = extract_property(QString("Exif.Image.Make"));

    if (make.startsWith(QString("NIKON"), Qt::CaseInsensitive)) mode = NIKON;
    else if (make.startsWith(QString("Canon"), Qt::CaseInsensitive)) mode = CANON;
    else if (make.startsWith(QString("OLYMPUS"), Qt::CaseInsensitive)) mode = OLYMPUS;
    else if (make.length() == 0) mode = NONE;
    else mode = OTHER;

    p_af_tune = extract_af_tune();
    p_comment = extract_comment();
    p_focus_distance = extract_focus_distance();
    p_focal_length = extract_focal_length();
    p_aperture = extract_aperture();
    p_focal_ratio = extract_focal_ratio();
}

char*   Exiv2_property::eat_non_whitespace(char* cp) {
    while (*cp != ' ' && *cp != '\t') {
        cp++;
    }
    return cp;
}

char*   Exiv2_property::eat_whitespace(char* cp) {
    while (*cp == ' ' || *cp == '\t') {
        cp++;
    }
    return cp;
}

QString Exiv2_property::extract_property(QString propname) {
    vector<char> buffer(4096);

    QProcess evp;
    evp.setProgram(exiv2_binary);
    evp.setStandardOutputFile(tfname);
    evp.setArguments(
        QStringList() << "-g" << propname << "-pt" << ifname
    );
    evp.start();
    evp.waitForFinished(-1);
    int rval2 = evp.exitStatus() == QProcess::NormalExit;

    if (!rval2) {
        logger.error("exiv2 call failed on file %s : prop %s, exitstatus=%d, exitcode=%d\n", 
            ifname.toLocal8Bit().constData(), propname.toLocal8Bit().constData(),
            evp.exitStatus(), evp.exitCode()
        );
        return QString("");
    }

    ifstream ifs(tfname.toLocal8Bit().constData());
    if (!ifs.fail()) {
        ifs.getline(buffer.data(), 4096);
        buffer[4095] = 0; // ensure strlen will stop

        if (strlen(buffer.data()) == 0) {
            return QString("N/A");
        }

        // seek over three whitespace regions
        char* cp = buffer.data();
        for (int i=0; i < 3; i++) {
            cp = eat_non_whitespace(cp);
            cp = eat_whitespace(cp);
        }

        QString rval = QString(cp);
        // now cp points to the last record, which is the one we want
        return rval;

    } else {
        logger.error("failed to open %s [input=%s]\n", tfname.toLocal8Bit().constData(), ifname.toLocal8Bit().constData());
    }

    return QString("not found");
}


QString Exiv2_property::extract_af_tune(void) {
    if (mode != NONE) {
        switch (mode) {
        case NIKON:
            {
                QString aft = extract_property(QString("Exif.NikonAFT.AFFineTuneAdj"));
                bool ok;
                int value = aft.toInt(&ok);
                if (value > 20) {
                    value -= 256;
                }
                return QString("%1").arg(value);
            }
            
            break;
        case CANON:
            return extract_property(QString("Exif.NikonAFT.AFFineTuneAdj"));
            break;
        case OTHER:
        case NONE:
        default:
            return QString("N/A");
        }
    } else {
        return QString("N/A");
    }
}

QString Exiv2_property::extract_comment(void) {
    if (mode != NONE) {
        return extract_property(QString("Exif.Photo.UserComment"));
    } else {
        return QString("N/A");
    }
}

QString Exiv2_property::extract_focus_distance(void) {
    if (mode != NONE) {
        switch (mode) {
        case NIKON:
            return extract_property(QString("Exif.NikonLd3.FocusDistance"));
            break;
        case CANON:
            return extract_property(QString("Exif.CanonSi.SubjectDistance")); // right tag?
            break;
        case OTHER:
        case NONE:
        default:
            return QString("N/A");
        }
    } else {
        return QString("N/A");
    }
}

QString Exiv2_property::extract_focal_length(void) {
    if (mode != NONE) {
        return extract_property(QString("Exif.Photo.FocalLength"));
    } else {
        return QString("N/A");
    }
}

QString Exiv2_property::extract_aperture(void) {
    if (mode != NONE) {
        return extract_property(QString("Exif.Photo.FNumber"));
    } else {
        return QString("N/A");
    }
}

QString Exiv2_property::extract_focal_ratio(void) {

    // try reading 35-mm equivalent focal length
    // this works for many cameras, including:    
    // Nikon D7000, D7100, 
    // Sony ILCA-77M2, ILCE-6300, ILCA-99M2
    // Panasonic DMC-GH4 
    // Pentax K-3 II, K-1
    // Fujifilm X-Pro2
    // f = Exif.Photo.FocalLengthIn35mmFilm
    // then fr = f/36
    
    QString f35 = extract_property(QString("Exif.Photo.FocalLengthIn35mmFilm"));
    
    if (f35.compare("N/A") != 0) {
        QStringList fparts = f35.split(' ', EMPTY_STRING_PARTS);
        if (fparts.size() < 1) {
            logger.debug("%s\n", "failed to parse 35 mm equivalent focal length field in EXIF data");
            return QString("-1");
        }
        double f = fparts[0].toDouble();
        
        // w=35.8 (canon 5d)
        // w=35.9 (pentax k-1)
        // w=36   (nikon D3, D5), or 35.97, 35.63 ?
        // just round up, call w=36 mm
        
        return f > 0 ? QString("%1").arg(f/36.0) : QString("-1");
    }
    
    // Olympus E-M5MarkII, E-M1, probably all Olympus ILCs?
    // a = Exif.OlympusIp.AspectRatio
    // f = Exif.Photo.FocalLength
    // d = Exif.OlympusEq.FocalPlaneDiagonal
    // w = sqrt(d*d/(1 + 1/(a*a)))
    // fr = f/w
    if (mode == OLYMPUS) {
        QString aspect = extract_property(QString("Exif.OlympusIp.AspectRatio"));
        QString focal = extract_property(QString("Exif.Photo.FocalLength"));
        QString diag = extract_property(QString("Exif.OlympusEq.FocalPlaneDiagonal"));
        
        if (aspect.compare("N/A") != 0 && focal.compare("N/A") != 0 && diag.compare("N/A") != 0) {
            QStringList arparts = aspect.split(':', EMPTY_STRING_PARTS);
            if (arparts.size() != 2) {
               return QString("-1");
            }
            double a = arparts[0].toDouble() / arparts[1].toDouble();
            
            QStringList dparts = diag.split('/', EMPTY_STRING_PARTS);
            if (dparts.size() != 2) {
               return QString("-1");
            }
            double d = dparts[0].toDouble() / dparts[1].toDouble();
            
            QStringList fparts = focal.split(' ', EMPTY_STRING_PARTS);
            if (fparts.size() < 1) {
                logger.debug("%s\n", "failed to parse  focal length field in EXIF data");
                return QString("-1");
            }
            
            double f = fparts[0].toDouble();
            double w = sqrt(d*d / (1 + 1.0/(a*a)));
            
            return QString("%1").arg(f/w);
        } else {
            return QString("-1");
        }
    }
    
    // Try the Canon method on everything else
    
    // Tested with Canon 70D, but probably works on all Canons
    // p = Exif.Photo.PixelXDimension
    // Exif.Photo.FocalPlaneResolutionUnit (just a sanity check, should be inch?)
    // r = Exif.Photo.FocalPlaneXResolution
    // f = Exif.Photo.FocalLength
    // then fr = f/(25.4*p/r)
    
    QString xres = extract_property(QString("Exif.Photo.FocalPlaneXResolution"));
    if (xres.compare("N/A") != 0) {
        double r = xres.toDouble();
        
        QString pix_width = extract_property(QString("Exif.Photo.PixelXDimension"));
        double p = pix_width.toDouble();
        
        QString focal = extract_property(QString("Exif.Photo.FocalLength"));
        
        QStringList fparts = focal.split(' ', EMPTY_STRING_PARTS);
        if (fparts.size() < 1) {
            logger.debug("%s\n", "failed to parse  focal length field in EXIF data");
            return QString("-1");
        }
        
        double f = fparts[0].toDouble();
        
        if (f > 0 && r > 0 && p > 0) {
            return QString("%1").arg(f/(25.4*p/r));
        } else {
            return QString("-1");
        }
    }
    
    return QString("-1");
}
