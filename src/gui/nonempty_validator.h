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
#ifndef NONEMPTY_VALIDATOR_H
#define NONEMPTY_VALIDATOR_H

#include <QDoubleValidator>
#include <QIntValidator>

class Nonempty_DoubleValidator : public QDoubleValidator {
  public:
    explicit Nonempty_DoubleValidator(QObject* parent= nullptr) :QDoubleValidator(parent){}
    
    Nonempty_DoubleValidator(double bottom, double top, int decimals, double fallback, QObject* parent=nullptr)
        :QDoubleValidator(bottom, top, decimals, parent), fallback(fallback) {}

    virtual void fixup(QString& input)const override {
        input = QString::number(fallback, 'f', decimals());
    }
    
    double fallback = 0;
};

class Nonempty_IntValidator : public QIntValidator {
  public:
    explicit Nonempty_IntValidator(QObject* parent= nullptr) :QIntValidator(parent){}
    
    Nonempty_IntValidator(double bottom, double top, double fallback, QObject* parent=nullptr)
        :QIntValidator(bottom, top, parent), fallback(fallback) {}

    virtual void fixup(QString& input)const override {
        input = QString::number(fallback);
    }
    
    int fallback = 0;
};

#endif

