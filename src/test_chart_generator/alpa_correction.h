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

#ifndef ALPA_CORRECTION_H
#define ALPA_CORRECTION_H

vector< pair<double, double> > alpa_correction_table = {
    {18,197.4524975785},
    {17,185.8135199602},
    {16,174.0755473925},
    {15,162.5355647235},
    {14,151.0804348683},
    {13,139.7242999625},
    {12,128.4530178703},
    {11,117.2807307276},
    {10,106.1932963986},
    {9,95.1907148833},
    {8,84.2729861818},
    {7,73.440110294},
    {6,62.7062293556},
    {5,52.0430590953},
    {4,41.4505995132},
    {3,30.9854191516},
    {2,20.5768073325},
    {1,10.2389061916},
    {0,0},
    {1,10.0974848353},
    {2,20.2656803488},
    {3,30.292454506},
    {4,40.2343758495},
    {5,50.1055865149},
    {6,59.8919443665},
    {7,69.60759154},
    {8,79.266670171},
    {9,88.8408959883},
    {10,98.3444111274},
    {11,107.7772155885},
    {12,117.1393093714},
    {13,126.4306924762},
    {14,135.6655070385},
    {15,144.8578951939},
    {16,153.9230041287},
    {17,162.9456866566},
    {18,171.9118006421},
    {19,180.821346085},
    {20,189.6460387142}
};

#endif
