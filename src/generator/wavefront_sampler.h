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
#ifndef WAVEFRONT_SAMPLER_H
#define WAVEFRONT_SAMPLER_H

#include "include/common_types.h"
#include "../../include/threadpool.h"

#include <cmath>
#include <complex>
using std::complex;

class Wavefront_sampler {
  public:
    Wavefront_sampler(void) {};
    Wavefront_sampler(double lambda, double pitch, double N, double w020=0, double w040=0) 
      : lambda(lambda), pitch(pitch), N(N), w020(w020), w040(w040),
        nsamples(6000), x(nsamples, 0), y(nsamples, 0), cached_psf(nsamples) {
        
        calculate_bessel_roots();
        
        // first point is replicated to allow interpolation later
        x[0] = 0;
        y[0] = 0;
        // sample densely near the center of the PSF
        size_t hr_samples = nsamples/3;
        double hr_diam = 6.0;
        double step = hr_diam/double(hr_samples-1);
        size_t n;
        double cx = 0;
        for (n=1; n < hr_samples; n++) {
            x[n] = cx;
            cx += step;
        }
        // and a bit more sparsely far out
        step = (diam - hr_diam)/double(nsamples - hr_samples - 1);
        for (; n < nsamples; n++) {
            x[n] = cx;
            cx += step;
        }
        
        printf("starting PSF render ... ");
        ThreadPool& tp = ThreadPool::instance();
        
        cached_psf[0] = 0;
        
        size_t nblocks = tp.size();
        size_t blocksize = cached_psf.size() / nblocks;
        vector< std::future<void> > futures;
        for (size_t b=0; b < nblocks; b++) {
            size_t lower = std::max(size_t(1), b * blocksize);
            size_t upper = std::min((b+1)*blocksize, cached_psf.size());
            
            futures.emplace_back( 
                tp.enqueue( [lower,upper,this] {
                    for (size_t k=lower; k < upper; k++) {
                        cached_psf[k] = psf(x[k] + 0.5*(x[k] - x[k-1]));
                    }
                })
            );
            
        }
        for (size_t i=0; i < futures.size(); i++) {
            futures[i].wait();
        }
        printf("done!\n");
        
        // Build cumulative distribution by integrating psf(r)*r piecewise over
        // each interval. We approximate psf(r) as a linear function over the
        // interval.
        y[1] = cached_psf[1]*x[1]*x[1];
        for (size_t k=2; k < cached_psf.size(); k++) {
            double m = (cached_psf[k] - cached_psf[k-1]) / (x[k] - x[k-1]);
            double c = cached_psf[k-1] - m*x[k-1];
            double I = m/3.0*(x[k]*x[k]*x[k] - x[k-1]*x[k-1]*x[k-1]) + 0.5*c*(x[k]*x[k] - x[k-1]*x[k-1]);
            y[k] = y[k-1] + I;
        }
        
        for (size_t n=0; n < nsamples; n++) {
            y[n] /= y[nsamples-1];
        }
        
        diam_prob = 1; // TODO: can we do better than this?
    }
    
    double rairy2d(double& sx, double& sy, normal_sampler& sampler) {
        double px;
        double py;
        
        do {
            sampler.runif2d(px, py, 0.5, 0.5);
            px += 0.5; // undo the scaling of runif2d
            py += 0.5;
            py *= 2*M_PI;
        } while (px > diam_prob);
        
        size_t xidx = upper_bound(y.begin(), y.end(), px) - y.begin();
        double xp;
        // interpolate to prevent stratification
        if (xidx < nsamples - 1) {
            double alpha = (px - y[xidx])/(y[xidx+1] - y[xidx]);
            xp = x[xidx]*(1-alpha) + x[xidx+1]*alpha;
        } else {
            xp = x[xidx];
        }
        
        double rad = (xp) * (lambda/pitch)*N;
        sx = rad*cos(py);
        sy = rad*sin(py);
        
        if (xidx == 0) return 1e-6;
        return (y[xidx] - y[xidx-1]) / ((M_PI*x[xidx]*x[xidx] - M_PI*x[xidx-1]*x[xidx-1]));
    }
    
    double get_psf(double r) const {
        if (r > x.back()) {
            return cached_psf.back();
        }
        size_t idx = upper_bound(x.begin(), x.end(), r) - x.begin();
        if (idx < nsamples - 1) {
            double alpha = (r - x[idx])/(x[idx+1] - x[idx]);
            return std::max(0.0, cached_psf[idx]*(1 - alpha) + cached_psf[idx+1]*alpha);
        } else {
            return cached_psf[idx];
        }
    }
    
    static constexpr double diam = 120.0;
    double diam_prob;
    
    double lambda;
    double pitch;
    double N;
    double w020;
    double w040;
    
    size_t nsamples;
    vector<double> x;
    vector<double> y;
    
  private:
    // Hankle-transform integral from Wyant and Creath (J.C. Wyant, K. Creath, Basic Wavefront Aberration 
    // Theory for Optical Metrology, Applied Optics and Optical Engineering, 11(29), 1992), page 20.
    // Both w040 and w020 are specified as fractions/multiples of the wavelength.
	std::complex<double> f(double rho, double r) {
        double rhosq = rho * rho;
        std::complex<double> abb = std::exp(std::complex<double>(0, 2*M_PI*(w040*rhosq*rhosq + w020*rhosq)));
        double bes = j0(r*M_PI*rho) * rho;
        return abb * bes;
    }

    double cnorm(const std::complex<double>& c) {
        return std::abs(c)*std::abs(c);
    }
    
    // Adapted from TOMS468
    int quad_cp(double a, double b, std::complex<double> *result, int& k, 
        double epsil, int& npts, int& icheck, double r) {
        
        /* 
           Paraphrased from TOMS468:  
           The result is obtained using a sequence of 1, 3, 7, 15, 31, 63, 
           127, and 255 point interlacing formulae (no integrand 
           evaluations are wasted) of respective degrees 1, 5, 11, 23, 
           47, 95, 191 and 383.  the formulae are based on the optimal
           extensions of the 3-point gauss formula.  details of
           the formulae are given in 'the optimum addition of points
           to quadrature formulae' by t.n.l. patterson, maths. comp.
           vol 22, 847-856, 1968.
        */
        
        // Gauss-Kronrod coefficients of different orders, flattened
        constexpr double p[381] = {
            0.77459666924148337704, 0.55555555555555555556, 0.88888888888888888889, 0.26848808986833344073, 0.96049126870802028342, 0.10465622602646726519, 0.434243749346802558, 0.40139741477596222291, 0.45091653865847414235,
            0.13441525524378422036, 0.051603282997079739697, 0.20062852938698902103, 0.99383196321275502221, 0.017001719629940260339, 0.88845923287225699889, 0.092927195315124537686, 0.62110294673722640294,
            0.17151190913639138079, 0.22338668642896688163, 0.2191568584015874964, 0.22551049979820668739, 0.06720775429599070354, 0.025807598096176653565, 0.10031427861179557877, 0.0084345657393211062463,
            0.046462893261757986541, 0.085755920049990351154, 0.10957842105592463824, 0.99909812496766759766, 0.0025447807915618744154, 0.98153114955374010687, 0.016446049854387810934, 0.92965485742974005667,
            0.035957103307129322097, 0.8367259381688687355, 0.056979509494123357412, 0.70249620649152707861, 0.076879620499003531043, 0.53131974364437562397, 0.093627109981264473617, 0.33113539325797683309,
            0.10566989358023480974, 0.11248894313318662575, 0.11195687302095345688, 0.11275525672076869161, 0.033603877148207730542, 0.012903800100351265626, 0.050157139305899537414, 0.0042176304415588548391,
            0.023231446639910269443, 0.042877960025007734493, 0.054789210527962865032, 0.0012651565562300680114, 0.0082230079572359296693, 0.017978551568128270333, 0.028489754745833548613, 0.038439810249455532039,
            0.046813554990628012403, 0.052834946790116519862, 0.055978436510476319408, 0.99987288812035761194, 3.6322148184553065969e-4, 0.99720625937222195908, 0.0025790497946856882724, 0.98868475754742947994,
            0.0061155068221172463397, 0.97218287474858179658, 0.010498246909621321898, 0.94634285837340290515, 0.015406750466559497802, 0.9103711569570042925, 0.020594233915912711149, 0.86390793819369047715,
            0.025869679327214746911, 0.80694053195021761186, 0.03107355111168796488, 0.73975604435269475868, 0.03606443278078257264, 0.66290966002478059546, 0.040714410116944318934, 0.57719571005204581484,
            0.044914531653632197414, 0.48361802694584102756, 0.048564330406673198716, 0.38335932419873034692, 0.051583253952048458777, 0.27774982202182431507, 0.053905499335266063927, 0.16823525155220746498,
            0.055481404356559363988, 0.056344313046592789972, 0.056277699831254301273, 0.056377628360384717388, 0.016801938574103865271, 0.0064519000501757369228, 0.025078569652949768707, 0.0021088152457266328793,
            0.011615723319955134727, 0.021438980012503867246, 0.027394605263981432516, 6.3260731936263354422e-4, 0.0041115039786546930472, 0.0089892757840641357233, 0.014244877372916774306, 0.019219905124727766019,
            0.023406777495314006201, 0.026417473395058259931, 0.027989218255238159704, 1.8073956444538835782e-4, 0.0012895240826104173921, 0.0030577534101755311361, 0.0052491234548088591251, 0.0077033752332797418482,
            0.010297116957956355524, 0.012934839663607373455, 0.01553677555584398244, 0.01803221639039128632, 0.020357755058472159467, 0.022457265826816098707, 0.024282165203336599358, 0.025791626976024229388,
            0.026952749667633031963, 0.027740702178279681994, 0.028138849915627150636, 0.99998243035489159858, 5.0536095207862517625e-5, 0.99959879967191068325, 3.7774664632698466027e-4, 0.99831663531840739253,
            9.3836984854238150079e-4, 0.99572410469840718851, 0.0016811428654214699063, 0.9914957211781061324, 0.0025687649437940203731, 0.98537149959852037111, 0.0035728927835172996494, 0.97714151463970571416,
            0.0046710503721143217474, 0.96663785155841656709, 0.0058434498758356395076, 0.95373000642576113641, 0.007072489995433555468, 0.93832039777959288365, 0.0083428387539681577056, 0.92034002547001242073,
            0.0096411777297025366953, 0.89974489977694003664, 0.010955733387837901648, 0.87651341448470526974, 0.012275830560082770087, 0.85064449476835027976, 0.01359157100976554679, 0.82215625436498040737,
            0.014893641664815182035, 0.79108493379984836143, 0.016173218729577719942, 0.75748396638051363793, 0.017421930159464173747, 0.72142308537009891548, 0.018631848256138790186, 0.68298743109107922809,
            0.019795495048097499488, 0.64227664250975951377, 0.020905851445812023852, 0.59940393024224289297, 0.021956366305317824939, 0.55449513263193254887, 0.022940964229387748761, 0.50768775753371660215,
            0.02385405210603854008, 0.45913001198983233287, 0.024690524744487676909, 0.40897982122988867241, 0.025445769965464765813, 0.35740383783153215238, 0.02611567337670609768, 0.30457644155671404334,
            0.026696622927450359906, 0.25067873030348317661, 0.027185513229624791819, 0.19589750271110015392, 0.027579749566481873035, 0.14042423315256017459, 0.027877251476613701609, 0.08445404008371088371,
            0.028076455793817246607, 0.028184648949745694339, 0.028176319033016602131, 0.028188814180192358694, 0.0084009692870519326354, 0.0032259500250878684614, 0.012539284826474884353, 0.0010544076228633167722,
            0.0058078616599775673635, 0.010719490006251933623, 0.013697302631990716258, 3.1630366082222647689e-4, 0.0020557519893273465236, 0.0044946378920320678616, 0.0071224386864583871532, 0.0096099525623638830097,
            0.011703388747657003101, 0.013208736697529129966, 0.013994609127619079852, 9.0372734658751149261e-5, 6.4476204130572477933e-4, 0.0015288767050877655684, 0.0026245617274044295626, 0.0038516876166398709241,
            0.0051485584789781777618, 0.0064674198318036867274, 0.00776838777792199122, 0.00901610819519564316, 0.010178877529236079733, 0.011228632913408049354, 0.012141082601668299679, 0.012895813488012114694,
            0.013476374833816515982, 0.013870351089139840997, 0.014069424957813575318, 2.5157870384280661489e-5, 1.8887326450650491366e-4, 4.6918492424785040975e-4, 8.4057143271072246365e-4, 0.0012843824718970101768,
            0.0017864463917586498247, 0.0023355251860571608737, 0.0029217249379178197538, 0.003536244997716777734, 0.0041714193769840788528, 0.0048205888648512683476, 0.005477866693918950824, 0.0061379152800413850435,
            0.0067957855048827733948, 0.0074468208324075910174, 0.008086609364788859971, 0.0087109650797320868736, 0.0093159241280693950932, 0.009897747524048749744, 0.010452925722906011926, 0.01097818315265891247,
            0.01147048211469387438, 0.01192702605301927004, 0.012345262372243838455, 0.012722884982732382906, 0.01305783668835304884, 0.013348311463725179953, 0.01359275661481239591, 0.013789874783240936517,
            0.013938625738306850804, 0.014038227896908623303, 0.014088159516508301065, 0.99999759637974846476, 6.937936432410326717e-6, 0.99994399620705437576, 5.3275293669780613125e-5, 0.99976049092443204733,
            1.3575491094922871973e-4, 0.99938033802502358193, 2.4921240048299729402e-4, 0.9987456144680951147, 3.8974528447328229322e-4, 0.99780535449595727456, 5.5429531493037471492e-4, 0.99651414591489027385,
            7.4028280424450333046e-4, 0.99483150280062100052, 9.4536151685852538246e-4, 0.99272134428278861533, 0.0011674841174299594077, 0.99015137040077015918, 0.0014049079956551446427, 0.98709252795403406719,
            0.0016561127281544526052, 0.98351865757863272876, 0.0019197129710138724125, 0.97940628167086268381, 0.0021944069253638388388, 0.97473445975240266776, 0.0024789582266575679307, 0.96948465950245923177,
            0.002772195764593450994, 0.96364062156981213252, 0.0030730184347025783234, 0.95718821610986096274, 0.0033803979910869203823, 0.95011529752129487656, 0.0036933779170256508183, 0.94241156519108305981,
            0.0040110687240750233989, 0.934068436157725788, 0.0043326409680929828545, 0.92507893290707565236, 0.0046573172997568547773, 0.91543758715576504064, 0.0049843645647655386012, 0.90514035881326159519,
            0.0053130866051870565663, 0.89418456833555902286, 0.0056428181013844441585, 0.88256884024734190684, 0.0059729195655081658049, 0.87029305554811390585, 0.0063027734490857587172, 0.85735831088623215653,
            0.0066317812429018878941, 0.84376688267270860104, 0.0069593614093904229394, 0.82952219463740140018, 0.0072849479805538070639, 0.81462878765513741344, 0.0076079896657190565832, 0.7990922909608414018,
            0.0079279493342948491103, 0.78291939411828301639, 0.0082443037630328680306, 0.76611781930376009072, 0.0085565435613076896192, 0.74869629361693660282, 0.0088641732094824942641, 0.73066452124218126133,
            0.0091667111635607884067, 0.71203315536225203459, 0.0094636899938300652943, 0.69281376977911470289, 0.0097546565363174114611, 0.6730188302304184792, 0.010039172044056840798, 0.6526616654100174961,
            0.010316812330947621682, 0.63175643771119423041, 0.010587167904885197931, 0.61031811371518640016, 0.010849844089337314099, 0.58836243444766254143, 0.011104461134006926537, 0.56590588542365442262,
            0.011350654315980596602, 0.54296566649831149049, 0.011588074033043952568, 0.51955966153745702199, 0.011816385890830235763, 0.49570640791876146017, 0.01203527078527956263, 0.47142506587165887693,
            0.012244424981611985899, 0.44673539866202847374, 0.012443560190714035263, 0.42165768662616330006, 0.012632403643542078765, 0.39621280605761593918, 0.012810698163877361967, 0.37042208795007823014,
            0.012978202239537399286, 0.34430734159943802278, 0.013134690091960152836, 0.31789081206847668318, 0.01327995174393053065, 0.29119514851824668196, 0.013413793085110098513, 0.26424337241092676194,
            0.013536035934956213614, 0.23705884558982972721, 0.013646518102571291428, 0.20966523824318119477, 0.013745093443001896632, 0.18208649675925219825, 0.013831631909506428676, 0.15434681148137810869,
            0.013906019601325461264, 0.12647058437230196685, 0.013968158806516938516, 0.09848239659811920209, 0.01401796803945660881, 0.070406976042855179063, 0.014055382072649964277, 0.042269164765363603212,
            0.014080351962553661325, 0.014093886410782462614, 0.014092845069160408355, 0.014094407090096179347
        };

        int i;
        double x, sum, diff;
        std::complex<double> acum;
        int iold, inew;
        std::complex<double> funct[127];
        std::complex<double> fzero;

        --result; // compensate for 1-based array. eeew.

        icheck = 0;
        // check for trivial case
        if (a == b) {
            k = 2;
            result[1] = 0;
            result[2] = 0;
            npts = 0;
            return 0;
        }
        
        // scale factors
        sum = (b + a) / 2.0;
        diff = (b - a) / 2.0;
        
        bool initialize = true;
        
        while (true) {
            if (initialize) {
                //  1-point Gauss
                fzero = f(sum, r);
                result[1] = fzero * 2.0 * diff;
                i = 0;
                iold = 0;
                inew = 1;
                k = 2;
                acum = 0;
                initialize = false;
            } else {
                if (k == 8) { // maximum order reached before convergence
                    icheck = 1;
                    npts = inew + iold;
                    return 0;
                }
                
                ++k;
                acum = 0;
                // contribution from function values already computed
                for (int j = 1; j <= iold; ++j) {
                    ++i;
                    acum += p[i - 1] * funct[j - 1];
                }
            }
            // contribution from new function values
            iold += inew;
            for (int j = inew; j <= iold; ++j) {
                ++i;
                x = p[i - 1] * diff;
                funct[j - 1] = f(sum + x, r) + f(sum - x, r);
                ++i;
                acum += p[i - 1] * funct[j - 1];
            }
            
            inew = iold + 1;
            ++i;
            result[k] = (acum + p[i - 1] * fzero) * diff;
            
            // check for convergence
            if (cnorm(result[k] - result[k - 1]) - epsil * cnorm(result[k]) <= 0.0) {
                npts = inew + iold;
                return 0;
            } 
        }
    } 
    
    int i_sign(int a, int b) {
        if (b < 0) {
            return -abs(a);
        } else {
            return abs(a);
        }
    }
    
    // Adapted from TOMS468
    std::complex<double> qsuba_cp(double a, double b, double epsil, int& npts, int& icheck, 
        double& relerr, double r) {
        
        constexpr int ismax = 100;

        /* System generated locals */
        std::complex<double> ret_val(0, 0);

        int k, ic, nf, is;
        double sub1, sub2, sub3;
        double comp, stack[ismax], estim;
        std::complex<double> result[8];

        quad_cp(a, b, result, k, epsil, npts, icheck, r);
        ret_val = result[k - 1];
        relerr = 0;
        if (cnorm(ret_val) != 0) {
            relerr = cnorm(result[k - 1] - result[k - 2]) / cnorm(ret_val);
        }
        // check if subdivision is needed
        if (icheck == 0) {
            return ret_val;
        }
        
        // subdivide
        estim = fabs(cnorm(ret_val) * epsil);
        relerr = 0;
        ret_val = std::complex<double>(0, 0);
        is = 1;
        ic = 1;
        sub1 = a;
        sub3 = b;
        
        while (true) { 
            sub2 = 0.5*(sub1 + sub3);
            quad_cp(sub1, sub2, result, k, epsil, nf, icheck, r);
            npts += nf;
            comp = cnorm(result[k - 1] - result[k - 2]);
            
            bool accumulate = true;
            if (icheck == 1) {
                if (comp <= estim) {
                    ic = i_sign(2, ic);
                } else {
                    // subdivide
                    if (is < ismax) {
                        stack[is - 1] = sub1;
                        ++is;
                        stack[is - 1] = sub2;
                        ++is;
                        accumulate = false;
                    } else {
                        printf("qsuba split depth reached\n");
                        ic = -abs(ic); 
                    }
                }
            } 
            if (accumulate) {
                ret_val += result[k - 1];
                relerr += comp;
            }
            
            quad_cp(sub2, sub3, result, k, epsil, nf, icheck, r);
            npts += nf;
            comp = cnorm(result[k - 1] - result[k - 2]);
            
            if (icheck == 1) {
                if (comp <= estim) {
                    ic = i_sign(2, ic);
                } else {
                    sub1 = sub2;
                    continue;
                }
            }
            
            ret_val += result[k - 1];
            relerr += comp;
            
            if (is == 1) {
                //  subdivision result 
                icheck = ic;
                relerr /= cnorm(ret_val);
                return ret_val;
            }
            
            // pop the last interval that failed to converge
            --is;
            sub3 = stack[is - 1];
            --is;
            sub1 = stack[is - 1];
        }
    }
	
	inline double psf(double r) {
        double a = 0.0;
        double b = 1.0;
        double epsil = 1e-5;
        int npts = 0;
        int icheck = 0;
        double relerr = 0;
	    double p = 0;
		std::complex<double> cp(0, 0);
		
		// with r < 1 we should not have any roots
		if (r < 1) {
		    cp = qsuba_cp(a, b, epsil, npts, icheck, relerr, r);
		    return std::abs(cp)*std::abs(cp);
		}
		
		size_t expected_roots = ceil(r) + 1;
		size_t ic_sum = 0;
		size_t ri = 0;
		bool done = false;
		// otherwise integrate over all the intervals defined by the roots
		for (; !done && ri < expected_roots; ri++) {
		    b = jr[ri] / (r*M_PI);
		    if (b > 1.0) { 
		        b = 1.0;
		        done = true;
            }
            if (a < 1.0) {
                cp += qsuba_cp(a, b, epsil, npts, icheck, relerr, r);
                ic_sum += icheck;
            }
		    
		    a = b;
        }
        p = std::abs(cp)*std::abs(cp);
        
        if (ic_sum > 0) {
            printf("failed to converge (%d) with r=%.8lf\n", int(ic_sum), r);
        }
        
		return p;
    }
    
    void calculate_bessel_roots(void) {
		
		jr = vector<double>(ceil(Wavefront_sampler::diam)+1);
		
		jr[0] = 2.404826;
		jr[1] = 5.520078;
		
		constexpr double reps = 1e-12;
		for (size_t k=2; k < jr.size(); k++) {
		    double er = jr[k-1] + (jr[k-1] - jr[k-2]);
		    double x0 = er;
		    double x1 = x0 +  j0(x0)/j1(x0);
		    size_t it = 0;
		    while (fabs(x1 - x0) > fabs(x1)*reps && it < 5000) {
		        x0 = x1;
		        x1 = x0 + j0(x0)/j1(x0);
		        it++;
		    }
		    jr[k] = x1;
		}
    }
    
    vector<double> jr; // roots of Bessel j0 that fall within diam
	vector<double> cached_psf; // PSF evaluated at x[i]
};

#endif

