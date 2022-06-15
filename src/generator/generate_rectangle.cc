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
#include <math.h>
#include <stdlib.h>
#include "render.h"
#include "render_integral.h"
#include "render_importance_sampling.h"
#include "noise_source.h"
#include "multipolygon_geom.h" // temporary?
#include "quadtree.h"
#include "png_icc_headers.h"
#include "../../include/stride_range.h"
#include "../../include/threadpool.h"
#include <string.h>

#include "config.h"

#include "tclap/CmdLine.h"
#include "tclap_constraints.h"

#include <vector>
using std::vector;

#include <iostream>
#include <string>

using std::string;
using std::stringstream;

using std::cout;
using std::endl;

#define SQR(x) ((x)*(x))

const int border = 30;
int dummy_crc = 0;

inline unsigned char reverse_gamma(double x) {
    const double C_linear = 0.0031308;
    const double S_linear = 12.9232102;
    const double SRGB_a = 0.055;
    
    if (x < C_linear) {
        return lrint(255 * x * S_linear);
    }
    return lrint(255 * ((1 + SRGB_a) * pow(x, 1.0/2.4) - SRGB_a));
}

inline double fwd_gamma(double x) { // x in [0,1]
    const double S_linear = 12.9232102;
    const double C_srgb = 0.04045;
    const double SRGB_a = 0.055;
    
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    
    if (x < C_srgb) {
        return x / S_linear;
    }
    
    return  pow( (x + SRGB_a)/(1+SRGB_a), 2.4 );
}

uint32_t read_uint32(FILE* fin) {
    unsigned char b0 = fgetc(fin) & 0xff;
    unsigned char b1 = fgetc(fin) & 0xff;
    unsigned char b2 = fgetc(fin) & 0xff;
    unsigned char b3 = fgetc(fin) & 0xff;
    return (uint32_t(b0) << 24) | (uint32_t(b1) << 16) |
        (uint32_t(b2) << 8) | uint32_t(b3);
}

void inject_png_icc_profile(string fname, bool linear) {
    int icc_chunk_size = png_gray_gamma_icc_size;
    unsigned char* icc_chunk = png_gray_gamma_icc_chunk;
    if (linear) {
        icc_chunk_size = png_gray_linear_icc_size;
        icc_chunk = png_gray_linear_icc_chunk;
    } 

    FILE* fin = fopen(fname.c_str(), "rb+");
    
    unsigned char magic[4];
    if (fin && fread(magic, 1, 4, fin) == 4 && memcmp(magic, "\x89\x50\x4e\x47", 4) == 0) {
        fseek(fin, 4, SEEK_CUR); // skip over rest of PNG signature
        uint32_t chunk_size = read_uint32(fin);
        if (fread(magic, 1, 4, fin) != 4 || memcmp(magic, "IHDR", 4) != 0) return; 
        
        printf("Output file is PNG, injecting ICC profile\n");
        
        fseek(fin, chunk_size + 4, SEEK_CUR);
        size_t write_pos = ftell(fin);
        
        const size_t bsize = 64*1024;
        vector<unsigned char> buffer(bsize, 0xaa);
        fseek(fin, 0, SEEK_END);
        size_t move_read = 0; 
        size_t read_size = bsize;
        
        if (int64_t(ftell(fin)) - int64_t(bsize) < int64_t(write_pos)) {
            read_size = int64_t(ftell(fin)) - int64_t(write_pos);
            move_read = write_pos;
        } else {
            read_size = bsize;
            move_read = int64_t(ftell(fin)) - int64_t(bsize);
        }
        
        bool done = false;
        while (!done) {
            // now move everything up
            fseek(fin, move_read, SEEK_SET);
            
            size_t nread = fread(buffer.data(), 1, read_size, fin);
            fseek(fin, move_read + icc_chunk_size, SEEK_SET);
            
            size_t nwritten = fwrite(buffer.data(), 1, nread, fin);
            done = move_read <= write_pos;
            
            if (nwritten != nread) {
                fprintf(stderr, "Error: write failed during ICC injection. Did you run out of disk space?\n");
                fclose(fin);
                return;
            }
            
            if (int64_t(move_read) - int64_t(bsize) < int64_t(write_pos)) {
                read_size = int64_t(move_read) - int64_t(write_pos);
                move_read = write_pos;
            } else {
                read_size = bsize;
                move_read = int64_t(move_read) - int64_t(bsize);
            }
        }
        
        fseek(fin, write_pos, SEEK_SET);
        fwrite(icc_chunk, 1, icc_chunk_size, fin);
        
        fclose(fin);
    }    
}

// functor for Stride_range
class Render_rows {
  public:
    Render_rows(cv::Mat& in_img, const Render_polygon& in_r, Noise_source& noise_source, 
        double blacklevel=0.025, double whitelevel=0.975, bool gamma_correct=true, bool use_16bit=false,
        int buffer_border=30, int& crc=dummy_crc, double analogue_scale=1.0,
        int offset_x=0, int offset_y=0)
     : img(in_img), rect(in_r), noise_source(noise_source),
       gamma_correct(gamma_correct), blacklevel(blacklevel), whitelevel(whitelevel),
       use_16bit(use_16bit), buffer_border(buffer_border),
       crc(crc), analogue_scale(analogue_scale),
       offset_x(offset_x), offset_y(offset_y) {
     
        
    }


    inline void putpixel(int row, int col, double value) const {

        value = noise_source.sample(value, row*img.cols + col);

        if (use_16bit) {
            img.at<uint16_t>(row, col) = lrint(value*65535);
        } else {
            if (gamma_correct) {
                img.at<uchar>(row, col) = reverse_gamma(value);
            } else {
                img.at<uchar>(row, col) = lrint(value*255);
            }
        }
    }
     
    void operator()(const Stride_range& r) const {
        for (size_t row=r.begin(); row != r.end(); r.increment(row)) {

            if (int(row) < buffer_border || int(row) > img.rows - buffer_border) {
                for (int col = 0; col < img.cols; col++) {
                    putpixel(row, col, whitelevel);
                }
            } else {
                for (int col = 0; col < buffer_border; col++) {
                    putpixel(row, col, whitelevel);
                }
                for (int col = img.cols - buffer_border; col < img.cols; col++) {
                    putpixel(row, col, whitelevel);
    
                }
                for (int col = buffer_border; col < img.cols-buffer_border; col++) {
                    
                    double rval = 0;
                    rval = rect.evaluate(col*analogue_scale + offset_x, row*analogue_scale + offset_y, blacklevel, whitelevel);
                    
                    putpixel(row, col, rval);
                }
            }
            
            crc++; // this will cause a race condition, but that does not matter ...
            if (crc % (img.rows/10) == 0) {
                //printf("completed %d rows (%.2lf%%)\n", crc, double(crc*100)/double(img.rows));
                printf("...%d%%", (int)lrint(double(crc*100)/double(img.rows)));
                fflush(stdout);
            }
        }
    } 
     
    cv::Mat& img;
    const Render_polygon& rect;
    Noise_source& noise_source;
    bool gamma_correct;
    double blacklevel;
    double whitelevel;
    bool use_16bit;
    int buffer_border;
    int& crc;
    double analogue_scale;
    int offset_x;
    int offset_y;
};


// functor for Stride_range
class Render_esf {
  public:
    Render_esf(const Render_polygon& in_r, vector< pair<double, double> >& esf, double length, double theta,
        int oversampling_factor, double xoff=0, double yoff=0)
     : rect(in_r), length(length), oversampling_factor(oversampling_factor),
       p(rect.t_geom.cx+xoff, rect.t_geom.cy+yoff), 
       sample_pos(Render_esf::n_samples(length, oversampling_factor)),
       esf(esf) {
        Point2d cur(p);
        Point2d d(cos(-theta)/oversampling_factor, sin(-theta)/oversampling_factor);
        cur = cur + (0.25+0.125)*length*oversampling_factor*d;
        for (int i=0; i < (int)sample_pos.size(); i++) {
            sample_pos[i] = cur;
            cur = cur + d;
        }
        
    }
    
    static int n_samples(double length, int oversampling_factor) {
        return (int)floor(length*oversampling_factor*0.25);
    }
    
    void write(const string& profile_fname) {
        FILE* fout = fopen(profile_fname.c_str(), "wt");
        for (size_t i=0; i < esf.size(); i++) {
            fprintf(fout, "%.12le %.12le\n", esf[i].first, esf[i].second);
        }
        fclose(fout);
    }


    void operator()(const Stride_range& r) const {
    
        for (size_t row=r.begin(); row != r.end(); r.increment(row)) {
            double rval = 0;
            rval = rect.evaluate(sample_pos[row].x, sample_pos[row].y, 0.05, 0.95);
            esf[row].first = norm(sample_pos[row] - p) - length*0.5;
            esf[row].second = rval;
        }
    } 
     
    const Render_polygon& rect;
    double length;
    int    oversampling_factor;
    Point2d  p;
    vector< Point2d > sample_pos;
    vector< pair<double, double> >& esf;
};


//------------------------------------------------------------------------------
int main(int argc, char** argv) {

    int width  = 300;
    int height = 300;
    
    double sigma = 0.3;
    double theta;
    int rseed;
    
    stringstream ss;
    ss << genrectangle_VERSION_MAJOR << "." << genrectangle_VERSION_MINOR;
    
    tclap_zero_one<double> zero_one;
    tclap_ranged<int> bit_depth_limits(2, 16);
    tclap_ranged<double> offset_limits(-1, 1);
    tclap_ranged<double> wavelength_limits(0.2, 0.9);
    tclap_ranged<double> fillfactor_limits(0.01, 1.0);
    tclap_nonneg<int> nonneg_int;
    tclap_nonneg<double> nonneg_double;
    
    TCLAP::CmdLine cmd("Generate rectangles with known MTF50 values", ' ', ss.str());
    TCLAP::ValueArg<std::string> tc_out_name("o", "output", "Output file name", false, "rect.png", "filename", cmd);
    TCLAP::ValueArg<double> tc_theta("a", "angle", "Orientation angle (degrees)", false, 4.0, "angle (degrees)", cmd);
    TCLAP::ValueArg<int> tc_seed("s", "seed", "Noise random seed", false, int(time((time_t*)nullptr)), "seed", cmd);
    TCLAP::ValueArg<double> tc_noise("n", "noise", "Noise magnitude (linear standard deviation, range [0,1])", false, 0.01, &zero_one, cmd);
    TCLAP::ValueArg<double> tc_blur("b", "blur", "Blur magnitude (linear standard deviation, range [0.185, +inf))", false, 0.374781, "std. dev", cmd);
    TCLAP::ValueArg<double> tc_mtf("m", "mtf50", "Desired MTF50 value (range (0, 1.0] cycles/pixel)", false, 0.3, &zero_one, cmd);
    TCLAP::ValueArg<double> tc_cr("c", "contrast", "Contrast reduction [0,1]", false, 0.1, &zero_one, cmd);
    TCLAP::ValueArg<double> tc_dim("d", "dimension", "Dimension of the image, in pixels", false, 100, &nonneg_double, cmd);
    TCLAP::ValueArg<double> tc_yoff("y", "yoffset", "Subpixel y offset [-1,1]", false, 0, &offset_limits, cmd);
    TCLAP::ValueArg<double> tc_xoff("x", "xoffset", "Subpixel x offset [-1,1]", false, 0, &offset_limits, cmd);
    TCLAP::ValueArg<double> tc_ar("r", "aspect-ratio", "Aspect ratio of rectangle [0,1]", false, 1.0, &zero_one, cmd);
    TCLAP::ValueArg<double> tc_read_noise("", "read-noise", "Read noise magnitude (linear standard deviation, in electrons, range [0,+inf))", false, 3.0, &nonneg_double, cmd);
    TCLAP::ValueArg<double> tc_pattern_noise("", "pattern-noise", "Fixed pattern noise magnitude (linear fraction of signal, range [0,1])", false, 0.02, &zero_one, cmd);
    TCLAP::ValueArg<double> tc_adc_gain("", "adc-gain", "ADC gain (linear electrons-per-DN, range (0,+inf))", false, 1.5, &nonneg_double, cmd);
    TCLAP::ValueArg<double> tc_psf_theta("", "psf-theta", "Angle between PSF major axis and horizontal", false, 0, "degrees", cmd);
    TCLAP::ValueArg<double> tc_psf_ratio("", "psf-ratio", "PSF minor axis fraction of major axis", false, 1.0, &nonneg_double, cmd);
    TCLAP::ValueArg<int> tc_adc_depth("", "adc-depth", "ADC depth (number of bits, range [2,16])", false, 14, &bit_depth_limits, cmd);
    TCLAP::ValueArg<double> tc_aperture("", "aperture", "Aperture f-number [0,1000]", false, 8.0, &nonneg_double, cmd);
    TCLAP::ValueArg<double> tc_aperture2("", "aperture2", "Aperture f-number of second axis of rectangular aperture stop (if -p rect-box) [0,1000]", false, 8.0, &nonneg_double, cmd);
    TCLAP::ValueArg<double> tc_pitch("", "pixel-pitch", "Pixel pitch (size) [0,20] micron", false, 4.73, &nonneg_double, cmd);
    TCLAP::ValueArg<double> tc_lambda("", "lambda", "Light wavelentgth (affects diffraction) [0.2,0.9] micron", false, 0.55, &wavelength_limits, cmd);
    TCLAP::ValueArg<double> tc_olpf_split("", "olpf-offset", "OLPF beam splitter offset (in units of pixels)", false, 0.375, &nonneg_double, cmd);
    TCLAP::ValueArg<int> tc_samples("", "airy-samples", "Number of half-samples (n) per axis per pixel for Airy PSFs [actual #samples = (2n+1)^2]", false, 30, &nonneg_int, cmd);
    TCLAP::ValueArg<std::string> tc_target_name("", "target-poly", "Target polygon file name", false, "poly.txt", "filename", cmd);
    TCLAP::ValueArg<double> tc_fillfactor("", "fill-factor", "Fill-factor of photosite [0.01,1]", false, 1.0, &fillfactor_limits, cmd);
    TCLAP::ValueArg<double> tc_ascale("", "analogue-scale", "Analogue scaling factor", false, 1.0, &nonneg_double, cmd);
    TCLAP::ValueArg<int> tc_roi_start_row("", "roi-row", "Region Of Interest (ROI) starting row", false, 0, "pixels", cmd);
    TCLAP::ValueArg<int> tc_roi_start_col("", "roi-col", "Region Of Interest (ROI) starting column", false, 0, "pixels", cmd);
    TCLAP::ValueArg<int> tc_roi_width("", "roi-width", "Region Of Interest (ROI) width", false, 0, &nonneg_int, cmd);
    TCLAP::ValueArg<int> tc_roi_height("", "roi-height", "Region Of Interest (ROI) height", false, 0, &nonneg_int, cmd);
    TCLAP::ValueArg<std::string> tc_photosite_poly("", "photosite-poly", "Photosite polygon file name", false, "photo.txt", "filename", cmd);
    TCLAP::ValueArg<double> tc_w020("", "w020", "Defocus w020 magnitude (wavelengths)", false, 0.0, "wavelengths", cmd);
    TCLAP::ValueArg<double> tc_w040("", "w040", "Spherical aberration w040 magnitude (wavelengths)", false, 0.0, "wavelengths", cmd);
    TCLAP::ValueArg<double> tc_blacklevel("", "blacklevel", "Simulated black level [0.0, 1)", false, 0.025, &zero_one, cmd);
    TCLAP::ValueArg<double> tc_whitelevel("", "whitelevel", "Simulated white level (0.0, 1.0]", false, 0.975, &zero_one, cmd);
    
    vector<string> psf_names;
    psf_names.push_back("gaussian");
    psf_names.push_back("gaussian-sampled");
    psf_names.push_back("airy");
    psf_names.push_back("airy-box");
    psf_names.push_back("airy-4dot-olpf");
    psf_names.push_back("wavefront");
    psf_names.push_back("wavefront-box");
    psf_names.push_back("rect");
    psf_names.push_back("rect-box");
    TCLAP::ValuesConstraint<string> psf_constraints(psf_names);
    TCLAP::ValueArg<std::string> tc_psf("p", "psf-type", "Point Spread Function (PSF) type", false, "gaussian", &psf_constraints );
    cmd.add(tc_psf);

    vector<string> photosite_names;
    photosite_names.push_back("square");
    photosite_names.push_back("circle");
    photosite_names.push_back("rounded-square");
    TCLAP::ValuesConstraint<string> photosite_constraints(photosite_names);
    TCLAP::ValueArg<std::string> tc_photosite_geom("", "photosite-geom", "Photosite aperture geometry", false, "square", &photosite_constraints );
    cmd.add(tc_photosite_geom);

    TCLAP::SwitchArg tc_linear("l","linear","Generate output image with linear intensities (default is sRGB gamma corrected)", cmd, false);
    TCLAP::SwitchArg tc_16("","b16","Generate linear 16-bit image", cmd, false);
    TCLAP::SwitchArg tc_profile("","esf-only","Sample the ESF, rather than render an image (default filename \"profile.txt\")", cmd, false);
    
    cmd.parse(argc, argv);
 
    rseed = tc_seed.getValue();
    srand(rseed);
    
    theta = tc_theta.getValue() / 180.0 * M_PI;
    double psf_theta = tc_psf_theta.getValue() / 180.0 * M_PI;
    
    Render_polygon_is::Render_type psf_type;
    if ( tc_psf.getValue().compare("gaussian") == 0) {
        psf_type = Render_polygon::GAUSSIAN;
    }
    if ( tc_psf.getValue().compare("gaussian-sampled") == 0) {
        psf_type = Render_polygon::GAUSSIAN_SAMPLED;
    }
    if ( tc_psf.getValue().compare("airy") == 0) {
        psf_type = Render_polygon::AIRY;
    }
    if ( tc_psf.getValue().compare("airy-box") == 0) {
        psf_type = Render_polygon::AIRY_PLUS_BOX;
    }
    if ( tc_psf.getValue().compare("airy-4dot-olpf") == 0) {
        psf_type = Render_polygon::AIRY_PLUS_4DOT_OLPF;
    }
    if ( tc_psf.getValue().compare("wavefront") == 0) {
        psf_type = Render_polygon::WAVEFRONT;
    }
    if ( tc_psf.getValue().compare("wavefront-box") == 0) {
        psf_type = Render_polygon::WAVEFRONT_PLUS_BOX;
    }
    if ( tc_psf.getValue().compare("rect") == 0) {
        psf_type = Render_polygon::RECT;
    }
    if ( tc_psf.getValue().compare("rect-box") == 0) {
        psf_type = Render_polygon::RECT_PLUS_BOX;
    }
    
    // perform some sanity checking
    if ( (tc_mtf.isSet() || tc_blur.isSet()) && 
        !(psf_type == Render_polygon_is::GAUSSIAN || psf_type == Render_polygon_is::GAUSSIAN_SAMPLED) ) {
        
        printf("Conflicting options selected. You cannot specify \"-m\" or \"-b\" options if PSF type is non-Gaussian.");
        printf("Aborting.\n");
        return -1;
    }
    
    double rwidth = 0;
    double rheight = 0;
    
    rwidth = tc_dim.getValue();
    if (rwidth < 35) {
        fprintf(stderr, "Warning: longest rectangle edge less than 35 pixels --- this may"
                        " produce inaccurate results when used with MTF Mapper\n");
    }
    
    double ar = tc_ar.getValue();
    if (rwidth*ar < 1) {
        fprintf(stderr, "Warning: specified aspect ratio too small, fixing rectangle at 1 pixel wide\n");
        ar = 1.0/rwidth;
    }
    if (rwidth*ar < 25) {
        fprintf(stderr, "Warning: specified aspect ratio will produce a thin rectangle (< 25 pixels wide).\n"
                        "This rectangle may produce inaccurate results when used with MTF Mapper.\n");
    }
    if (ar > 1) {
        fprintf(stderr, "Warning: aspect ratio > 1 specified, clipping to 1.0\n");
        ar = 1.0;
    }
    
    rheight = rwidth*ar;

    double diag = sqrt(rwidth*rwidth + rheight*rheight);
    width = lrint(diag + 32) + 2*border;
    height = width += width % 2;
    
    if (tc_mtf.isSet() && tc_blur.isSet()) {
        printf("Warning: you can not specify both blur and mtf50 values; choosing mtf50 value, and proceeding ...\n");
    }
    
    double mtf = tc_mtf.getValue();
    if (tc_mtf.isSet()) {
        sigma = sqrt( log(0.5)/(-2*M_PI*M_PI*mtf*mtf) );
    } else {
        sigma = tc_blur.getValue();
        mtf   = sqrt( log(0.5)/(-2*M_PI*M_PI*sigma*sigma) );
    }
    
    if (sigma < 0.185) {
        printf("It does not make sense to set blur below 0.185; you are on your own ...\n");
    }
    
    double blacklevel = tc_blacklevel.getValue();
    double whitelevel = tc_whitelevel.getValue();
    if (tc_cr.isSet()) {
        if (tc_blacklevel.isSet() || tc_whitelevel.isSet()) {
            fprintf(stderr, "Error: -c / --contrast cannot be used together with --blacklevel or --whitelevel.\nAborting.\n");
            return 1;
        }
        blacklevel = tc_cr.getValue() / 2.0;
        whitelevel = 1.0 - blacklevel;
    }
    
    if (tc_samples.getValue() < 1) {
        printf("Error. The argument to --airy-samples must be positive");
        return -1;
    }
    if (tc_samples.getValue() > 100 && !tc_profile.getValue()) {
        printf("Note: You have specified a large number of samples per pixel in 2D rendering mode.\n");
        printf("      Rendering may be very slow, unless you have an amazingly fast machine.\n\n");
    }
    if (tc_samples.getValue() < 5) {
        printf("Note: You have specified a very small number of samples per pixel (--airy-samples).\n");
        printf("      Output is going to be very noisy. You have been warned ...\n\n");
    }
    if (tc_samples.isSet() && 
        psf_type != Render_polygon::AIRY &&
        psf_type != Render_polygon::AIRY_PLUS_BOX &&
        psf_type != Render_polygon::AIRY_PLUS_4DOT_OLPF &&
        psf_type != Render_polygon::WAVEFRONT &&
        psf_type != Render_polygon::WAVEFRONT_PLUS_BOX &&
        psf_type != Render_polygon::RECT &&
        psf_type != Render_polygon::RECT_PLUS_BOX) {
        
        printf("Warning: You have specified the number of Airy samples (--airy-samples), but you\n");
        printf("         are not rendering with an Airy-based PSF.\n");
    }
    
    bool use_gamma = !tc_linear.getValue();
    bool use_16bit = tc_16.getValue();
    
    if (use_gamma && use_16bit) {
        printf("Setting both gamma and 16-bit output not supported. Choosing to disable gamma, and enable 16-bit output\n");
        use_gamma = false;
    }
    
    bool use_sensor_model = false;
    if (tc_read_noise.isSet() || tc_pattern_noise.isSet() ||
        tc_adc_depth.isSet()  || tc_adc_depth.isSet() ) {
        use_sensor_model = true;
    }

    if (tc_noise.isSet() && use_sensor_model) {
        printf("Noise sigma set, but full sensor noise model parameters also specified.\n"
               "Ignoring noise sigma, and going with full sensor model instead.\n");
    }
    
    printf("output filename = %s, theta = %lg degrees, seed = %d,\n ", 
        tc_out_name.getValue().c_str(), theta/M_PI*180, rseed
    );
    if (psf_type >= Render_polygon::AIRY) {
        printf("\t aperture = f/%.1lg, pixel pitch = %.3lg, lambda = %.3lg, ",
             tc_aperture.getValue(), tc_pitch.getValue(), tc_lambda.getValue()
		);
		if (psf_type == Render_polygon::AIRY_PLUS_4DOT_OLPF) {
			printf("OLPF split = %.3f pixels", tc_olpf_split.getValue());
		}
    } else {
        printf("\t sigma = %lg (or mtf50 = %lg), ", 
            sigma, mtf
        );
    }

    
    if (use_sensor_model) {
        printf("\n\t full sensor noise model, with read noise = %.1lg electrons, fixed pattern noise fraction = %.3lg,\n"
               "\t adc gain = %.2lg e/DN, adc depth = %d bits\n",
               tc_read_noise.getValue(), tc_pattern_noise.getValue(), 
               tc_adc_gain.getValue(), tc_adc_depth.getValue()
        );
    } else {
        printf("\n\t additive Gaussian noise with sigma = %lg\n", tc_noise.getValue());
    }
    
    printf("\t output in sRGB gamma = %d, intensity range [%lg, %lg], 16-bit output:%d, dimension: %dx%d\n",
        use_gamma, blacklevel, whitelevel, use_16bit, width, height
    );

    Geometry* target_geom = new Polygon_geom(
        width*0.5 + tc_xoff.getValue(), 
        height*0.5 + tc_yoff.getValue(),
        rwidth*sqrt(2.0),
        rheight*sqrt(2.0),
        M_PI/2 - theta,
        4
    );
    

    if (tc_target_name.isSet()) {
        delete target_geom;
        target_geom = new Quadtree (
            tc_xoff.getValue(),
            tc_yoff.getValue(),
            tc_target_name.getValue(),
            tc_ascale.getValue()
        );
        
       //((Quadtree*)target_geom)->print_bounds(0);
       width = 40 + target_geom->bounds.max_x / tc_ascale.getValue();
       height = 40 + target_geom->bounds.max_y / tc_ascale.getValue();
       width += width % 8;
       height += height % 8;
       printf("setting image dimensions to: %d, %d\n", 
           static_cast<int>(lrint(width)), 
           static_cast<int>(lrint(height))
       );
    }

    Render_polygon default_target(
        *target_geom,
        sigma,
        sigma*tc_psf_ratio.getValue(),
        M_PI/2 - psf_theta
    );

    Geometry* photosite_geom = new Polygon_geom(
        0, 0,
        sqrt(2*tc_fillfactor.getValue()), 
        sqrt(2*tc_fillfactor.getValue()),
        0, 4
    );

    bool display_mtf_equation = false;
    if (tc_photosite_geom.getValue().compare("square") == 0) {
        // do nothing
        display_mtf_equation = true;
    }
    if (tc_photosite_geom.getValue().compare("circle") == 0) {
        double eff = tc_fillfactor.getValue() * (M_PI/4.0);
        photosite_geom = new Polygon_geom(
            0, 0,
            2*sqrt(eff/M_PI), 2*sqrt(eff/M_PI),
            0, 60
        );
    }
    if (tc_photosite_geom.getValue().compare("rounded-square") == 0) {
        
        // a hard-coded shape that looks a bit like a blend between a box and a circle ...
        const int points_per_side = 20;
        vector<cv::Vec2d> verts(4*(points_per_side-1));
        int oidx = 0;
        
        // build top row
        double x1 = 1.0/(points_per_side/2);
        double scale = 0.5/(sqrt(1-x1*x1)*(1-pow(fabs(x1),1.8))/7.0 + 1) * tc_fillfactor.getValue();
        double x = 1;
        for (int i=0; i < points_per_side/2; i++) {
            double y = sqrt(1-x*x)*(1-pow(fabs(x),1.8))/7.0 + 1; // arbitrary empirical function that "looks ok"
            verts[oidx][0] = scale*x;
            verts[oidx][1] = scale*y;
            oidx++;
            x -= 1.0/(points_per_side/2);
            
        }
        x = 1 - 1.0/(points_per_side/2);;
        for (int i=0; i < (points_per_side/2-1); i++) {
            double y = sqrt(1-x*x)*(1-pow(fabs(x),1.8))/7.0 + 1;
            verts[oidx+(points_per_side/2-2)-i][0] = -scale*x;
            verts[oidx+(points_per_side/2-2)-i][1] = scale*y;
            x -= 1.0/(points_per_side/2);
        }
        oidx += points_per_side/2 - 1;
        
        // build left column by transposing top row
        for (int i=0; i < points_per_side-1; i++) {
            verts[oidx][0] = -verts[i][1];
            verts[oidx][1] = verts[i][0];
            oidx++;
        }
        
        // build bottom row by flipping top row
        for (int i=0; i < points_per_side-1; i++) {
            verts[oidx][0] = -verts[i][0];
            verts[oidx][1] = -verts[i][1];
            oidx++;
        }
        
        // build right column by transposing top row
        for (int i=0; i < points_per_side-1; i++) {
            verts[oidx][0] = verts[i][1];
            verts[oidx][1] = -verts[i][0];
            oidx++;
        }
        
        photosite_geom = new Polygon_geom(verts);
        
        printf("rounded-square photosite area = %lg\n", (dynamic_cast<Polygon_geom*>(photosite_geom))->compute_area());
    }
    
    if (tc_photosite_poly.isSet()) {
        if (tc_photosite_geom.isSet()) {
            printf("Warning: --photosite-poly overrides the --photosite-geom option %s that was specified.\n", tc_photosite_geom.getValue().c_str());
        }
        FILE* poly = fopen(tc_photosite_poly.getValue().c_str(), "rt");
        if (!poly) {
            printf("Error: Could not open specified photosite polygon file %s. Aborting\n", tc_photosite_poly.getValue().c_str());
            exit(-1);
        }
        int nverts;
        int nread = fscanf(poly, "%d", &nverts);
        const char poly_format[] =
            "Could not parse specified photosite polygon file. Expected format is:\n "
            "\t<number_of_vertices>\n"
            "\t<x1> <y1>\n"
            "\t<x2> <y2>\n"
            "\t...\n"
            "\t<xn> <yn>\n"
            "where vertices are specified in clockwise order. Only one polygon allowed\n";
        if (nread != 1 || nverts <= 2) {
            printf("%s", poly_format);
            exit(-1);
        }
        vector<cv::Vec2d> verts(nverts);
        for (int i = 0; i < nverts; i++) {
            nread = fscanf(poly, "%lf %lf", &verts[i][0], &verts[i][1]);
            if (nread != 2) {
                printf("%s", poly_format);
                exit(-1);
            }
        }
        fclose(poly);
        printf("successfully read in user-specified photosite polygon with %d vertices\n", nverts);
        
        photosite_geom = new Polygon_geom(verts);
    }
    
    printf("photosite area = %lg\n", (dynamic_cast<Polygon_geom*>(photosite_geom))->compute_area());


    // decide which PSF rendering algorithm to use
    Render_polygon* rect=0;
    switch (psf_type) {
        case Render_polygon::AIRY:
        case Render_polygon::AIRY_PLUS_BOX:
        case Render_polygon::AIRY_PLUS_4DOT_OLPF:
        case Render_polygon::WAVEFRONT:
        case Render_polygon::WAVEFRONT_PLUS_BOX:
        case Render_polygon::RECT:
        case Render_polygon::RECT_PLUS_BOX:
            rect = build_psf(psf_type, 
                *target_geom,
                *photosite_geom,
                tc_aperture.getValue(),
                tc_pitch.getValue(),
                tc_lambda.getValue(),
                tc_olpf_split.getValue(),
                tc_samples.isSet() ? tc_samples.getValue() : 0,
                tc_w020.getValue(),
                tc_w040.getValue(),
                tc_aperture2.getValue()
            );
            break;
        case Render_polygon::GAUSSIAN:
            if (tc_target_name.isSet()) {
                fprintf(stderr, "Error! Gaussian integral PSF not compatible with --target-poly (at the moment).\n");
                fprintf(stderr, "       Please specify one of the Airy PSFs, or gaussian-sampled (see -p option).\n");
                return -1;
            }
            if (tc_photosite_geom.isSet() || tc_photosite_poly.isSet()) {
                printf("Warning: photosite geometry (--photosite-geom or --photosite-poly) NOT used when rendering Gaussian PSF. Use airy/airy-box/airy-4dot-olpf PSF instead.\n");
                printf("Continuing to render using pure Gaussian PSF.\n");
            }
            rect = new Render_rectangle_integral(*dynamic_cast<Polygon_geom*>(target_geom), sigma);
            break;
        case Render_polygon::GAUSSIAN_SAMPLED:
            if (tc_photosite_geom.isSet() || tc_photosite_poly.isSet()) {
                printf("Warning: photosite geometry (--photosite-geom or --photosite-poly) NOT used when rendering Gaussian PSF. Use airy/airy-box/airy-4dot-olpf PSF instead.\n");
                printf("Continuing to render using pure Gaussian PSF.\n");
            }
            rect = &default_target;
            if (tc_samples.isSet()) {
                rect->initialize_samples(tc_samples.getValue());
            }
            break;
    }
    rect->set_img_dimensions(height, width, 
        psf_type != Render_polygon::RECT_PLUS_BOX && !tc_psf_theta.isSet(), 
        tc_psf_theta.isSet() ? tc_psf_theta.getValue()*M_PI/180.0 + M_PI/2.0 : M_PI/2.0
    );
    if (tc_psf_ratio.isSet()) {
        rect->set_psf_ratio(tc_psf_ratio.getValue());
    }

    if ((tc_roi_start_row.isSet() || tc_roi_start_col.isSet()) &&
        !(tc_roi_width.isSet() && tc_roi_height.isSet())) {
        printf("Error: You must specify --roi-width and --roi-height if you use the --roi-row or --roi-col options.\n");
        exit(-1);
    }
    if ((tc_roi_width.isSet() && !tc_roi_height.isSet()) ||
        (!tc_roi_width.isSet() && tc_roi_height.isSet())) {
        printf("Error: You must specify both --roi-width and --roi-height.\n");
        exit(-1);
    }
    int x_origin = 0;
    int y_origin = 0;
    if (tc_roi_width.isSet() && tc_roi_height.isSet()) {
        width = tc_roi_width.getValue();
        height = tc_roi_height.getValue();
        x_origin = tc_roi_start_col.getValue();
        y_origin = tc_roi_start_row.getValue();
        if (x_origin < 0 || y_origin < 0) {
            printf("Error: ROI origin must be non-negative.\n");
            exit(-1);
        }
        printf("User-specified ROI (row,col): (%d, %d) to (%d, %d)\n.", y_origin, x_origin, y_origin + height, x_origin + width);
    }
    
    cv::Mat img;
    if (use_16bit) {
        img = cv::Mat(height, width, CV_16UC1);
    } else {
        img = cv::Mat(height, width, CV_8UC1);
    }

    
    Noise_source* ns = 0;

    if (use_sensor_model) {
        ns = new Sensor_model_noise(
            img.rows*img.cols,
            tc_read_noise.getValue(),
            tc_pattern_noise.getValue(),
            tc_adc_gain.getValue(),
            tc_adc_depth.getValue()
        );
    } else {
        ns = new Additive_gaussian_noise(img.rows*img.cols, tc_noise.getValue());
    }
    
    size_t nthreads = std::thread::hardware_concurrency();
    ThreadPool tp (nthreads);
    
    if (!tc_profile.getValue()) {
        int crc = 0;
        Render_rows rr(img, *rect, *ns, blacklevel, whitelevel, use_gamma, use_16bit, 
            tc_roi_width.isSet() ? 0 : border,
            crc, tc_ascale.getValue(), x_origin, y_origin
        );
        printf("progress: 0%% ");
        fflush(stdout);
        Stride_range::parallel_for(rr, tp, height);
        printf("\n");
        imwrite(tc_out_name.getValue(), img);
        inject_png_icc_profile(tc_out_name.getValue(), tc_16.getValue() || tc_linear.getValue());
    } else {
        string profile_fname("profile.txt");
        if (tc_out_name.isSet()) {
            profile_fname = tc_out_name.getValue();
        }
        // render call
        const int oversample = 32*4;
        vector< pair<double, double> > esf(Render_esf::n_samples(rwidth, oversample));
        Render_esf re(*rect, esf, rwidth, theta, oversample, tc_xoff.getValue(), tc_yoff.getValue());
        Stride_range::parallel_for(re, tp, esf.size());
        re.write(profile_fname);
    }
    
    if (ns && ns->clip_count > 0) {
        printf("\nWarning: %lg%% of pixel intensities were clipped!\n  Check your black/white levels relative to noise magnitude.\n\n",
            100.0*double(ns->clip_count) / double(ns->size)
        );
    }
    
    if (display_mtf_equation) { // for now we only output equations for square photosites
        printf("MTF curve:  %s\n", rect->get_mtf_curve().c_str());
        printf("PSF : %s\n", rect->get_psf_curve().c_str());
        printf("Note: above MTF/PSF functions apply strictly to an ESF taken perpendicular to a step edge.\n");
        printf("Note: edge orientation set using '-a' is NOT taken into account here\n");
        printf("MTF50 = %lf\n", rect->get_mtf50_value());    
    }

    delete ns;
    if (rect != &default_target) {
        delete rect;
    }
    
    return 0;
}
