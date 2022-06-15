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
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <string.h>
#include <locale.h>

#include <tclap/CmdLine.h>

#include <opencv2/imgcodecs/imgcodecs.hpp>

using std::string;
using std::stringstream;

#include "include/logger.h"
Logger logger;

#include "include/common_types.h"
#include "include/thresholding.h"
#include "include/mtf_core.h"
#include "include/mtf_core_tbb_adaptor.h"
#include "include/output_version.h"
#include "include/ca_core.h"
#include "include/ca_core_tbb_adaptor.h"
#include "include/mtf_renderer_annotate.h"
#include "include/mtf_renderer_profile.h"
#include "include/mtf_renderer_mfprofile.h"
#include "include/mtf_renderer_grid.h"
#include "include/mtf_renderer_print.h"
#include "include/mtf_renderer_stats.h"
#include "include/mtf_renderer_sfr.h"
#include "include/mtf_renderer_esf.h"
#include "include/mtf_renderer_edges.h"
#include "include/mtf_renderer_lensprofile.h"
#include "include/mtf_renderer_chart_orientation.h"
#include "include/mtf_renderer_focus.h"
#include "include/ca_renderer_print.h"
#include "include/ca_renderer_grid.h"
#include "include/scanline.h"
#include "include/distance_scale.h"
#include "include/auto_crop.h"
#include "include/imatest_crop.h"
#include "include/bayer.h"
#include "include/demosaic.h"
#include "include/stride_range.h"
#include "include/distortion_optimizer.h"
#include "include/undistort_rectilinear.h"
#include "include/undistort_equiangular.h"
#include "include/undistort_stereographic.h"
#include "include/esf_sampler.h"
#include "include/tiffsniff.h"
#include "include/display_profile.h"
#include "include/esf_model_kernel.h"
#include "include/esf_model_loess.h"
#include "include/job_metadata.h"
#include "config.h"

//-----------------------------------------------------------------------------
void print_version_info(void) {
    printf("MTF mapper version %d.%d.%d\n", mtfmapper_VERSION_MAJOR, mtfmapper_VERSION_MINOR, mtfmapper_VERSION_SUB);
}

//-----------------------------------------------------------------------------
int main(int argc, char** argv) {
    setlocale(LC_ALL, "C");

    stringstream ss;
    ss << mtfmapper_VERSION_MAJOR << "." << mtfmapper_VERSION_MINOR << "." << mtfmapper_VERSION_SUB;
    
    TCLAP::CmdLine cmd("Measure MTF50 values across edges of rectangular targets", ' ', ss.str());
    TCLAP::UnlabeledValueArg<std::string>  tc_in_name("<input_filename>", 
        "Input image file name (many extensions supported)", true, "input.png", "image_filename", cmd
    );
    TCLAP::UnlabeledValueArg<std::string>  tc_wdir("<working_directory>", 
        "Working directory (for output files); \".\" is fine", true, ".", "directory", cmd
    );
    TCLAP::SwitchArg tc_profile("p","profile","Generate MTF50 profile", cmd, false);
    TCLAP::SwitchArg tc_mf_profile("","mf-profile","Generate MTF50 profile (Manual focus chart)", cmd, false);
    TCLAP::SwitchArg tc_annotate("a","annotate","Annotate input image with MTF50 values", cmd, false);
    TCLAP::SwitchArg tc_surface("s","surface","Generate MTF50 surface plots", cmd, false);
    TCLAP::SwitchArg tc_linear("l","linear","Input image is linear 8-bit (default for 8-bit is assumed to be sRGB gamma corrected)", cmd, false);
    TCLAP::SwitchArg tc_print("r","raw","Print raw MTF50 values", cmd, false);
    TCLAP::SwitchArg tc_edges("q","edges","Print raw MTF50 values, grouped by edge location", cmd, false);
    TCLAP::SwitchArg tc_sfr("f","sfr","Store raw SFR curves for each edge", cmd, false);
    TCLAP::SwitchArg tc_esf("e","esf","Store raw ESF and PSF curves for each edge", cmd, false);
    TCLAP::SwitchArg tc_lensprof("","lensprofile","Render M/S lens profile plot", cmd, false);
    TCLAP::SwitchArg tc_chart_orientation("","chart-orientation","Visualize chart orientation relative to camera", cmd, false);
    TCLAP::SwitchArg tc_border("b","border","Add a border of 20 pixels to the image", cmd, false);
    TCLAP::SwitchArg tc_absolute("","absolute-sfr","Generate absolute SFR curve (MTF) i.s.o. relative SFR curve", cmd, false);
    TCLAP::SwitchArg tc_smooth("","nosmoothing","Disable SFR curve (MTF) smoothing", cmd, false);
    TCLAP::SwitchArg tc_autocrop("","autocrop","Automatically crop image to the chart area", cmd, false);
    TCLAP::SwitchArg tc_focus("","focus","Compute focus depth using special 'focus' type chart", cmd, false);
    TCLAP::SwitchArg tc_log_append("", "log-append", "Append to log file in stead of overwriting log file", cmd, false);
    TCLAP::SwitchArg tc_debug("", "debug", "Enable debug output messages", cmd, false);
    TCLAP::SwitchArg tc_single_roi("", "single-roi", "Treat the entire input image as the ROI", cmd, false);
    TCLAP::SwitchArg tc_distort_opt("", "optimize-distortion", "Optimize lens distortion coefficients", cmd, false);
    TCLAP::SwitchArg tc_rectilinear("", "rectilinear-equivalent", "Measure MTF in rectilinear equivalent projection", cmd, false);
    TCLAP::SwitchArg tc_distort_crop("", "no-undistort-crop", "Do not crop undistorted image (equiangular, stereographic)", cmd, false);
    TCLAP::SwitchArg tc_full_sfr("", "full-sfr", "Output the full SFR/MTF curve (up to 4 c/p) when combined with -q or -f", cmd, false);
    TCLAP::SwitchArg tc_ima_mode("", "imatest-chart", "Treat the input image as an Imatest chart (crop black bars out)", cmd, false);
    TCLAP::SwitchArg tc_lensprofile_fixed("", "lensprofile-fixed-size", "Lens profile output is scaled to maximum sensor radius", cmd, false);
    TCLAP::SwitchArg tc_monotonic_filter("", "monotonic-esf-filter", "Force the application of a monontonic ESF noise filter (use at own risk!)", cmd, false);
    TCLAP::SwitchArg tc_ca("", "ca", "Estimate chromatic aberration", cmd, false);
    TCLAP::SwitchArg tc_ca_fraction("", "ca-fraction", "Chromatic aberration image generated in units of radial distance fraction", cmd, false);
    TCLAP::SwitchArg tc_jpeg("", "jpeg", "Annotated image saved in JPEG format to gain speed", cmd, false);
    TCLAP::SwitchArg tc_checkerboard("", "checkerboard", "Process the input image as a checkerboard pattern", cmd, false);
    TCLAP::SwitchArg tc_ca_all("", "ca-all-edges", "Chromatic aberration is calculated on all edges, not just tangential edges", cmd, false);
    TCLAP::SwitchArg tc_allow_partial("", "allow-partial", "Allow partial targets (cut by image boundary) to be processed", cmd, false);
    TCLAP::SwitchArg tc_invert("", "invert", "Invert the image brightness before processing (white targets on black backgrounds)", cmd, false);
    #ifdef MDEBUG
    TCLAP::SwitchArg tc_bradley("", "bradley", "Use Bradley thresholding i.s.o Sauvola thresholding", cmd, false);
    #endif
    TCLAP::ValueArg<double> tc_angle("g", "angle", "Angular filter [0,360)", false, 0, "angle", cmd);
    TCLAP::ValueArg<double> tc_snap("", "snap-angle", "Snap-to angle modulus [0,90)", false, 1000, "angle", cmd);
    TCLAP::ValueArg<double> tc_thresh("t", "threshold", "Dark object threshold (0,1), default 0.55", false, 0.55, "threshold", cmd);
    TCLAP::ValueArg<string> tc_gnuplot("", "gnuplot-executable", "Full path (including filename) to gnuplot executable ", false, "gnuplot", "filepath", cmd);
    TCLAP::ValueArg<double> tc_pixelsize("", "pixelsize", "Pixel size in microns. This also switches units to lp/mm", false, 1.0, "size", cmd);
    TCLAP::ValueArg<double> tc_lp1("", "lp1", "Lens profile resolution 1 (lp/mm or c/p)", false, 10.0, "lp/mm", cmd);
    TCLAP::ValueArg<double> tc_lp2("", "lp2", "Lens profile resolution 2 (lp/mm or c/p)", false, 30.0, "lp/mm", cmd);
    TCLAP::ValueArg<double> tc_lp3("", "lp3", "Lens profile resolution 3 (lp/mm or c/p)", false, 50.0, "lp/mm", cmd);
    TCLAP::ValueArg<string> tc_logfile("", "logfile", "Output written to <logfile> in stead of standard out", false, "", "filename", cmd);
    TCLAP::ValueArg<int> tc_gpwidth("", "gnuplot-width", "Width of images rendered by gnuplot", false, 1024, "pixels", cmd);
    TCLAP::ValueArg<double> tc_focal("", "focal-ratio", "Specify focal ratio for use in chart orientation estimation", false, -2, "ratio", cmd);
    TCLAP::ValueArg<double> tc_equiangular("", "equiangular", "Treat input image as equi-angular mapping (fisheye) with the specified focal length", false, 16.0, "focal length(mm)", cmd);
    TCLAP::ValueArg<double> tc_stereographic("", "stereographic", "Treat input image as stereographic mapping (fisheye) with the specified focal length", false, 8.0, "focal length(mm)", cmd);
    TCLAP::ValueArg<double> tc_zscale("", "zscale", "Z-axis scaling of '-s' outputs [0,1]. A value of 0 means z-axis scale starts at zero, and 1.0 means z-axis starts from minimum measurement", false, 0.0, "scale factor", cmd);
    TCLAP::ValueArg<double> tc_thresh_win("", "threshold-window", "Fraction of min(img width, img height) to use as window size during thresholding; range (0,1]", false, 0.33333, "fraction", cmd);
    TCLAP::ValueArg<double> tc_mtf_contrast("", "mtf", "Specify target contrast, e.g., --mtf 30 yields MTF30 results. Range [10, 90], default is 50", false, 50.0, "percentage", cmd);
    TCLAP::ValueArg<double> tc_alpha("", "alpha", "Standard deviation of smoothing kernel [1,20]", false, 13, "unitless", cmd);
    TCLAP::ValueArg<double> tc_surface_max("", "surface-max", "Specify maximum value in MTF50 surface plots", false, -1, "units depend on other settings", cmd);
    TCLAP::ValueArg<string> tc_roi_file("", "roi-file", "Only process ROIs defined in <roifile>, rather than using automatic target selection", false, "", "<roifile>", cmd);
    TCLAP::ValueArg<int> tc_checkerboard_radius("", "checkerboard-radius", "Radius of dilation structuring element when processing checkerboard images", false, 2, "pixels", cmd);
    #ifdef MDEBUG
    TCLAP::ValueArg<double> tc_ridge("", "ridge", "Specify ridge regression parameter [0,+infy)", false, 5e-8, "unitless", cmd);
    TCLAP::ValueArg<double> tc_noise_seed("", "noise-seed", "Image noise seed", false, 10, "unitless", cmd);
    TCLAP::ValueArg<double> tc_noise_sd("", "noise-sd", "Image noise sd", false, 0, "unitless", cmd);
    TCLAP::SwitchArg tc_single("","single-threaded","Force single-threaded operation", cmd, false);
    #endif
    
    TCLAP::ValuesConstraint<int> output_version_constraints(Output_version::instance().get_versions());
    TCLAP::ValueArg<int> tc_output_version("v", "output-version", "Specify output format version (affects -q)", false, 1, &output_version_constraints);
    cmd.add(tc_output_version);

    vector<string> allowed_bayer_subsets;
    allowed_bayer_subsets.push_back("red");
    allowed_bayer_subsets.push_back("green");
    allowed_bayer_subsets.push_back("blue");
    allowed_bayer_subsets.push_back("none");
    TCLAP::ValuesConstraint<string> bayer_constraints(allowed_bayer_subsets);
    TCLAP::ValueArg<std::string> tc_bayer("", "bayer", "Select Bayer subset", false, "none", &bayer_constraints );
    cmd.add(tc_bayer);
    
    vector<string> allowed_esf_samplers;
    for (const auto& s: Esf_sampler::esf_sampler_names) {
        allowed_esf_samplers.push_back(s);
    }
    TCLAP::ValuesConstraint<string> esf_sampler_constraints(allowed_esf_samplers);
    TCLAP::ValueArg<std::string> tc_esf_sampler("", "esf-sampler", "Select ESF sampler type", false, "piecewise-quadratic", &esf_sampler_constraints);
    cmd.add(tc_esf_sampler);
    
    vector<string> allowed_cfa_patterns;
    allowed_cfa_patterns.push_back("rggb");
    allowed_cfa_patterns.push_back("bggr");
    allowed_cfa_patterns.push_back("grbg");
    allowed_cfa_patterns.push_back("gbrg");
    TCLAP::ValuesConstraint<string> cfa_pattern_constraints(allowed_cfa_patterns);
    TCLAP::ValueArg<std::string> tc_cfa_pattern("", "cfa-pattern", "Select CFA pattern", false, "rggb", &cfa_pattern_constraints);
    cmd.add(tc_cfa_pattern);
    
    vector<string> allowed_esf_models;
    for (const auto& s: Esf_model::esf_model_names) {
        allowed_esf_models.push_back(s);
    }
    TCLAP::ValuesConstraint<string> esf_model_constraints(allowed_esf_models);
    TCLAP::ValueArg<std::string> tc_esf_model("", "esf-model", "Select ESF model type", false, "kernel", &esf_model_constraints);
    cmd.add(tc_esf_model);
    
    cmd.parse(argc, argv);
    
    string esf_sampler_name = tc_esf_sampler.getValue();

    if (tc_logfile.isSet()) {
        logger.redirect(tc_logfile.getValue(), tc_log_append.getValue());
    }
    if (tc_debug.getValue()) {
        logger.enable_level(Logger::LOGGER_DEBUG);
    }
    
    bool lpmm_mode = false;
    double pixel_size = 1;
    if (tc_pixelsize.isSet()) {
        logger.info("%s\n", "Info: Pixel size has been specified, measurement will be reported in lp/mm, rather than c/p");
        lpmm_mode = true;
        pixel_size = 1000 / tc_pixelsize.getValue();
        logger.info("working with=%lf pixels per mm\n", pixel_size);
    }
    
    Job_metadata job_metadata(
        pixel_size,
        0.5, // will be updated later after clamping
        Bayer::from_string(tc_bayer.getValue())
    );

    if (!tc_profile.isSet() && !tc_annotate.isSet() && !tc_surface.isSet() && !tc_print.isSet() && !tc_sfr.isSet() && !tc_edges.isSet()) {
        logger.info("%s\n", "Warning: No output specified. You probably want to specify at least one of the following flags: [-r -p -a -s -f -q]");
    }

    cv::Mat cvimg;
    try {
        cvimg = cv::imread(tc_in_name.getValue(),-1);
        job_metadata.channels = cvimg.channels();
        job_metadata.set_image_dims(cvimg.cols, cvimg.rows);
    } catch (const cv::Exception& ex) {
        cout << ex.what() << endl;
    }

    if (!cvimg.data) {
        logger.error("Fatal error: could not open input file <%s>.\nFile is missing, or not where you said it would be, or you do not have read permission.\n", tc_in_name.getValue().c_str());
        return 2;
    }

    if (!(cvimg.depth() == CV_8U || cvimg.depth() == CV_16U)) {
        logger.error("%s\n", "Fatal error: Invalid image type. Only 8-bit unsigned and 16-bit unsigned integer images supported.");
        return 5;
    }
    
    struct STAT sb;
    if (STAT(tc_wdir.getValue().c_str(), &sb) != 0) {
        logger.error("Fatal error: specified output directory <%s> does not exist\n", tc_wdir.getValue().c_str());
        return 3;
    } else {
        if (!S_ISDIR(sb.st_mode)) {
            logger.error("Fatal error: speficied output directory <%s> is not a directory\n", tc_wdir.getValue().c_str());
            return 3;
        }
    }
    
    Tiffsniff tiff(tc_in_name.getValue(), cvimg.elemSize1() == 1);
    Display_profile display_profile;
    if (tiff.profile_found()) {
        display_profile = tiff.profile();
    } else {
        if (cvimg.elemSize1() == 1 && !tc_linear.getValue()) {
            display_profile.force_sRGB();
        }
    }
	
    if (tc_linear.getValue()) {
        display_profile.force_linear();
    }
    
    if (cvimg.channels() == 4) {
        logger.info("%s\n", "Input image had 4 channels. Only the first 3 will be used.");
        cv::Mat reduced(cvimg.rows, cvimg.cols, (cvimg.elemSize1() == 1) ? CV_8UC3 : CV_16UC3);
        int from_to[] = {0,0, 1,1, 2,2};
        cv::mixChannels(&cvimg, 1, &reduced, 1, from_to, 3);
        cvimg = reduced;
    }
    
    int in_num_channels = cvimg.channels();
    cv::Mat rgb_img;
    
    if (tc_ca.getValue() && cvimg.channels() == 3) {
        rgb_img = cvimg; // TODO: proper linearization ?
    }
    
    cvimg = display_profile.to_luminance(cvimg);
    
    assert(cvimg.type() == CV_16UC1);

    const int border_width = 100;
    if (tc_border.getValue()) {
        logger.info("The -b option has been specified, adding a %d-pixel border to the image\n", border_width);
        double max_val = 0;
        double min_val = 0;
        cv::minMaxLoc(cvimg, &min_val, &max_val);
        cv::Mat border;
        cv::copyMakeBorder(cvimg, border, border_width, border_width, border_width, border_width, cv::BORDER_CONSTANT, cv::Scalar((int)max_val));
        cvimg = border;
    }
    
    if (tc_equiangular.isSet() && !tc_pixelsize.isSet()) {
        logger.error("%s\n", "Fatal error: You must specify the pixel size (pitch) with the --pixelsize option when using --equiangular option. Aborting.");
        return 1;
    }
    
    if (tc_stereographic.isSet() && !tc_pixelsize.isSet()) {
        logger.error("%s\n", "Fatal error: You must specify the pixel size (pitch) with the --pixelsize option when using --stereographic option. Aborting.");
        return 1;
    }
    
    if (tc_stereographic.isSet() && tc_equiangular.isSet()) {
        logger.error("%s\n", "Fatal error: You may only specify one of --stereographic and --equiangular. Aborting.");
        return 1;
    }
    
    int gnuplot_width = std::max(1024, tc_gpwidth.getValue());
    
    // process working directory
    std::string wdir(tc_wdir.getValue());
    if (wdir[wdir.length()-1]) {
        wdir = tc_wdir.getValue() + "/";
    }

    char slashchar='/';
    #ifdef _WIN32
    // on windows, mangle the '/' into a '\\'
    std::string wdm;
    for (size_t i=0; i < wdir.length(); i++) {
        if (wdir[i] == '/') {
            wdm.push_back('\\');
            wdm.push_back('\\');
        } else {
            wdm.push_back(wdir[i]);
        }
    }
    wdir = wdm;
    slashchar='\\';
    #endif
    
    // strip off supposed extention suffix,
    // and supposed path prefix
    std::string in_img_name = tc_in_name.getValue();
    std::replace(in_img_name.begin(), in_img_name.end(), '/', slashchar);
    int ext_idx=-1;
    int path_idx=0;
    for (int idx= in_img_name.length()-1; idx >= 0 && path_idx == 0; idx--) {
        if (in_img_name[idx] == '.' && ext_idx < 0) {
            ext_idx = idx;
        }
        if (in_img_name[idx] == slashchar && path_idx == 0) {
            path_idx = idx;
        }
    }
    if (ext_idx < 0) {
        ext_idx = in_img_name.length();
    }
    std::string img_filename;
    for (int idx=path_idx; idx < ext_idx; idx++) {
        char c = in_img_name[idx];
        if (c == slashchar) continue;
        if (c == '_') {
            img_filename.push_back('\\');
            img_filename.push_back('\\');
        }
        img_filename.push_back(c);
    }
	
    cv::Mat masked_img;
    
    cv::Rect img_dimension_correction(0,0, cvimg.cols, cvimg.rows);
    
    if (tc_autocrop.getValue()) {
        Auto_cropper ac(cvimg);
        cvimg = ac.subset(cvimg, &img_dimension_correction);
    }
    if (tc_ima_mode.getValue()) {
        Imatest_cropper ic(cvimg);
        ic.fill_bars(cvimg);
    }
    
    cv::Mat rawimg = cvimg;
    if (tc_bayer.isSet()) {
        simple_demosaic(cvimg, rawimg, 
            Bayer::from_cfa_string(tc_cfa_pattern.getValue()), 
            Bayer::from_string(tc_bayer.getValue()), tc_single_roi.getValue()
        );
        //imwrite(string("prewhite.png"), rawimg);
        //imwrite(string("white.png"), cvimg);
    }
    
    if (tc_invert.getValue()) {
        invert(cvimg);
        if (rawimg.data != cvimg.data) {
            invert(rawimg);
        }
    }
    
    Undistort* undistort = nullptr;
    if (tc_equiangular.isSet()) {
        logger.info("Treating input image as equi-angular with focal length %.2lf, unmapping\n", tc_equiangular.getValue());
        undistort = new Undistort_equiangular(img_dimension_correction, tc_equiangular.getValue(), tc_pixelsize.getValue()/1000.0);
        undistort->set_rectilinear_equivalent(tc_rectilinear.getValue());
    }
    if (tc_stereographic.isSet()) {
        logger.info("Treating input image as stereographic with focal length %.2lf, unmapping\n", tc_stereographic.getValue());
        undistort = new Undistort_stereographic(img_dimension_correction, tc_stereographic.getValue(), tc_pixelsize.getValue()/1000.0);
        undistort->set_rectilinear_equivalent(tc_rectilinear.getValue());
    }
    if (undistort) {
        undistort->set_allow_crop(!tc_distort_crop.getValue());
        if (esf_sampler_name.compare("deferred") != 0) {
            logger.info("Warning: Because you specified an undistortion model,"
                " your ESF sampler choice of '%s' has been changed to 'deferred'\n", 
                esf_sampler_name.c_str()
            );
        }
    }
    
    if (esf_sampler_name.compare("deferred") == 0 && !undistort) {
        if (!tc_distort_opt.getValue()) {
            logger.error("%s\n", "Error: Deferred ESF sampler cannot be used if no undistortion model is specified."
                "See '--optimize-distortion, --equiangular, or --stereographic' options");
            return 1;
        }
    }
    
    if (tc_focus.isSet() || tc_mf_profile.getValue()) {
        esf_sampler_name = "line";
        logger.info("%s\n", "Note: because --focus output option was selected, the --esf_sampler option has been changed to \"line\".");
    }
    
    bool finished;
    bool distortion_applied = false;
    do {
        finished = true;
    
        if (undistort) {
            cvimg = undistort->unmap(cvimg, rawimg);
        }
        
        logger.info("%s\n", "Thresholding image ...");
        int brad_S = tc_border.getValue() ? 
            max(cvimg.cols, cvimg.rows) : 
            max(tc_thresh_win.isSet() ? 20 : 500, int(min(cvimg.cols, cvimg.rows)*tc_thresh_win.getValue()));
        double brad_threshold = tc_thresh.getValue();
        #ifdef MDEBUG
            if (tc_bradley.getValue()) {
                printf("using Bradley thresholding\n");
                bradley_adaptive_threshold(cvimg, masked_img, brad_threshold, brad_S);
            } else {
                sauvola_adaptive_threshold(cvimg, masked_img, brad_threshold/0.55*0.85, brad_S);
            }
        #else
            sauvola_adaptive_threshold(cvimg, masked_img, brad_threshold/0.55*0.85, brad_S); // fudge the threshold to maintain backwards compatibility
        #endif
        
        const int erosion_size = std::min(5, std::max(1, tc_checkerboard_radius.getValue()));
        if (tc_checkerboard.getValue()) {
            // first do a small-scale closing operation to prevent
            // small interior / boundary holes from growing in the subsequent dilation
            cv::Mat open_element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(3, 3));
            cv::morphologyEx(masked_img, masked_img, cv::MORPH_OPEN, open_element);
            
            // dilate masked image to break checkerboard corners
            cv::Mat element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(2*erosion_size + 1, 2*erosion_size + 1));
            cv::morphologyEx(masked_img, masked_img, cv::MORPH_DILATE, element);
        }
        
        logger.info("%s\n", "Computing gradients ...");
        Gradient gradient(cvimg);
        
        logger.info("%s\n", "Component labelling ...");
        Component_labeller::zap_borders(masked_img);
        // largest component boundary length determined empirically
        const int64_t boundary_long_side = 2*std::max(cvimg.rows, cvimg.cols)*0.4;
        const int64_t boundary_short_side = 2*std::min(cvimg.rows, cvimg.cols)*0.4;
        const int64_t max_boundary_length = std::max(int64_t(8000), boundary_long_side + boundary_short_side);
        Component_labeller cl(masked_img, 60, false, max_boundary_length);

        if (cl.get_boundaries().size() == 0 && !(tc_single_roi.getValue() || tc_roi_file.isSet())) {
            logger.error("%s\n", "Error: No black objects found. Try a lower threshold value with the -t option.");
            return 4;
        }
        
        // fudge the image gradient on the edges
        if (tc_allow_partial.getValue()) {
            logger.info("%s\n", "Fudging image gradient where target objects intersect the image edges.\n");
            for (int col=0; col < cvimg.cols; col++) {
                if (cl(col, 2) > 0) {
                    gradient.write_grad_y(col, 2, -3.0);
                    gradient.write_grad_x(col, 2, 0.0);
                }
                if (cl(col, cvimg.rows - 3) > 0) {
                    gradient.write_grad_y(col, cvimg.rows - 3, 3.0);
                    gradient.write_grad_x(col, cvimg.rows - 3, 0.0);
                }
            }
            for (int row=0; row < cvimg.rows; row++) {
                if (cl(2, row) > 0) {
                    gradient.write_grad_y(2, row, 0.0);
                    gradient.write_grad_x(2, row, -3.0);
                }
                if (cl(cvimg.cols - 3, row) > 0) {
                    gradient.write_grad_y(cvimg.cols - 3, row, 0.0);
                    gradient.write_grad_x(cvimg.cols - 3, row, 3.0);
                }
            }
        }
        
        // try to restore the object boundaries to compensate for dilation of thresholded image
        if (tc_checkerboard.getValue()) {
            cl.inflate_boundaries(erosion_size);
        }
        
        // now we can destroy the thresholded image
        masked_img = cv::Mat(1,1, CV_8UC1);
        
        Mtf_core mtf_core(
            cl, gradient, cvimg, rawimg, tc_bayer.getValue(), tc_cfa_pattern.getValue(),
            undistort ? "deferred" : (tc_distort_opt.getValue() ? "line" : esf_sampler_name), // force deferred sampler if an undistortion model is specified
            undistort, 
            tc_border.getValue() ? border_width+1 : 0 
        );
        mtf_core.set_absolute_sfr(tc_absolute.getValue());
        mtf_core.set_sfr_smoothing(!tc_smooth.getValue());
        mtf_core.set_allow_partial(tc_allow_partial.getValue());
        if (tc_border.getValue()) {
            logger.debug("setting border to %d\n", border_width);
        }
        if (tc_esf_model.getValue().compare("kernel") == 0) {
            mtf_core.set_esf_model(std::unique_ptr<Esf_model>(new Esf_model_kernel()));
        } else { // only alternative at the moment is "loess"
            #ifdef MDEBUG
            mtf_core.set_esf_model(
                std::unique_ptr<Esf_model>(
                    new Esf_model_loess(
                        tc_alpha.isSet() ? tc_alpha.getValue() : 5.5,
                        tc_ridge.getValue()
                    )
                )
            );
            #else
            mtf_core.set_esf_model(std::unique_ptr<Esf_model>(new Esf_model_loess()));
            #endif
        }
        mtf_core.get_esf_model()->set_monotonic_filter(tc_monotonic_filter.getValue());
        
        if (tc_alpha.isSet()) {
            mtf_core.set_esf_model_alpha_parm(tc_alpha.getValue());
        }
        
        if (tc_snap.isSet()) {
            mtf_core.set_snap_angle(tc_snap.getValue()/180*M_PI);
        }
        if (tc_focus.isSet() || tc_mf_profile.getValue()) {
            mtf_core.set_sliding(true);
            if (tc_mf_profile.getValue()) {
                mtf_core.set_samples_per_edge(5);
            }
        }
        if (tc_chart_orientation.isSet()) {
            mtf_core.set_find_fiducials(true);
        }
        
        if (tc_distort_opt.getValue() && !distortion_applied) {
            mtf_core.set_ridges_only(true);
        }

        if (tc_full_sfr.getValue()) {
            mtf_core.use_full_sfr();
        }
        
        if (tc_mtf_contrast.isSet()) {
            double contrast = tc_mtf_contrast.getValue();
            if (contrast < 10) {
                if (contrast < 1) {
                    logger.error("Warning: Requested MTF%02d, clamped to MTF01 instead\n", int(contrast));
                    contrast = 1;
                } else {
                    logger.error("Warning: Requested MTF%02d, which is highly likely to be affected by noise, and may cause some edges not be detected.\n", int(contrast));
                }
            }
            if (contrast > 90) {
                logger.error("Warning: Requested MTF%02d, clamped to MTF90 instead\n", int(contrast));
                contrast = 90;
            }
            mtf_core.set_mtf_contrast(contrast / 100.0);
            job_metadata.mtf_contrast = contrast / 100.0;
        }
        #ifdef MDEBUG
        mtf_core.noise_seed = tc_noise_seed.getValue();
        mtf_core.noise_sd = tc_noise_sd.getValue();
        #endif
        
        Mtf_core_tbb_adaptor ca(&mtf_core);
        
        if (tc_single_roi.getValue()) {
            mtf_core.process_image_as_roi(cv::Rect2i(0, 0, cvimg.cols, cvimg.rows));
        } else {
            if (tc_roi_file.isSet()) {
                mtf_core.process_manual_rois(tc_roi_file.getValue());
            } else {
                #ifdef MDEBUG
                if (tc_single.getValue()) {
                    ca(Stride_range(size_t(0), mtf_core.num_objects()-1, 1));
                } else {
                    logger.debug("Parallel MTF%2d calculation\n", int(mtf_core.get_mtf_contrast()*100));
                    Stride_range::parallel_for(ca, ThreadPool::instance(), mtf_core.num_objects());
                }
                #else
                Stride_range::parallel_for(ca, ThreadPool::instance(), mtf_core.num_objects());
                #endif
            }
        }
        
        if (mtf_core.get_blocks().size() == 0 && !(tc_focus.getValue() || tc_mf_profile.getValue())) {
            logger.error("%s\n", "Error: No suitable target objects found.");
            return 4;
        }
        
        if (tc_distort_opt.getValue() && !distortion_applied) { 
            Distortion_optimizer dist_opt(mtf_core.get_blocks(), Point2d(rawimg.cols/2, rawimg.rows/2));
            dist_opt.solve();
            logger.info("Optimal distortion coefficients: %lg %lg\n", dist_opt.best_sol[0], dist_opt.best_sol[1]);
            
            vector<double> coeffs(2);
            for (int i=0; i < 2; i++) {
                coeffs[i] = dist_opt.best_sol[i];
            }
            undistort = new Undistort_rectilinear(img_dimension_correction, coeffs);
            undistort->set_max_val(dist_opt.get_max_val());
            
            finished = false;
            distortion_applied = true;
            logger.info("%s\n", "Performing second pass on undistorted image.");
            continue; // effectively jump back to the start
        }
        
        if (tc_ca.getValue()) {
            Ca_core chromatic(mtf_core);
            if (tc_ca_all.getValue()) {
                chromatic.set_allow_all_edges();
            }
        
            if (in_num_channels == 3) {
                logger.info("%s\n", "Using original RGB input image to estimate CA.");
                vector<cv::Mat> channels = display_profile.to_linear_rgb(rgb_img);

                if (undistort) {
                    undistort->apply_padding(channels);
                }

                chromatic.set_rgb_channels(channels);
            }
            
            Ca_core_tbb_adaptor ca_adaptor(chromatic);
            
            size_t num_blocks = mtf_core.get_blocks().size();
            if (num_blocks >= 1) {
                #ifdef MDEBUG
                if (tc_single.getValue()) {
                    ca_adaptor(Stride_range(size_t(0), num_blocks - 1, 1));
                } else {
                    logger.info("Parallel CA calculation with %ld blocks\n", num_blocks);
                    Stride_range::parallel_for(ca_adaptor, ThreadPool::instance(), num_blocks);
                }
                #else
                logger.info("Parallel CA calculation with %ld blocks\n", num_blocks);
                Stride_range::parallel_for(ca_adaptor, ThreadPool::instance(), num_blocks);
                #endif
            }
        }
        
        Distance_scale distance_scale;
        if (tc_mf_profile.getValue() || tc_focus.getValue() || tc_chart_orientation.getValue()) {
            distance_scale.construct(mtf_core, true, &img_dimension_correction, tc_focal.getValue(), wdir + string("fiducial_correspondence.txt"));
        }
        
        // release most of the resources we no longer need
        masked_img.release();
        gradient.release();
        cl.release();
        
        // now render the computed MTF values
        if (tc_annotate.getValue()){
            Mtf_renderer_annotate annotate(cvimg, wdir + string("annotated"), lpmm_mode, pixel_size, tc_jpeg.getValue());
            annotate.render(mtf_core.get_blocks());
        }
        
        bool gnuplot_warning = true;
        bool few_edges_warned = false;
        
        if (tc_profile.getValue()) {
            if (mtf_core.get_blocks().size() < 10) {
                logger.info("Warning: fewer than 10 edges found, so MTF%2d surfaces/profiles will not be generated. Are you using suitable input images?\n", int(mtf_core.get_mtf_contrast()*100));
                few_edges_warned = true;
            } else {
                Mtf_renderer_profile profile(
                    img_filename,
                    wdir, 
                    string("profile.txt"),
                    tc_gnuplot.getValue(),
                    cvimg,
                    gnuplot_width,
                    lpmm_mode,
                    pixel_size,
                    int(mtf_core.get_mtf_contrast()*100)
                );
                profile.render(mtf_core.get_blocks());
                gnuplot_warning = !profile.gnuplot_failed();
            }
        }
        
        if (tc_mf_profile.getValue()) {
            Mtf_renderer_mfprofile profile(
                distance_scale,
                wdir, 
                string("focus_peak.png"),
                cvimg,
                lpmm_mode,
                pixel_size
            );
            profile.render(mtf_core.get_samples());
        }

        if (tc_surface.getValue()) {
            if (mtf_core.get_blocks().size() < 10) {
                if (!few_edges_warned) {
                    logger.info("Warning: fewer than 10 edges found, so MTF%2d surfaces/profiles will not be generated. Are you using suitable input images?\n", int(mtf_core.get_mtf_contrast()*100));
                }
            } else {
                Mtf_renderer_grid grid(
                    img_filename,
                    wdir, 
                    string("grid.txt"),
                    tc_gnuplot.getValue(),
                    cvimg,
                    gnuplot_width,
                    lpmm_mode,
                    pixel_size,
                    tc_zscale.getValue(),
                    tc_surface_max.getValue(),
                    lrint(mtf_core.get_mtf_contrast()*100)
                );
                grid.set_gnuplot_warning(gnuplot_warning);
                grid.set_sparse_chart(tc_ima_mode.getValue());
                grid.render(mtf_core.get_blocks());
            }
        }
        
        if (tc_print.getValue()) {
            Mtf_renderer_print printer(
                wdir + string("raw_mtf_values.txt"), 
                tc_angle.isSet(), 
                tc_angle.getValue()/180.0*M_PI,
                lpmm_mode,
                pixel_size
            );
            printer.render(mtf_core.get_blocks());
        }
        
        if (tc_edges.getValue()) {
            Mtf_renderer_edges printer(
                wdir + string("edge_mtf_values.txt"), 
                wdir + string("edge_sfr_values.txt"),
                wdir + string("edge_line_deviation.txt"),
                wdir + string("serialized_edges.bin"),
                Output_version::type(tc_output_version.getValue()),
                job_metadata,
                lpmm_mode, pixel_size
            );
            printer.render(mtf_core.get_blocks());
        }
        
        if (tc_lensprof.getValue()) {
        
            vector<double> resolutions;
            // try to infer what resolutions the user wants
            if (!(tc_lp1.isSet() || tc_lp2.isSet() || tc_lp3.isSet())) {
                if (lpmm_mode) {
                    // if nothing is specified explicitly, use the first two defaults
                    resolutions.push_back(tc_lp1.getValue());
                    resolutions.push_back(tc_lp2.getValue());
                } else {
                    // otherwise just pick some arbitrary values
                    resolutions.push_back(0.1);
                    resolutions.push_back(0.2);
                }
            } else {
                if (tc_lp1.isSet()) {
                    resolutions.push_back(tc_lp1.getValue());
                }
                if (tc_lp2.isSet()) {
                    resolutions.push_back(tc_lp2.getValue());
                }
                if (tc_lp3.isSet()) {
                    resolutions.push_back(tc_lp3.getValue());
                }
            }

            sort(resolutions.begin(), resolutions.end());
            
            Mtf_renderer_lensprofile printer(
                img_filename,
                wdir, 
                string("lensprofile.txt"),
                tc_gnuplot.getValue(),
                cvimg,
                resolutions,
                gnuplot_width,
                lpmm_mode,
                pixel_size
            );
            printer.set_sparse_chart(tc_ima_mode.getValue());
            printer.set_fixed_size(tc_lensprofile_fixed.getValue());
            printer.render(mtf_core.get_blocks());
        }
        
        if (tc_chart_orientation.getValue()) {
            Mtf_renderer_chart_orientation co_renderer(
                img_filename,
                wdir, 
                string("chart_orientation.png"),
                cvimg,
                gnuplot_width,
                distance_scale,
                &img_dimension_correction
            );
            co_renderer.render(mtf_core.get_blocks());
        }

        if (tc_sfr.getValue() || tc_absolute.getValue()) {
            Mtf_renderer_sfr sfr_writer(
                wdir + string("raw_sfr_values.txt"), 
                lpmm_mode,
                pixel_size
            );
            sfr_writer.render(mtf_core.get_blocks());
        }
        
        if (tc_esf.getValue()) {
            Output_version::type ver = Output_version::type(tc_output_version.getValue());
            Mtf_renderer_esf esf_writer(
                wdir + string("raw_esf_values.txt"), 
                wdir + (ver >= Output_version::V2 ? string("raw_lsf_values.txt") : string("raw_psf_values.txt")),
                ver
            );
            esf_writer.render(mtf_core.get_blocks());
        }
        
        if (tc_focus.getValue()) {
            Mtf_renderer_focus profile(
                distance_scale,
                mtf_core.get_sliding_edges(),
                wdir, 
                string("focus_peak.png"),
                cvimg,
                lpmm_mode,
                pixel_size
            );
            profile.render(mtf_core.get_samples(), mtf_core.bayer, &mtf_core.ellipses, &img_dimension_correction);
        }
        
        Mtf_renderer_stats stats(lpmm_mode, pixel_size);
        if (tc_focus.getValue() || tc_mf_profile.getValue()) {
            stats.render(mtf_core.get_samples());
        } else {
            stats.render(mtf_core.get_blocks());
        }
        
        if (tc_ca.getValue()) {
            Ca_renderer_print ca_print(wdir + string("chromatic_aberration.txt"), cvimg, tc_ca_all.getValue());
            ca_print.render(mtf_core.get_blocks());
            
            Ca_renderer_grid grid(
                img_filename,
                wdir, 
                string("ca_grid.txt"),
                tc_gnuplot.getValue(),
                cvimg,
                gnuplot_width,
                lpmm_mode,
                pixel_size,
                tc_ca_fraction.getValue(),
                tc_ca_all.getValue()
            );
            grid.set_sparse_chart(tc_ima_mode.getValue());
            grid.render(mtf_core.get_blocks());
        }
        
    } while (!finished);
    
    return 0;
}
