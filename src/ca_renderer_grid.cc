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
#include "include/mtf_renderer_grid.h"
#include "include/ca_renderer_grid.h"
#include "include/grid_interpolator.h"
#include <memory>

class Grid_functor_ca : public Grid_functor {
  public:
    Grid_functor_ca(bool red_ca, double scale_factor = 1.0, bool sparse_mode = false) 
    : red_ca(red_ca), scale_factor(scale_factor), sparse_mode(sparse_mode) {}
    
    virtual ~Grid_functor_ca(void) {}
    
    void set_red_ca(bool val) {
        red_ca = val;
    }
    
    inline bool get_red_ca(void) const {
        return red_ca;
    }
    
    virtual double value(const Block& block, size_t edge) const {
        return (get_red_ca() ? block.get_ca(edge).x : block.get_ca(edge).y) * scale_factor;
    }
    
    virtual bool in_range(double val) const {
        return val > -16.0*scale_factor && val < 16.0*scale_factor; 
    }
    
    virtual double clamp(double value, double upper, double lower) const {
        return value < lower ? lower : (value > upper ? upper : value);
    }
    
    virtual double nodata(void) const {
        return -100.0*scale_factor;
    }
    
    double scale(double x) {
        return x*scale_factor;
    }

    virtual double smoothing_factor(void) const {
        return sparse_mode ? 1e-1 : 1e-2;
    }

    virtual int pruning_threshold(void) const {
        return sparse_mode ? 3 : 2;
    }

  protected:
    bool red_ca = true;
    double scale_factor = 1.0;
    bool sparse_mode = false;
};

class Grid_functor_ca_radial : public Grid_functor_ca {
  public:
    Grid_functor_ca_radial(bool red_ca, cv::Size img_dims, bool sparse_mode) 
    : Grid_functor_ca(red_ca, 1.0, sparse_mode), 
      centre(img_dims.width/2.0, img_dims.height/2.0) {}
      
    virtual double value(const Block& block, size_t edge) const {
        double radial_distance = cv::norm(block.get_edge_centroid(edge) - centre);
        double radial_scale = 100.0/std::max(100.0, radial_distance);
        return (get_red_ca() ? block.get_ca(edge).x : block.get_ca(edge).y) * radial_scale;
    }
    
  private:
    cv::Point2d centre;
};

Ca_renderer_grid::Ca_renderer_grid(
    const std::string& img_filename,
    const std::string& wdir, 
    const std::string& fname, 
    const std::string& gnuplot_binary, 
    const cv::Mat& img, 
    int gnuplot_width,
    bool lpmm_mode, double pixel_size,
    bool fraction_mode, 
    bool allow_all_edges)
    : img_filename(img_filename),
      wdir(wdir), fname(fname), 
      gnuplot_binary(gnuplot_binary), 
      img(img), gnuplot_width(gnuplot_width),
      lpmm_mode(lpmm_mode), pixel_size(pixel_size),
      fraction_mode(fraction_mode),
      img_centre(img.cols/2, img.rows/2),
      img_dims(img.cols, img.rows),
      allow_all_edges(allow_all_edges) {

    const int coarse_grid_size = 40;
    const int fine_grid_size = 200;
    
    size_t grid_y_coarse;
    size_t grid_x_coarse;
    size_t grid_y_fine;
    size_t grid_x_fine;
    
    if (img.rows > img.cols) {
        grid_y_coarse = coarse_grid_size;
        grid_x_coarse = coarse_grid_size * img.cols / img.rows;
        grid_y_fine = fine_grid_size;
        grid_x_fine = fine_grid_size * img.cols / img.rows;
    } else {
        grid_x_coarse = coarse_grid_size;
        grid_y_coarse = coarse_grid_size * img.rows / img.cols;
        grid_x_fine = fine_grid_size;
        grid_y_fine = fine_grid_size * img.rows / img.cols;
    }
    
    for (int i=0; i < 2; i++) {
        grid_coarse.push_back(cv::Mat(grid_y_coarse, grid_x_coarse, CV_32FC1, 0.0));
        grid_fine.push_back(cv::Mat(grid_y_fine, grid_x_fine, CV_32FC1, 0.0));
    }
}

void Ca_renderer_grid::render(const vector<Block>& blocks) {

    if (blocks.size() < 6) {
        logger.error("%s\n", "Too few valid blocks found. No CA surface can be generated");
        return;
    }
    
    std::unique_ptr<Grid_functor_ca> ca_ftor; 
    
    if (fraction_mode) {
        ca_ftor = std::unique_ptr<Grid_functor_ca>(new Grid_functor_ca_radial(true, img_dims, sparse_mode));
    } else {
        ca_ftor = std::unique_ptr<Grid_functor_ca>(new Grid_functor_ca(true, lpmm_mode ? 1000.0/pixel_size : 1, sparse_mode));
    }
    
    interpolate_grid(*ca_ftor, allow_all_edges ? NEITHER : MERIDIONAL, grid_coarse[0], grid_fine[0], img_dims, blocks, ca_ftor->scale(-200), ca_ftor->smoothing_factor(), ca_ftor->pruning_threshold());
    ca_ftor->set_red_ca(false);
    interpolate_grid(*ca_ftor, allow_all_edges ? NEITHER : MERIDIONAL, grid_coarse[1], grid_fine[1], img_dims, blocks, ca_ftor->scale(-200), ca_ftor->smoothing_factor(), ca_ftor->pruning_threshold());
    
    vector<double> zmax(2, -1e50);
    vector<double> zmin(2,  1e50);
    
    FILE* file = fopen((wdir+fname).c_str(), "wt");
    for (size_t k=0; k < 2; k++) {
        fprintf(file, "#coarse %s grid\n", k == 0 ? "red" : "blue");
        for (int y=0; y < grid_coarse[k].rows; y++) {
            for (int x=0; x < grid_coarse[k].cols; x++) {
                fprintf(file, "%lf %lf %.5lf\n", 
                    x*img_dims.width/grid_coarse[k].cols/pixel_size, 
                    y*img_dims.height/grid_coarse[k].rows/pixel_size, 
                    grid_coarse[k].at<float>(y,x)
                );
                zmax[k] = std::max(zmax[k], (double)grid_coarse[k].at<float>(y, x));
                zmin[k] = std::min(zmin[k], (double)grid_coarse[k].at<float>(y, x));
            }
            fprintf(file, "\n");
        }
        fprintf(file, "\n\n");
    }

    #if 0
    // we no longer use the fine grid, but rather let gnuplot interpolate on its own
    for (size_t k=0; k < 2; k++) {
        fprintf(file, "#fine %s grid\n", k == 0 ? "red" : "blue");
        for (int y=0; y < grid_fine[k].rows; y++) {
            for (int x=0; x < grid_fine[k].cols; x++) {
                fprintf(file, "%lf %lf %.5lf\n", 
                    x*img_dims.width/grid_fine[k].cols/pixel_size, 
                    y*img_dims.height/grid_fine[k].rows/pixel_size, 
                    grid_fine[k].at<float>(y,x)
                );
                zmax[k] = std::max(zmax[k], (double)grid_fine[k].at<float>(y,x));
                zmin[k] = std::min(zmin[k], (double)grid_fine[k].at<float>(y,x));
            }
            fprintf(file, "\n");
        }
        fprintf(file, "\n\n");
    }
    #endif

    fclose(file);
    
    // ensure zero is included in the range for each of R/B
    for (size_t k=0; k < 2; k++) {
        if (zmax[k]*zmin[k] > 0) {
            if (zmax[k] < 0) {
                zmax[k] = 0;
            } else {
                zmin[k] = 0;
            }
        }
    }
    
    double g_zmax = std::max(zmax[0], zmax[1]);
    double g_zmin = std::min(zmin[0], zmin[1]);

    if (g_zmax < 0.1 && g_zmin < 0 && g_zmin > -0.1) {
        // We are probably dealing with a synthetic grayscale image, so there is no real CA
        // Bump the limits a bit to produce a more pleasing plot
        g_zmax = 0.1;
        g_zmin = -0.1;
    }
    
    const int width_in_pixels = int(gnuplot_width*600.0/1024);
    int fontsize = std::max(long(10), lrint(10.0*gnuplot_width/1024.0));
    int title_fontsize = fontsize + 2;
    
    FILE* gpf = fopen((wdir + std::string("ca_grid.gnuplot")).c_str(), "wt");
    fprintf(gpf, "%s\n", diverging_palette.c_str());
    fprintf(gpf, "set cbrange [%lf:%lf]\n", g_zmin, g_zmax);
    fprintf(gpf, "set xlab \"column (%s)\"\n", lpmm_mode ? "mm" : "pixels");
    fprintf(gpf, "set ylab \"row (%s)\"\n",  lpmm_mode ? "mm" : "pixels");
    fprintf(gpf, "set xtics out nomirror\n");
    fprintf(gpf, "set xtics scale 0.3\n");
    fprintf(gpf, "set ytics out nomirror\n");
    fprintf(gpf, "set ytics scale 0.3\n");
    if (!lpmm_mode) {
        fprintf(gpf, "set ytics rotate by 45\n");
    }
    fprintf(gpf, "set pm3d map impl interpolate 3,3\n");
    fprintf(gpf, "set cntrlabel onecolor\n");
    fprintf(gpf, "set contour surface\n");
    fprintf(gpf, "set cntrparam order 8\n");
    fprintf(gpf, "set cntrparam bspline\n");
    fprintf(gpf, "set hidden3d\n"); // we need this hidden3d option to hide the cell boundary lines
    fprintf(gpf, "set autoscale xfix\n");
    fprintf(gpf, "set autoscale yfix\n");
    fprintf(gpf, "\n");
    fprintf(gpf, "set term pngcairo dashed transparent enhanced size %d, %d font '%s,%d'  background rgb \"white\"\n",
        width_in_pixels, 
        (int)lrint(width_in_pixels*2.3*grid_fine[0].rows/double(grid_fine[0].cols)), 
        #ifdef _WIN32
        "Verdana",
        #else
        "Arial",
        #endif
        fontsize
    );
    
    fprintf(gpf, "set output \"%sca_image.png\"\n", wdir.c_str());
    if (img_filename.length() > 0) {
        fprintf(gpf, "set multiplot title \"%s\" font \",%d\"\n", img_filename.c_str(), title_fontsize);
        fprintf(gpf, "set tmargin 4\n");
    } else {
        fprintf(gpf, "set multiplot\n");
    }
    
    fprintf(gpf, "set size 1,0.5\n");
    fprintf(gpf, "set origin 0.0,0.5\n");
    fprintf(gpf, "set yrange [*:*] reverse\n");
    if (!fraction_mode) {
        fprintf(gpf, "set title \"Red vs Green (shift in %s)\"\n", lpmm_mode ? "micron" : "pixels");
    } else {
        fprintf(gpf, "set title \"Red vs Green (%% of radial distance)\"\n");
    }
    fprintf(gpf, "splot \"%s\" i 0 notitle w l lc rgb \"#77303030\"\n", (wdir+fname).c_str());
    fprintf(gpf, "set origin 0.0,0.0\n");
    if (!fraction_mode) {
        fprintf(gpf, "set title \"Blue vs Green (shift in %s)\"\n", lpmm_mode ? "micron" : "pixels");
    } else {
        fprintf(gpf, "set title \"Blue vs Green (%% of radial distance)\"\n");
    }
    fprintf(gpf, "splot \"%s\" i 1 notitle w l lc rgb \"#77303030\"\n",  (wdir+fname).c_str());
    fprintf(gpf, "unset multiplot\n");
    
    
    fclose(gpf);
    
    char* buffer = new char[1024];
    #ifdef _WIN32
    sprintf(buffer, "\"\"%s\" \"%sca_grid.gnuplot\"\"", gnuplot_binary.c_str(), wdir.c_str());
    #else
    sprintf(buffer, "\"%s\" \"%sca_grid.gnuplot\"", gnuplot_binary.c_str(), wdir.c_str());
    #endif
    int rval = system(buffer);
    if (rval != 0) {
        logger.error("Failed to execute gnuplot (error code %d)\n", rval);
        logger.info("You can try to execute [%s] to render the plots manually\n", buffer);
    } else {
        logger.debug("%s\n", "Gnuplot plot completed successfully. Look for grid_image.png and grid_surface.png");
    }
    
    delete [] buffer;
}

