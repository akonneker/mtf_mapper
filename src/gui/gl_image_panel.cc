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
#include "gl_image_panel.h"
#include <QOpenGLShaderProgram>
#include <QOpenGLTexture>
#include <QMouseEvent>
#include <opencv2/imgcodecs/imgcodecs.hpp>

#include <cmath>

int GL_image_panel::program_counter = 0;

GL_image_panel::GL_image_panel(QWidget *parent)
    : QOpenGLWidget(parent),
      program(0) {
      
    logger.debug("GL_image_panel ctor: OpenGL version: %d.%d, samples=%d\n", format().majorVersion(), format().minorVersion(), format().samples());
}

GL_image_panel::~GL_image_panel() {
    makeCurrent();
    vbo.destroy();
    for (int i=0; i < (int)textures.size(); i++) {
        delete textures[i];
    }
    delete program;
    doneCurrent();
}

QSize GL_image_panel::minimumSizeHint() const {
    return QSize(50, 50);
}

QSize GL_image_panel::sizeHint() const {
    return QSize(800, 800);
}

void GL_image_panel::initializeGL() {
    initializeOpenGLFunctions();
    static int panel_counter = 0;
    logger.debug("initializeGL (panel = %d): OpenGL version: %d.%d, samples=%d\n", ++panel_counter, format().majorVersion(), format().minorVersion(), format().samples());
    
    if (textures.size() == 0 && default_image != nullptr) {
        load_image(*default_image);
    }

    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_MULTISAMPLE);

    // shader program for main texture image
    QOpenGLShader *vshader = new QOpenGLShader(QOpenGLShader::Vertex, this);
    const char *vsrc =
        "#version 120\n"
        "attribute highp vec4 vertex;\n"
        "attribute highp vec4 texCoord;\n"
        "varying highp vec4 texc;\n"
        "uniform mediump mat4 projectionMatrix;\n"
        "uniform mediump mat4 modelMatrix;\n"
        "void main(void)\n"
        "{\n"
        "    gl_Position = projectionMatrix * modelMatrix * vertex;\n"
        "    texc = texCoord;\n"
        "}\n";
    vshader->compileSourceCode(vsrc);

    QOpenGLShader *fshader = new QOpenGLShader(QOpenGLShader::Fragment, this);
    const char *fsrc =
        "#version 120\n"
        "uniform sampler2D texture;\n"
        "varying highp vec4 texc;\n"
        "uniform float gamma;\n"
        "void main(void)\n"
        "{\n"
        "    gl_FragColor =  pow(texture2D(texture, texc.st), vec4(1.0 / gamma, 1.0 / gamma, 1.0 / gamma, 0.0));\n"
        "}\n";
    fshader->compileSourceCode(fsrc);

    program = new QOpenGLShaderProgram;
    program->addShader(vshader);
    program->addShader(fshader);
    program->bindAttributeLocation("vertex", prog_vert_att);
    program->bindAttributeLocation("texCoord", prog_texcoord_att);
    program->link();

    program->bind();
    program->setUniformValue("texture", 0);
    
    // call any initialization code in derived classes
    initialize_overlay();
}

void GL_image_panel::paintGL() {
    
    int w = size().width();
    int h = size().height();
    
    reset_scroll_range();
    
    glClearColor(1,1,1,0);
    glClear(GL_COLOR_BUFFER_BIT);
    
    view = QMatrix4x4();
    view.translate(-vp.centre.x + int(w/2), -vp.centre.y + int(h/2));
    view.scale(scale_factor);
    
    projection = QMatrix4x4();
    projection.ortho(0, w, h, 0, -1, 1);
    
    program->bind();
    program->setUniformValue("modelMatrix", view);
    program->setUniformValue("projectionMatrix", projection);
    
    vbo.bind();
    program->enableAttributeArray(prog_vert_att);
    program->enableAttributeArray(prog_texcoord_att);
    program->setAttributeBuffer(prog_vert_att, GL_FLOAT, 0, 3, 5 * sizeof(GLfloat));
    program->setAttributeBuffer(prog_texcoord_att, GL_FLOAT, 3 * sizeof(GLfloat), 2, 5 * sizeof(GLfloat));
    
    program->setUniformValue("gamma", float(gamma_value));
    
    for (int i = 0; i < (int)textures.size(); i++) {
        textures[i]->bind();
        glDrawArrays(GL_TRIANGLE_FAN, i * 4, 4);
    }
    
    // call any paint code in derived classes
    paint_overlay();
}


void GL_image_panel::resizeGL([[maybe_unused]] int width, [[maybe_unused]] int height) {
    double w = size().width();
    double h = size().height();
    if (imgsize.width() > w || imgsize.height() > h) {
        scale_factor = std::min(w / imgsize.width(), h / imgsize.height());
    }
    
    reset_scroll_range();
    
    update();
}

static int next_pow2(uint32_t x) {
    uint32_t k=1;
    while (k < 31 && (uint32_t(1) << k) < x) {
        k++;
    }
    return uint32_t(1) << k;
}

bool GL_image_panel::load_image(const QString& fname) {
    if (cache_enabled && fname.toStdString().compare(current_fname) == 0) {
        return true;
    }
    
    cv::Mat cvimg;
    if (cache_enabled && image_cache.find(fname.toStdString()) != image_cache.end()) {
        cvimg = image_cache[fname.toStdString()].fetch();
    } else { 
        cvimg = cv::imread(fname.toStdString(), cv::IMREAD_ANYCOLOR | cv::IMREAD_ANYDEPTH | cv::IMREAD_IGNORE_ORIENTATION);
        if (cvimg.empty()) {
            return false;
        }
        
        img_channels = cvimg.channels();
        img_depth = 8;
        img_assumed_linear = false;
        
        // perform normalization / scaling before color conversion?
        if (cvimg.depth() != CV_8U) {
            double alpha = 1.0;
            double beta = 0;
            img_depth = 16;
            img_assumed_linear = true;
            
            if (cvimg.depth() == CV_16U) {
                logger.info("GL_image_panel::load_image: 16-bit histogram percentile scaling path\n");
                
                vector<uint64_t> histo(65536, 0);
                size_t total_pixels = size_t(cvimg.channels())*size_t(cvimg.rows)*size_t(cvimg.cols);
                uint16_t* ptr = (uint16_t*)cvimg.data;
                uint16_t* sentinel = ptr + total_pixels;
                while (ptr < sentinel) {
                    histo[*ptr++]++;
                }
                // now we can extract the true (clipping) max and min values
                int minval = histo.size();
                int maxval = 0;
                for (size_t i=0; i < histo.size(); i++) {
                    int revi = histo.size() - 1 - i;
                    minval = histo[revi] > 0 ? revi : minval;
                    maxval = histo[i] > 0 ? i : maxval;
                }
                
                // find upper percentile
                size_t upper_c_target = total_pixels * 0.01;
                uint64_t acc = 0;
                size_t hist_upper = maxval;
                while (hist_upper > 0 && acc < upper_c_target) {
                    acc += histo[hist_upper];
                    hist_upper--;
                }
                
                // find lower percentile
                size_t lower_c_target = total_pixels * 0.005;
                acc = 0;
                size_t hist_lower = minval;
                while (hist_lower < size_t(maxval) && acc < lower_c_target) {
                    acc += histo[hist_lower];
                    hist_lower++;
                }
                
                // prevent complete collapse on synthetic/degenerate images
                if (fabs(double(hist_upper) - double(hist_lower)) < 0.05*(maxval - minval)) {
                    hist_lower = minval;
                    hist_upper = maxval;
                }
                
                // scale intenisty so that minval -> 1 and maxval -> 254
                if (hist_upper > hist_lower) {
                    alpha = 253.0/double(hist_upper - hist_lower);
                    beta = 1.0 - alpha * hist_lower;
                }
                
                vector<uint8_t> lut(65536, 0);
                size_t crush_thresh = std::min(2048, std::min(minval + 1, (int)hist_lower));
                size_t saturation_thresh = maxval;
                size_t pot = 1;
                while (pot < 16 && (1 << pot) < maxval) {
                    pot++;
                }
                saturation_thresh = std::max((1 << pot) - 1, (int)hist_upper);
                
                size_t lut_idx = 0;
                while (lut_idx < crush_thresh) {
                    lut[lut_idx++] = 0;
                }
                while (lut_idx < hist_lower) {
                    lut[lut_idx++] = 1;
                }
                while (lut_idx < hist_upper) {
                    lut[lut_idx] = std::max(1, std::min(254, int(lut_idx*alpha + beta)));
                    lut_idx++;
                }
                while (lut_idx < saturation_thresh) {
                    lut[lut_idx++] = 254;
                }
                while (lut_idx < lut.size()) {
                    lut[lut_idx++] = 255;
                }
                
                cv::Mat dst(cvimg.rows, cvimg.cols, cvimg.channels() == 1 ? CV_8UC1 : CV_8UC3);
                uint16_t* in_ptr = (uint16_t*)cvimg.data;
                uint8_t* out_ptr = (uint8_t*)dst.data;
                while (in_ptr < sentinel) {
                    *out_ptr++ = lut[*in_ptr++];
                }
                cvimg = dst;
                
                img_minmax.first = minval;
                img_minmax.second = maxval;
            } else {
                logger.info("GL_image_panel::load_image: unspecified (not 8- or 16-bit) min-max scaling path\n");
                double minval, maxval;
                cv::minMaxIdx(cvimg.reshape(1), &minval, &maxval);
                
                if (maxval > minval) {
                    alpha = 255.0/(maxval - minval);
                    beta = -alpha * minval;
                }
                
                img_minmax.first = minval;
                img_minmax.second = maxval;
                
                cv::Mat dst;
                cvimg.convertTo(dst, CV_8U, alpha, beta);
                cvimg = dst;
            }
        } else {
            double minval, maxval;
            cv::minMaxIdx(cvimg.reshape(1), &minval, &maxval);
            img_minmax.first = minval;
            img_minmax.second = maxval;
        }
        
        if (cvimg.channels() == 1) {
            cv::cvtColor(cvimg, cvimg, cv::COLOR_GRAY2BGR);
        }
        
        // manual conversion, cvtColor seems a bit slow
        uint8_t* sptr = cvimg.data;
        uint8_t* sentinel = sptr + size_t(cvimg.rows)*size_t(cvimg.cols*3);
        while (sptr < sentinel) {
            std::swap(*sptr, *(sptr+2));
            sptr += 3;
        }
        
        if (cache_enabled) {
            image_cache[fname.toStdString()] = Cache_entry(cvimg);
            image_cache_size += image_cache[fname.toStdString()].size();
        }
        
        trim_cache();
    }
    
    current_fname = fname.toStdString();
    
    return load_image(cvimg);
}

bool GL_image_panel::load_image(QImage& qimg) {
    cv::Mat cvimg(qimg.height(), qimg.width(), CV_8UC3, qimg.bits());
    current_fname = "mtf_mapper_logo";
    return load_image(cvimg);
}

bool GL_image_panel::load_image(cv::Mat cvimg) {
    // Some other GLWidget could have the current context, so it is vital that
    // we switch the context before we attempt to change resources like textures
    makeCurrent();
    
    current_cv_image = cvimg;
    
    bool keep_zoom = false;

    vbo.release();
    vbo.destroy();
    for (int i=0; i < (int)textures.size(); i++) {
        delete textures[i];
    }
    textures.clear();
    
    static GLsizei hw_texw = 0;
    if (hw_texw == 0) { // we only really have to call this once 
        glGetIntegerv(GL_MAX_TEXTURE_SIZE, &hw_texw);
    }
    
    int texw = std::min(hw_texw, 2048); // we do not want the last block to become too large
    
    if (cvimg.cols == imgsize.width() && cvimg.rows == imgsize.height()) {
        keep_zoom = true;
    }
    
    imgsize = QSize(cvimg.cols, cvimg.rows);
    
    // now trim down the texture block size if the image is small
    texw = std::min(texw, std::max(next_pow2(imgsize.width()), next_pow2(imgsize.height())));
    
    QVector<GLfloat> vertData;
    
    int cblocks = (int)ceil(imgsize.width()/double(texw));
    int rblocks = (int)ceil(imgsize.height()/double(texw));
    
    int x_off = imgsize.width() / 2;
    int y_off = imgsize.height() / 2;
    
    for (int r=0; r < rblocks; r++) {
        int rstart = r*texw;
        int rend = std::min((r+1)*texw, imgsize.height());
        for (int c=0; c < cblocks; c++) {
            int cstart = c*texw;
            int cend = std::min((c+1)*texw, imgsize.width());
            
            // obtain pow2 dimensions for this block
            int tex_width = next_pow2(cend-cstart);
            int tex_height = next_pow2(rend-rstart);
            
            vertData.append(cstart - x_off);
            vertData.append(rstart - y_off);
            vertData.append(0);
            vertData.append(0);
            vertData.append(0);
            
            vertData.append(cend - x_off);
            vertData.append(rstart - y_off);
            vertData.append(0);
            vertData.append(double(cend - cstart)/tex_width);
            vertData.append(0);
            
            vertData.append(cend - x_off);
            vertData.append(rend - y_off);
            vertData.append(0);
            vertData.append(double(cend - cstart)/tex_width);
            vertData.append(double(rend - rstart)/tex_height);
            
            vertData.append(cstart - x_off);
            vertData.append(rend - y_off);
            vertData.append(0);
            vertData.append(0);
            vertData.append(double(rend - rstart)/tex_height);
            
            // we copy parts of the cv::Mat image into a QImage
            // this way we can have large images (>10k dims were problematic QT 5.7), and
            // at the same time have auto-mipmaps (which we only get from QImage)
            if ((cend - cstart) == tex_width && (rend - rstart) == tex_height) {
                size_t offset = cstart * 3 + rstart * cvimg.elemSize()*cvimg.cols;
                QImage tex_block(cvimg.data + offset, tex_width, tex_height, cvimg.elemSize()*cvimg.cols, QImage::Format_RGB888);
                textures.push_back(new QOpenGLTexture(tex_block));
            } else {
                // only create a new image (power-of-2) if we have a partial block
                QImage tex_block(tex_width, tex_height, QImage::Format_RGB888);
                cv::Mat srcimg(cvimg, cv::Rect(cstart, rstart, cend-cstart, rend-rstart)); // the ROI we want
                cv::Mat teximg(tex_height, tex_width, CV_8UC3, tex_block.bits());
                teximg = cv::Scalar(255, 255, 255); // fill the image to avoid borders (during mipmap building)
                srcimg.copyTo(teximg(cv::Rect(0, 0, cend-cstart, rend-rstart)));
                textures.push_back(new QOpenGLTexture(tex_block));
            }
            
            // but this looks better (other than the artifacts, of course)
            textures.back()->setMinificationFilter(QOpenGLTexture::LinearMipMapLinear);
            textures.back()->setMagnificationFilter(QOpenGLTexture::Linear);
            textures.back()->setWrapMode(QOpenGLTexture::ClampToEdge);
        }
    }
    
    double w = size().width();
    double h = size().height();
    if (!keep_zoom) {
        scale_factor = std::min(1.0, std::min(w / imgsize.width(), h / imgsize.height()));
        
        // reset centre position (in image)
        vp.centre = cv::Point2d(0,0);
    }
    
    reset_scroll_range();

    vbo.create();
    vbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);
    vbo.bind();
    vbo.allocate(vertData.constData(), vertData.count() * sizeof(GLfloat));
    update();
    return true;
}

void GL_image_panel::move(int nx, int ny) {
    vp.centre = cv::Point2d(nx, ny);
}

QPoint GL_image_panel::zoom(int step, int mx, int my) {
    double iw = imgsize.width();
    double ih = imgsize.height();
    int w = size().width();
    int h = size().height();
    
    // scale to next power of two, but skip to the next larger/smaller scale
    // if the relative change in scale is less than 25% when near the minimum scale
    double min_scale_factor = std::min(1.0, std::min(w/iw, h/ih));
    double next_scale_factor = scale_factor;
    if (step < 0) {
        double p = floor(log(scale_factor*0.5)/log(2.0));
        next_scale_factor = std::min(max_scale_factor, std::max(min_scale_factor, pow(2.0, p)));
        if (fabs(next_scale_factor - min_scale_factor)/min_scale_factor < 0.25) {
            next_scale_factor = min_scale_factor;
        }
    } else {
        if (step >= 0) {
            double p = floor(log(scale_factor*2)/log(2.0));
            next_scale_factor = std::max(min_scale_factor, std::min(max_scale_factor, pow(2.0, p)));
            if (fabs(next_scale_factor - min_scale_factor)/min_scale_factor < 0.25) {
                next_scale_factor = std::max(min_scale_factor, std::min(max_scale_factor, pow(2.0, p+1)));
            }
        }
    }
    scale_factor = next_scale_factor;
    
    QPointF m(mx, my);
    QPointF mp = view.inverted().map(m);
    
    QMatrix4x4 nview;
    nview.translate(-vp.centre.x + w/2, -vp.centre.y + h/2);
    nview.scale(scale_factor);
    
    QPointF mpi = nview.map(mp);
    QPointF delta = m - mpi;
    
    vp.centre.x -= delta.x();
    vp.centre.y -= delta.y();
    
    reset_scroll_range();
    
    return QPoint(vp.centre.x, vp.centre.y);
}

QPoint GL_image_panel::locate(QPoint pos) {
    double iw = imgsize.width();
    double ih = imgsize.height();
    
    QPointF p = view.inverted().map(QPointF(pos));
    double ix = p.x() + iw/2;
    double iy = p.y() + ih/2;
    
    located_pos = QPoint(ix, iy);
    
    return located_pos;
}

// A bit of a hack, but this is required to ensure that small images remain
// centered. Resetting this is necessary after scale changes, or window
// resizing, or loading a new image
void GL_image_panel::reset_scroll_range(void) {
    int w = size().width();
    int h = size().height();
    double iw = imgsize.width();
    double ih = imgsize.height();
    
    vp.set_x_range(-iw/2*scale_factor + w/2, iw/2*scale_factor - w/2);
    vp.set_y_range(-ih/2*scale_factor + h/2, ih/2*scale_factor - h/2);
}

void GL_image_panel::set_cache_size(uint64_t size) { 
    max_image_cache_size = size;
    trim_cache();
}

void GL_image_panel::trim_cache(void) {
    if (image_cache_size > max_image_cache_size) {
        // evict entries if necessary (we cannot evict pre-emptively because 
        // image dimensions are unknown until after imread())
        vector< pair<uint32_t, string> > cache_age;
        for (const auto& e: image_cache) {
            cache_age.push_back(make_pair(e.second.seq(), e.first));
        }
        sort(cache_age.begin(), cache_age.end());
        
        
        size_t ci = 0;
        while (image_cache_size > max_image_cache_size && ci < cache_age.size()) {
            auto cit = image_cache.find(cache_age[ci].second);
            image_cache_size -= cit->second.size();
            image_cache.erase(cit);
            ci++;
        }
    }
}

void GL_image_panel::mouseMoveEvent(QMouseEvent* event) {
    event->ignore();
}

void GL_image_panel::mousePressEvent(QMouseEvent* event) {
    event->ignore();
}

void GL_image_panel::mouseReleaseEvent(QMouseEvent* event) {
    event->ignore();
}



