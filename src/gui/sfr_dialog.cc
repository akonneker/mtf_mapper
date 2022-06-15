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
#include <QtWidgets> 
#include "sfr_dialog.h"

#include "common.h"
#include "config.h"
#include <locale.h>

#include <QtGlobal>
#if QT_VERSION < QT_VERSION_CHECK(5,11,0)
#define horizontalAdvance width
#endif

Sfr_dialog::Sfr_dialog(QWidget* parent ATTRIBUTE_UNUSED, const QRect initial_geom) 
: QDialog(parent, Qt::WindowSystemMenuHint | Qt::WindowTitleHint | Qt::WindowCloseButtonHint),
  cursor_domain_value(0), repainting(0) {
  
    if (initial_geom != QRect(-1, -1, 0, 0)) {
        setGeometry(initial_geom);
    }
    
    chart = new QChart();
    chart->legend()->hide();
    
    chart->setAnimationOptions(QChart::NoAnimation);
    
    x_axis = new QValueAxis();
    chart->addAxis(x_axis, Qt::AlignBottom);
    
    y_axis = new QValueAxis();
    chart->addAxis(y_axis, Qt::AlignLeft);
    
    x_scale_levels = {2,4,6,16,32};
    x_scale_slider = new QSlider(Qt::Horizontal);
    x_scale_slider->setTickPosition(QSlider::TicksBelow);
    x_scale_slider->setMinimum(0);
    x_scale_slider->setMaximum(x_scale_levels.size() - 1);
    x_scale_slider->setValue(x_scale_levels.size() - 2); // start at the equivalent of 16 pixels
    x_scale_slider->setTickInterval(1);
    
    lp_mm_cb = new QCheckBox(this);
    update_lp_mm_mode();
    lp_mm_cb->setCheckState(view.get_lp_mm_mode() ? Qt::Checked : Qt::Unchecked);
    
    chart_view = new Sfr_chartview(chart, this);
    chart_view->setRenderHint(QPainter::Antialiasing);
    
    mtfmapper_logo = std::shared_ptr<QIcon>(new QIcon);
    mtfmapper_logo->addFile(":/Icons/AppIcon256");
    
    QFontMetrics fm(QWidget::fontMetrics());
    int double_height = fm.height()*2 + 8;
    int em_width = fm.horizontalAdvance("m");
    
    QString button_style(
        "QPushButton:flat{ background-color: rgba(0, 0, 0, 7%); border: 1px solid rgba(0,0,0,20%); border-radius: 2px; } "
        "QPushButton:flat::pressed{ background-color: rgba(0, 0, 0, 15%); border: 1px solid rgba(0,0,0,20%); border-radius: 2px; } "
    );
    
    save_img_button = new QPushButton("Save\nimage", this);
    save_img_button->setMinimumWidth(fm.horizontalAdvance(save_img_button->text()));
    save_img_button->setMaximumWidth(fm.horizontalAdvance(save_img_button->text()));
    save_img_button->setMinimumHeight(double_height);
    save_img_button->setMaximumHeight(double_height);
    save_img_button->setFlat(true);
    save_img_button->setStyleSheet(button_style);
    
    QFont fixed_font = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    fixed_font.setStyleHint(QFont::TypeWriter);
    QFont default_font = QFontDatabase::systemFont(QFontDatabase::GeneralFont);
    default_font.setBold(true);
    QFontMetrics default_fm(default_font);
    vector<string> column_names = {" contrast ", " CNR ", " red-CA ", " blue-CA ", " angle ", " OSF ", " length ", " chan. "};
    vector<string> column_tooltips = {
        "",
        "CNR = Contrast-to-Noise Ratio\nUnits: unitless ratio\nHigher is better, > 30 is Ok, > 100 is excellent",
        "red-CA = Lateral Chromatic Aberration, shift between red channel and green channel\nUnits: pixels in c/p mode, micron in lp/mm mode\nSmaller magnitude is better",
        "blue-CA = Lateral Chromatic Aberration, shift between blue channel and green channel\nUnits: pixels in c/p mode, micron in lp/mm mode\nSmaller magnitude is better",
        "angle = Relative edge orientation angle (modulo 45 degrees)\nUnits: degrees\nValues < 1, close to 26.565, and > 44 are undesirable",
        "OSF = Over-Sampling Factor, measures how reliable a measurement is\nUnits: unitless factor\nValues <= 4 are undesirable",
        "length = Edge length, excluding trimmed zones near corners\nUnits: pixels\nValues <= 25 are undesirable, about 100 to 400 is considered exellent, > 400 is of questionable value",
        "chan. = Channel type\nL = luminance of RGB or Gray image\nR, G, B = Red, Green, Blue channel"
    };
    int min_col_width = default_fm.boundingRect(column_names[3].c_str()).width(); // use "blue-CA" as the minimum
    
    vector<QVBoxLayout*> table_col_layouts;
    QHBoxLayout* table_row_layouts = new QHBoxLayout;
    
    table_labels.resize(4);
    for (size_t r=0; r < table_labels.size(); r++) {
        table_labels[r].resize(column_names.size());
        for (size_t c=0; c < table_labels[r].size(); c++) {
            table_labels[r][c] = new QLabel(this);
            table_labels[r][c]->setAlignment(Qt::AlignCenter);
            if (r > 0) {
                table_labels[r][c]->setFont(fixed_font); // so that the decimal point lines up down the columns
            }
            
            if (r == 0) {
                table_col_layouts.push_back(new QVBoxLayout());
                table_col_layouts.back()->setStretch(0, 0);
            }
            table_col_layouts[c]->addWidget(table_labels[r][c]);
            
            table_labels[r][c]->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
            table_labels[r][c]->setFixedWidth(std::max(default_fm.boundingRect(column_names[c].c_str()).width(), min_col_width));
            if (c >= 1) {
                table_labels[r][c]->setToolTip(column_tooltips[c].c_str());
            }
        }
    }
    for (size_t c=0; c < table_labels[0].size(); c++) {
        table_labels[0][c]->setStyleSheet("font-weight: bold;");
        table_labels[0][c]->setText(column_names[c].c_str());
        table_labels[0][c]->setAlignment(Qt::AlignCenter);
    }
    
    x_label = new QLabel(this);
    x_label->setText("Frequency");
    x_label->setStyleSheet("font-weight: bold;");
    x_label->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    x_label->setAlignment(Qt::AlignCenter);
    x_label->setFixedWidth(default_fm.boundingRect(QString(" ") + x_label->text() + QString(" ")).width());
    x_label_value = new QLabel(this);
    x_label_value->setText("");
    x_label_value->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    x_label_value->setFixedWidth(7*em_width);
    x_label_value->setAlignment(Qt::AlignCenter);
    
    QVBoxLayout* first_col = new QVBoxLayout;
    first_col->addWidget(x_label);
    first_col->addWidget(x_label_value);
    first_col->setStretch(0, 0);
    first_col->setSpacing(1);
    first_col->addStretch();
    
    table_row_layouts->addLayout(first_col);
    for (size_t c=0; c < table_labels[0].size(); c++) {
        table_row_layouts->addLayout(table_col_layouts[c]);
        table_col_layouts[c]->setStretch(0, 0);
        table_col_layouts[c]->setSpacing(1);
        table_col_layouts[c]->addStretch();
    }
    table_row_layouts->addStretch();
    
    save_data_button = new QPushButton("Save\ndata", this);
    save_data_button->setMinimumWidth(fm.horizontalAdvance(save_data_button->text()));
    save_data_button->setMaximumWidth(fm.horizontalAdvance(save_data_button->text()));
    save_data_button->setMinimumHeight(double_height);
    save_data_button->setMaximumHeight(double_height);
    save_data_button->setFlat(true);
    save_data_button->setStyleSheet(button_style);
    
    // generate an alpha-blended version of the logo
    QPixmap logo_pixmap(mtfmapper_logo->pixmap(50));
    QImage logo_image(logo_pixmap.size(), QImage::Format_ARGB32_Premultiplied);
    logo_image.fill(Qt::transparent);
    QPainter p(&logo_image);
    p.setOpacity(0.5);
    p.drawPixmap(0, 0, logo_pixmap);
    p.end();
    QPixmap blended_pixmap = QPixmap::fromImage(logo_image);
    
    QLabel* logo = new QLabel;
    logo->setAlignment(Qt::AlignRight);
    logo->setPixmap(blended_pixmap);
    logo->setIndent(blended_pixmap.size().width());
    logo->setMinimumSize(QSize(blended_pixmap.size().width()*2, blended_pixmap.size().height()));
    logo->setMaximumSize(QSize(blended_pixmap.size().width()*2, blended_pixmap.size().height()));
    
    QHBoxLayout* button_layout = new QHBoxLayout;
    button_layout->addWidget(save_img_button);
    button_layout->addWidget(save_data_button);
    
    QVBoxLayout* button_logo_layout = new QVBoxLayout;
    button_logo_layout->addWidget(logo);
    button_logo_layout->addLayout(button_layout);
    
    box_view = new QComboBox(this);
    box_view->addItem("SFR");
    box_view->addItem("ESF <ctrl>");
    box_view->addItem("LSF <alt>");
    box_view->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    
    QLabel* view_label = new QLabel("Plot type:", this);
    view_label->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    
    QString view_tooltip("Hold down <ctrl> to display the ESF,\n or <alt> to display the LSF");
    box_view->setToolTip(view_tooltip);
    view_label->setToolTip(view_tooltip);
    
    QLabel* lp_mm_label = new QLabel("lp/mm units", this);
    
    x_scale_label = new QLabel("distance scale:");
    x_scale_label->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    
    QSettings settings("mtfmapper", "mtfmapper");
    view.set_lp_mm_mode((Qt::CheckState)settings.value("setting_lpmm").toInt() == Qt::Checked);
    
    QHBoxLayout* view_layout = new QHBoxLayout;
    view_layout->addWidget(view_label);
    view_layout->addWidget(box_view);
    view_layout->addSpacing(50);
    view_layout->addWidget(lp_mm_cb);
    view_layout->addWidget(lp_mm_label);
    view_layout->addStretch();
    view_layout->addWidget(x_scale_label);
    view_layout->addWidget(x_scale_slider);
    view_layout->addStretch();
    
    auto vl_margins = view_layout->contentsMargins();
    view_layout->setContentsMargins(vl_margins.left(), vl_margins.top(), vl_margins.right(), vl_margins.bottom() + 10);
    
    QVBoxLayout* vmain_layout = new QVBoxLayout;
    vmain_layout->addLayout(view_layout);
    vmain_layout->addLayout(table_row_layouts);
    
    QHBoxLayout* main_layout = new QHBoxLayout;
    main_layout->addLayout(vmain_layout);
    main_layout->addLayout(button_logo_layout);
    
    QGroupBox* gbox = new QGroupBox("");
    gbox->setLayout(main_layout);
    
    QVBoxLayout* hlayout = new QVBoxLayout;
    hlayout->addWidget(chart_view);
    hlayout->addWidget(gbox);
    hlayout->setStretchFactor(chart_view, 1.0);
    hlayout->setStretchFactor(gbox, 0.0);
    setLayout(hlayout);
    
    for (int i=0; i < 3; i++) {
        mtf50_text.push_back(new QGraphicsSimpleTextItem(chart));
    }
    
    mtf50_rect = new QGraphicsRectItem(chart);
    mtf50_rect->setRect(0,0,0,0);
    QColor charcoal(54, 69, 79, 100);
    mtf50_rect->setBrush(QBrush(charcoal));
    mtf50_rect->setPen(charcoal);
    
    setAutoFillBackground(true); 
    QPalette palette;
    palette.setColor(backgroundRole(), Qt::white);
    setPalette(palette);
    
    hide_xscale(); // because we start in SFR mode
    x_scale_changed(x_scale_slider->value());
    
    chart->resize(650, 350);
    setMinimumHeight(400);
    setMinimumWidth(750);
    setWindowTitle(view.title().c_str());
    
    connect(save_img_button, SIGNAL(clicked()), this, SLOT(save_image()));
    connect(save_data_button, SIGNAL(clicked()), this, SLOT(save_data()));
    connect(box_view, SIGNAL(currentIndexChanged(int)), this, SLOT(plot_type_changed(int)));
    connect(lp_mm_cb, SIGNAL(clicked()), this, SLOT(lp_mm_toggled()));
    connect(x_scale_slider, SIGNAL(valueChanged(int)), this, SLOT(x_scale_changed(int)));
}

void Sfr_dialog::clear(void) {
    chart->removeAllSeries();
    series.clear();
    for (size_t i=1; i < table_labels.size(); i++) {
        for (auto ll: table_labels[i]) {
            ll->setText("");
            set_label_background(ll, "none");
        }
    }
    for (auto m: mtf50_text) {
        m->setText("");
    }
    // TODO: Not sure why we have to forcibly delete all the entries; maybe we can just flag them
    for (auto e: entries) {
        e.clear();
    }
    entries.clear();
    y_axis->setMax(0.1);
}

void Sfr_dialog::reject(void) {
    clear();
    emit sfr_dialog_closed();
    emit send_geometry(geometry());
    QDialog::reject();
}

void Sfr_dialog::paintEvent(QPaintEvent* event) {
    repainting.testAndSetAcquire(0, 1);
    QDialog::paintEvent(event);
    
    bool changed = false;
    
    Qt::KeyboardModifiers current_modifier = QGuiApplication::queryKeyboardModifiers() & ~Qt::ShiftModifier;
    if (current_modifier != last_modifier) {
        // update
        if (current_modifier & Qt::ControlModifier && !(current_modifier & Qt::AltModifier)) {
            view.push_view();
            view.set_view(Entry_view::VIEW_ESF);
        } else {
            if (current_modifier & Qt::AltModifier && !(current_modifier & Qt::ControlModifier)) {
                view.push_view();
                view.set_view(Entry_view::VIEW_LSF);
            } else {
                view.pop_view();
            }
        }
        box_view->setCurrentIndex((int)view.get_view());
        changed = true;
    }
    
    if (changed) {
        view.update(entries, series, *chart, *x_axis, *y_axis, xscale());
    }
    
    last_modifier = current_modifier;

    vector<double> contrast_list;
    if (series.size() > 0) {

        for (size_t si=0; si < series.size(); si++) {
        
            QVector<QPointF> pts = series[si]->pointsVector();
            double contrast = pts[0].y();
            for (int i = 0; i < pts.size() && pts[i].x() < cursor_domain_value; i++) {
                contrast = pts[i].y();
            }
            
            contrast_list.push_back(contrast);
        }
        

        QFontMetrics fm(QWidget::fontMetrics());
        int th = fm.height();
        int text_end = chart->size().width() - 2*fm.horizontalAdvance("m");
        
        // add mtf50 tags in reverse
        for (int mi=series.size()-1; mi >= 0; mi--) {
            
            string mtf50_str = view.mtf_xx(entries[mi].info);
            int tw = fm.horizontalAdvance(mtf50_str.c_str()) + 2*fm.horizontalAdvance("m");
            
            mtf50_text[mi]->setBrush(series[mi]->pen().color());
            mtf50_text[mi]->setPen(QPen(series[mi]->pen().color()));
            mtf50_text[mi]->setPos(text_end - tw, th);
            mtf50_text[mi]->setText(mtf50_str.c_str());
            text_end -= tw;
        }
    }
    
    QPointF tpos(cursor_domain_value, y_axis->max() - 0.001);
    QPointF bpos(cursor_domain_value, y_axis->min() + 0.001);
    tpos = chart->mapToPosition(tpos);
    bpos = chart->mapToPosition(bpos);
    tpos.rx() = std::max(tpos.x()-1, chart->mapToPosition(QPointF(x_axis->min(), 0)).x());
    bpos.rx() = tpos.rx() + 2;
    
    mtf50_rect->setRect(tpos.x(), tpos.y(), bpos.x() - tpos.x(), bpos.y() - tpos.y());

    
    for (size_t mi=0; mi < contrast_list.size(); mi++) {
        auto fmt = view.format(entries[mi].info, cursor_domain_value, contrast_list[mi]);
        table_labels[mi+1][0]->setText(fmt["y_value"].c_str());
        table_labels[mi+1][1]->setText(fmt["cnr"].c_str());
        table_labels[mi+1][2]->setText(fmt["r_ca"].c_str());
        table_labels[mi+1][3]->setText(fmt["b_ca"].c_str());
        table_labels[mi+1][4]->setText(fmt["angle"].c_str());
        table_labels[mi+1][5]->setText(fmt["ox"].c_str());
        table_labels[mi+1][6]->setText(fmt["length"].c_str());
        table_labels[mi+1][7]->setText(fmt["channel"].c_str());
        
        if (mi == 0) {
            table_labels[0][0]->setText(fmt["y_label"].c_str());
            x_label->setText(fmt["x_label"].c_str());
            x_label_value->setText(fmt["x_value"].c_str());
        }
        
        for (size_t col=0; col < table_labels[0].size(); col++) {
            QPalette palette = table_labels[mi+1][col]->palette();
            palette.setColor(table_labels[mi+1][col]->foregroundRole(), series[mi]->pen().color());
            palette.setColor(table_labels[mi+1][col]->backgroundRole(), this->palette().color(QPalette::Window));
            table_labels[mi+1][col]->setPalette(palette);
            table_labels[mi+1][col]->setAutoFillBackground(false);
            table_labels[mi+1][col]->show();
        }

        set_label_background(table_labels[mi+1][1], fmt["cnr_col"]);
        set_label_background(table_labels[mi+1][4], fmt["angle_col"]);
        set_label_background(table_labels[mi+1][5], fmt["ox_col"]);
        set_label_background(table_labels[mi+1][6], fmt["length_col"]);
    }
    
    chart->update(); // this causes some lag, but it elliminates exposed cruft. a better solution would be nice
    repainting.testAndSetRelease(1, 0);
}

void Sfr_dialog::replace_entry(const Sfr_entry& entry) {
    if (entries.size() < 3) {
        entries.push_back(entry);
    } else {
        entries.back() = entry;
    }
    view.update(entries, series, *chart, *x_axis, *y_axis, xscale());
    
    update();
}

void Sfr_dialog::add_entry(const Sfr_entry& entry) {
    replace_entry(entry);
}

void Sfr_dialog::notify_mouse_position(double value, bool click) { // we still need this to update the bottom label, but we could merge this in the chart?
    lock_cursor ^= click;
    if (!lock_cursor) {
        cursor_domain_value = std::min(x_axis->max(), std::max(x_axis->min(), value));
        update();
    }
}

void Sfr_dialog::save_image(void) {
    QString savename = QFileDialog::getSaveFileName(
        this,
        tr("Save plot image"),
        QString(),
        QString("*.png")
    );
    
    // must append extension if none were specified
    if (!savename.contains('.')) {
        savename += ".png";
    }
    
    QImage img(size(), QImage::Format_RGB888);
    QPainter painter(&img);
    render(&painter);
    img.save(savename);
}

void Sfr_dialog::save_data(void) {
    setlocale(LC_ALL, "C");

    QString savename = QFileDialog::getSaveFileName(
        this,
        tr("Save CSV data"),
        QString(),
        QString("*.csv")
    );
    
    // must append extension if none were specified
    if (!savename.contains('.')) {
        savename += ".csv";
    }
    
    FILE* fout = fopen(savename.toLocal8Bit().data(), "wt");
    if (fout) {
        fprintf(fout, "frequency,");
        if (series.size() == 1) {
            fprintf(fout, "contrast\n");
            for (int i=0; i < series[0]->count()/20; i++) {
                fprintf(fout, "%.4lf,%.8lf\n", series[0]->at(i*20).x(), series[0]->at(i*20).y());
            }
        } else {
            size_t max_size = 0;
            for (size_t i = 0; i < series.size(); i++) {
                max_size = std::max(max_size, size_t(series[i]->count()));
            }
            max_size /= 20;
            int j;
            for (j=0; j < (int)series.size() - 1; j++) {
                fprintf(fout, "contrast_%d,", j);
            }
            fprintf(fout, "contrast_%d\n", j);
            
            for (size_t i=0; i < max_size; i++) {
                fprintf(fout, "%.4lf,", (int)i*20 < series[0]->count() ? series[0]->at(i*20).x() : 0.0);
                for (j=0; j < (int)series.size() - 1; j++) {
                    fprintf(fout, "%.8lf,", (int)i*20 < series[j]->count() ? series[j]->at(i*20).y() : 0.0);
                }
                fprintf(fout, "%.8lf\n", (int)i*20 < series[j]->count() ? series[j]->at(i*20).y() : 0.0);
            }
        }
        fclose(fout);
    } else {
        // display failure dialog
    }
}

bool Sfr_dialog::update_lp_mm_mode(void) {
    QSettings settings("mtfmapper", "mtfmapper");
    
    view.set_lp_mm_mode((Qt::CheckState)settings.value("setting_lpmm").toInt() == Qt::Checked);
    view.set_default_pixel_pitch(settings.value("setting_pixelsize").toFloat());
    lp_mm_cb->setCheckState(view.get_lp_mm_mode() ? Qt::Checked : Qt::Unchecked);
    
    view.update(entries, series, *chart, *x_axis, *y_axis, xscale());
    update();
    
    return true;
}

void Sfr_dialog::lp_mm_toggled(void) {
    view.set_lp_mm_mode(lp_mm_cb->checkState() == Qt::Checked);
    view.update(entries, series, *chart, *x_axis, *y_axis, xscale());
    update();
}

void Sfr_dialog::plot_type_changed(int index) {
    
    switch(index) {
    case 0: view.set_view(Entry_view::VIEW_SFR); 
        hide_xscale();
        break;
    case 1: view.set_view(Entry_view::VIEW_ESF); 
        show_xscale();
        break;
    case 2: view.set_view(Entry_view::VIEW_LSF); 
        show_xscale();
        break;
    }
    
    view.update(entries, series, *chart, *x_axis, *y_axis, xscale());
    setWindowTitle(view.title().c_str());
}

void Sfr_dialog::keyPressEvent(QKeyEvent* /*event*/) {
    update();
}

void Sfr_dialog::keyReleaseEvent(QKeyEvent* /*event*/) {
    update();
}

void Sfr_dialog::set_label_background(QLabel* label, const string& condition) {
    QColor bg_red(255, 128, 128, 32);
    QColor bg_yellow(255, 255, 128, 48);

    if (condition != "none") {
        label->setAutoFillBackground(true);
        QPalette pal = label->palette();
        pal.setColor(
            label->backgroundRole(),
            condition == "red" ? bg_red : bg_yellow
        );
        label->setPalette(pal);
    } else {
        label->setAutoFillBackground(false);
    }
}

void Sfr_dialog::x_scale_changed([[maybe_unused]] int val) {
    view.update(entries, series, *chart, *x_axis, *y_axis, xscale());
    update();
}

void Sfr_dialog::show_xscale(void) {
    x_scale_label->show();
    x_scale_slider->show();
}

void Sfr_dialog::hide_xscale(void) {
    x_scale_label->hide();
    x_scale_slider->hide();
}

int Sfr_dialog::xscale(void) { 
    if (x_scale_slider == nullptr) {
        return std::max(0, (int)x_scale_levels.size() - 2);
    }
    int val = std::max(0, std::min(x_scale_slider->value(), (int)x_scale_levels.size() - 1));
    return x_scale_levels[val];
}

