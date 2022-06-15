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
#ifndef BLOCK_H
#define BLOCK_H

#include "include/rectangle.h"
#include "include/snr.h"
#include "include/edge_model.h"
#include <memory>

#include <map>
using std::map;

#include "include/edge_info.h"

class Block {
  public:
    typedef enum {TOP, LEFT, RIGHT, BOTTOM} edge_position;    

    Block(void) : rect(Mrectangle()), mtf50(4,0.0), quality(4, 0.0), 
        sfr(4), esf(4), centroid(0,0), area(0.0), valid(true), 
        line_deviation(4, cv::Point3d(0,0,1)), snr(4), scansets(4), 
        chromatic_aberration(4, cv::Point2d(Edge_info::nodata, Edge_info::nodata)), 
        edge_model(4), valid_edge(4, false), edge_length(4, 0.0) {
    }

    Block(const Mrectangle& in_rect) : rect(in_rect), mtf50(4,0.0), 
        quality(4, 0.0), sfr(4), esf(4), centroid(0,0), area(0.0), valid(true), 
        line_deviation(in_rect.line_deviation), snr(4), scansets(4), 
        chromatic_aberration(4, cv::Point2d(Edge_info::nodata, Edge_info::nodata)),
        edge_model(4), valid_edge(4, false), edge_length(4, 0.0) {
    
        size_t top_edge_idx = 0;
        size_t bot_edge_idx = 0;
        size_t left_edge_idx = 0;
        size_t right_edge_idx = 0;
        
        if (rect.centroids.size() != 4) {
            valid = false;
            return;
        }
     
        for (size_t i=0; i < 4; i++) {
            centroid.x += get_edge_centroid(i).x;
            centroid.y += get_edge_centroid(i).y;
        
            if (get_edge_centroid(i).y < get_edge_centroid(top_edge_idx).y) {
                top_edge_idx = i;
            }
            
            if (get_edge_centroid(i).y > get_edge_centroid(bot_edge_idx).y) {
                bot_edge_idx = i;
            }
            
            if (get_edge_centroid(i).x < get_edge_centroid(left_edge_idx).x) {
                left_edge_idx = i;
            }
            
            if (get_edge_centroid(i).x > get_edge_centroid(right_edge_idx).x) {
                right_edge_idx = i;
            }
        }
        centroid.x /= 4;
        centroid.y /= 4;
        
        edge_lut[TOP]    = top_edge_idx;
        edge_lut[BOTTOM] = bot_edge_idx;
        edge_lut[LEFT]   = left_edge_idx;
        edge_lut[RIGHT]  = right_edge_idx;
        
        Point2d e1(
            get_edge_centroid(edge_lut[BOTTOM]).x - get_edge_centroid(edge_lut[TOP]).x,
            get_edge_centroid(edge_lut[BOTTOM]).y - get_edge_centroid(edge_lut[TOP]).y
        );
        
        Point2d e2(
            get_edge_centroid(edge_lut[RIGHT]).x - get_edge_centroid(edge_lut[LEFT]).x,
            get_edge_centroid(edge_lut[RIGHT]).y - get_edge_centroid(edge_lut[LEFT]).y
        );
        
        area = sqrt(SQR(e1.x) + SQR(e1.y)) * sqrt(SQR(e2.x) + SQR(e2.y));
    }
    
    void set_snr(size_t edge_number, const Snr& in_snr) {
        snr[edge_number] = in_snr;
    }
    
    const Snr& get_snr(size_t edge_number) const {
        return snr[edge_number];
    }

    void set_sfr(size_t edge_number, const vector<double>& in_sfr) {
        sfr[edge_number] = std::shared_ptr<vector<double>>(new vector<double>(in_sfr));
    }
    
    const vector<double>& get_sfr(size_t edge_number) const {
        return *sfr[edge_number];
    }
    
    void set_esf(size_t edge_number, const vector<double>& in_esf) {
        esf[edge_number] = std::shared_ptr<vector<double>>(new vector<double>(in_esf));
    }

    const vector<double>& get_esf(size_t edge_number) const {
        return *esf[edge_number];
    }
    
    const vector<Point2d>& get_ridge(size_t edge_number) const {
        return edge_model[edge_number]->ridge;
    }

    void set_normal(size_t edge_number, const Point2d& rgrad) {
        rect.normals[edge_number] = rgrad;
    }

    Point2d get_normal(size_t edge_number) const {
        return rect.normals[edge_number];
    }

    Point2d get_edge(size_t edge_number) const {
        assert(edge_number < 4);
        return rect.edges[edge_number];
    }
    
    Point2d get_corner(size_t edge_number) const {
        assert(edge_number < 4);
        return rect.corners[edge_number];
    }
    
    double get_edge_angle(size_t edge_number) const {
        assert(edge_number < 4);
        return rect.thetas[edge_number];
    }
    
    Point2d get_edge_centroid(size_t edge_number) const {
        assert(edge_number < 4);
        return rect.centroids[edge_number];
    }
    
    Point2d get_edge_centroid(edge_position ep) const {
        map<edge_position, size_t>::const_iterator it = edge_lut.find(ep);
        return rect.centroids[it->second];
    }
    
    // quality of 1.0 means good quality, 0.0 means unusably poor
    void set_mtf50_value(size_t edge_number, double mtf50_value, double quality_value) {
        assert(edge_number < 4);
        mtf50[edge_number] = mtf50_value;
        quality[edge_number] = quality_value;
    }
    
    double get_mtf50_value(size_t edge_number) const {
        assert(edge_number < 4);
        return mtf50[edge_number];
    }
    
    double get_mtf50_value(edge_position ep) const {
        map<edge_position, size_t>::const_iterator it = edge_lut.find(ep);
        return get_mtf50_value(it->second);
    }
    
    Point2d get_centroid(void) const {
        return centroid;
    }
    
    double get_area(void) const {
        return area;
    }
    
    double get_quality(size_t edge_number) const {
        return quality[edge_number];
    }
    
    double get_quality(edge_position ep) const {
        map<edge_position, size_t>::const_iterator it = edge_lut.find(ep);
        return quality[it->second];
    }
    
    int get_edge_index(edge_position ep) const {
        map<edge_position, size_t>::const_iterator it = edge_lut.find(ep);
        return it->second;
    }
    
    cv::Point3d get_line_deviation(size_t edge_number) const {
        assert(edge_number < 4);
        return line_deviation[edge_number];
    }
    
    void set_line_deviation(size_t edge_number, const cv::Point3d& deviation) {
        assert(edge_number < 4);
        line_deviation[edge_number] = deviation;
    }
    
    void set_scanset(size_t edge_number, const map<int, scanline>& scanset) {
        assert(edge_number < 4);
        scansets[edge_number] = std::shared_ptr<map<int, scanline>>(new map<int, scanline>(scanset));
    }
    
    const map<int, scanline>& get_scanset(size_t edge_number) const {
        assert(edge_number < 4);
        return *scansets[edge_number];
    }
    
    void set_ca(size_t edge_number, Point2d ca) {
        assert(edge_number < 4);
        chromatic_aberration[edge_number] = ca;
    }
    
    const Point2d& get_ca(size_t edge_number) const {
        assert(edge_number < 4);
        return chromatic_aberration[edge_number];
    }
    
    void set_edge_model(size_t edge_number,  const std::shared_ptr<Edge_model>& em) {
        assert(edge_number < 4);
        edge_model[edge_number] = em;
    }
    
    Edge_model& get_edge_model(size_t edge_number) {
        assert(edge_number < 4);
        return *edge_model[edge_number];
    }
    
    bool get_edge_valid(size_t edge_number) const {
        assert(edge_number < 4);
        return valid_edge[edge_number];
    }
    
    void set_edge_valid(size_t edge_number) {
        assert(edge_number < 4);
        valid_edge[edge_number] = true;
    }

    double get_edge_length(size_t edge_number) const {
        assert(edge_number < 4);
        return edge_length[edge_number];
    }

    void set_edge_length(size_t edge_number, double length) {
        assert(edge_number < 4);
        edge_length[edge_number] = length;
    }

    bool serialize(FILE* fout) const {
        vector<Point2d> edge_centroids(4);
        vector<double> edge_angle(4);
        vector<Point2d> edge_snr(4);
        vector<double> geometric_length(4);
        
        for (size_t k=0; k < 4; k++) {
            edge_centroids[k] = rect.get_centroid(k);
            edge_angle[k] = atan2(-get_normal(k).x, get_normal(k).y);
            edge_snr[k] = cv::Point2d(snr[k].mean_cnr(), snr[k].oversampling());
            vector<std::pair<double, size_t>> corner_dist(4);
            for (size_t j=0; j < 4; j++) {
                corner_dist[j] = std::pair<double, size_t>(cv::norm(edge_centroids[k] - get_corner(j)), j);
            }
            sort(corner_dist.begin(), corner_dist.end());
            geometric_length[k] = cv::norm(get_corner(corner_dist[0].second) - get_corner(corner_dist[1].second));
        }
        
        return Edge_info::serialize(
            fout,
            edge_centroids,
            edge_angle,
            mtf50,
            quality,
            sfr,
            esf,
            edge_snr,
            chromatic_aberration,
            valid_edge,
            edge_length,
            geometric_length
        );
    }
    
    Mrectangle rect;
    vector<double> mtf50;
    vector<double> quality;
    vector<std::shared_ptr<vector<double>>> sfr;
    vector<std::shared_ptr<vector<double>>> esf;
    map<edge_position, size_t> edge_lut;
    Point2d centroid;
    double area;
    bool valid;
    vector<cv::Point3d> line_deviation;
    vector<Snr> snr;
    vector<std::shared_ptr<map<int, scanline>>> scansets; 
    vector<Point2d> chromatic_aberration;
    vector<std::shared_ptr<Edge_model>> edge_model;
    vector<bool> valid_edge;
    vector<double> edge_length;
};

#endif
