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
#ifndef FOCUS_SURFACE_H
#define FOCUS_SURFACE_H

#include "include/logger.h"
#include "ratpoly_fit.h"
#include "mtf_profile_sample.h"
#include "distance_scale.h"
#include "mtf50_edge_quality_rating.h"

class Focus_surface  {
  public:
  
    Focus_surface(vector< Mtf_profile_sample >& data, int order_n, int order_m, Distance_scale& distance_scale) 
    : data(data), order_n(order_n), order_m(order_m), maxy(-1e50), maxx(-1e50) {
    
        double miny = 1e50;
        for (size_t i=0; i < data.size(); i++) {
            maxx = max(fabs(data[i].p.x), maxx);
            maxy = max(fabs(data[i].p.y), maxy);
            miny = min(fabs(data[i].p.y), miny);
        }
        
        logger.debug("chart extents: y:(%lf, %lf) mm, x:(%lf) mm\n", miny, maxy, maxx);
        
        // TODO: this loop is a great candidate for OpenMP, but some care must be taken
        // to protect the (few) shared structures (like "peaks")
        
        double min_fit_err = 1e50;
        VectorXd best_sol;
        vector<Sample> dummy_data;
        Ratpoly_fit best_fit(dummy_data, 5, 5);
        
        FILE* ffout = fopen("peaks.txt", "wt");
        
        maxy = 136; // empirical?
        
        vector<Sample> peak_pts;
        // |15 to 110| in steps of 2.5, width=5 ??
        for (int s=-1; s <= 1; s+=2) {
            //for (double d=max(10.0, 1.2*miny*cscale); d <= 0.95*maxy*cscale; d += 1.0) {
            for (double d=10.0; d <= 0.9*maxy; d += 2.0) {
                double midy = s*d;
                
                double mean_x = 0;
                double wsum = 0;
                
                vector<Sample> pts_row;
                for (size_t i=0; i < data.size(); i++) { // TODO: this may become a little slow with many points
                    double dy = midy - data[i].p.y;
                    if (fabs(dy) < 15 && fabs(data[i].p.y) > 5 && fabs(data[i].p.x) < 175 && fabs(data[i].p.y) < 136) { // at least 5 mm from centre of chart
                        
                        double yw = exp(-dy*dy/(2*3*3)); // sdev of 3 mm in y direction
                        pts_row.push_back( Sample(data[i].p.x, data[i].mtf, yw, data[i].quality <= poor_quality ? 0.25 : 1) );
                        mean_x += pts_row.back().weight * data[i].p.y;
                        wsum += pts_row.back().weight;
                    } 
                }
                
                if (pts_row.size() < 3*14) {
                    continue; 
                }
                
                // now filter out the really bad outliers
                const int sh = 4;
                const double sgw[] = {-21/231.0, 14/231.0, 39/231.0, 54/231.0, 59/231.0, 54/231.0, 39/231.0, 14/231.0, -21/231.0};
                sort(pts_row.begin(), pts_row.end());
                
                // just pretend our samples are equally spaced
                vector<Sample> ndata;
                for (size_t i=sh; i < pts_row.size() - sh; i++) {
                    double val = 0;
                    for (int w=-sh; w <= sh; w++) {
                        val += sgw[w+sh] * pts_row[i+w].y;    
                    }
                    
                    if (fabs(val - pts_row[i].y)/val < 0.05) {
                        ndata.push_back(pts_row[i]);
                    }
                    
                }
                pts_row = ndata;
                
                VectorXd sol;
                double lpeak;
                
                mean_x /= wsum;
                
                Ratpoly_fit cf(pts_row, 4, order_m);
                sol = rpfit(cf, true, true); 
                while (cf.order_n > 1 && cf.order_m > 0 && cf.has_poles(sol)) {
                    cf.order_m--;
                    sol = rpfit(cf, true, true);
                    if (cf.has_poles(sol)) {
                        cf.order_n--;
                        sol = rpfit(cf, true, true);
                    }
                }
                if (cf.has_poles(sol)) { 
                    // no solution without poles, give up, skip this sample?
                    logger.debug("Warning: no viable RP fit. Skipping curve centred at y=%lf\n", mean_x);
                    continue;
                }
                
                double err = cf.evaluate(sol);
                lpeak = cf.peak(sol); 
                
                double merr = err/double(pts_row.size());
                if (merr < min_fit_err) {
                    logger.debug("min fit err %lf at dist %lf\n", merr, midy);
                    min_fit_err = merr;
                }
                if (fabs(midy) < 20 && fabs(merr - min_fit_err)/merr < 1) {
                    best_sol = sol;
                    dummy_data = pts_row;
                    best_fit.order_n = cf.order_n;
                    best_fit.order_m = cf.order_m;
                    best_fit.xs_min = cf.xs_min;
                    best_fit.xs_scale = cf.xs_scale;
                    best_fit.ysf = cf.ysf;
                }
                
                peak_pts.push_back( Sample(mean_x, lpeak, 1, 1.0) );
                ridge_peaks.push_back(Point2d(lpeak, mean_x));
                
                fprintf(ffout, "%lf %lf\n", mean_x, lpeak);
            }
        }
        
        fclose(ffout);
        
        if (peak_pts.size() < 10) {
            logger.error("%s\n", "Not enough peak points to construct peak focus curve.");
            return;
        }
        
        Ratpoly_fit cf(peak_pts, 2,2);
        cf.base_value = 1;
        cf.pscale = 0;
        VectorXd sol = rpfit(cf, true, true);
        while (cf.order_m > 0 && cf.has_poles(sol)) {
            cf.order_m--;
            logger.debug("reducing order_m to %d\n", cf.order_m);
            sol = rpfit(cf, true, true);
        }
        if (cf.has_poles(sol)) { 
            // no solution without poles, give up, skip this sample?
            logger.debug("%s\n", "Warning: no viable RP fit to fpeaks data");
        }
        
        #if 1
        // attempt IRLS step
        vector<double> iweights(peak_pts.size(), 0);
        for (size_t i=0; i < peak_pts.size(); i++) {
            iweights[i] = peak_pts[i].yweight;
        }
        double prev_err = 1e50;
        for (int iter=0; iter < 50; iter++) {
            double errsum = 0;
            double wsum = 0;
            for (size_t i=0; i < peak_pts.size(); i++) {
                double y = cf.rpeval(sol, cf.scale(peak_pts[i].x)) / cf.ysf;
                double e = fabs(y - peak_pts[i].y);
                peak_pts[i].yweight = iweights[i] / max(0.0001, e);
                
                if (iter > 3 && e > 20) {
                    logger.debug("outlier at %lf mm, suppressing\n", e);
                    peak_pts[i].yweight = 0;
                }
                
                double w = peak_pts[i].yweight;
                errsum += e*w;
                wsum += w;
                
            }
            errsum /= wsum;
            logger.debug("iter %d err: %lf\n", iter, errsum);
            VectorXd oldsol = sol;
            sol = rpfit(cf, true, true);
            if (iter > 10 && (prev_err - errsum)/prev_err < 0.0001) {
                logger.debug("bailing out at iter %d\n", iter);
                if (errsum > prev_err) {
                    logger.debug("%s\n", "reverting to older solution");
                    sol = oldsol;
                }
                break;
            }
            prev_err = errsum;
        }
        #endif
        
        // now perform some bootstrapping to obtain bounds on the peak focus curve:
        vector<double> mc_pf;
        map<double, vector<double> > mc_curve;
        for (int iters=0; iters < 30; iters++) {
            vector<Sample> sampled_peak_pts;
            for (int j=0; j < peak_pts.size()*0.5; j++) {
                int idx = (int)floor(peak_pts.size()*double(rand())/double(RAND_MAX));
                sampled_peak_pts.push_back(peak_pts[idx]);
            }
            Ratpoly_fit mc_cf(sampled_peak_pts, cf.order_n, cf.order_m);
            mc_cf.base_value = 1;
            mc_cf.pscale = 0;
            VectorXd mc_sol = rpfit(mc_cf, true, true);
            mc_pf.push_back(mc_cf.rpeval(mc_sol, 0)/mc_cf.ysf);
            
            for (double y=-maxy; y < maxy; y += 10) {
                double x = mc_cf.rpeval(mc_sol, mc_cf.scale(y))/mc_cf.ysf;
                mc_curve[y].push_back(x);
            }
            
        }
        sort(mc_pf.begin(), mc_pf.end());
        for (map<double, vector<double> >::iterator it = mc_curve.begin(); it != mc_curve.end(); it++) {
            sort(it->second.begin(), it->second.end());
            ridge_p05.push_back(Point2d(it->second[0.05*it->second.size()], it->first));
            ridge_p95.push_back(Point2d(it->second[0.95*it->second.size()], it->first));
        }
        
        for (double y=-maxy; y < maxy; y += 1) {
            double x = cf.rpeval(sol, cf.scale(y))/cf.ysf;
            ridge.push_back(Point2d(x, y));
        }
        
        double x_inter = cf.rpeval(sol, 0)/cf.ysf;
        
        int x_inter_index = lower_bound(mc_pf.begin(), mc_pf.end(), x_inter) - mc_pf.begin();
        logger.debug("x_inter percentile: %.3lf\n", x_inter_index * 100 / double(mc_pf.size()));
        logger.debug("x_inter 95%% confidence interval: [%lf, %lf]\n", mc_pf[0.05*mc_pf.size()], mc_pf[0.95*mc_pf.size()]);
        
        distance_scale.estimate_depth_world_coords(x_inter, 0.0, focus_peak);
        distance_scale.estimate_depth_world_coords(mc_pf[0.05*mc_pf.size()], 0.0, focus_peak_p05);
        distance_scale.estimate_depth_world_coords(mc_pf[0.95*mc_pf.size()], 0.0, focus_peak_p95);
        
        logger.debug("focus_mm (on chart x axis) %lg\n", x_inter);
        
        logger.debug("focus_plane %lg\n", focus_peak);
        logger.debug("fp_interval: [%lf, %lf]\n", focus_peak_p05, focus_peak_p95);
        
        double curve_min = 1e50;
        double curve_max = -1e50;
        for (size_t i=0; i < dummy_data.size(); i++) {
            curve_min = min(curve_min, dummy_data[i].x);
            curve_max = max(curve_max, dummy_data[i].x);
        }
        FILE* profile = fopen("nprofile.txt", "wt");
        double curve_peak = best_fit.peak(best_sol);
        double curve_offset = curve_peak - x_inter;
        logger.debug("curve peak = %lf, x_inter = %lf\n", curve_peak, x_inter);
        for (double cx=curve_min; cx <= curve_max; cx += 1) {
            double mtf = best_fit.rpeval(best_sol, best_fit.scale(cx))/best_fit.ysf;
            double depth = 0;
            distance_scale.estimate_depth_world_coords(cx - curve_offset, 0.0, depth);
            fprintf(profile, "%lf %lf\n", depth, mtf);
        }
        fclose(profile);
    }
    
    VectorXd rpfit(Ratpoly_fit& cf, bool scale=true, bool refine=true) {
        const vector<Sample>& pts_row = cf.get_data();
        
        if (scale) {
            double xmin = 1e50;
            double xmax = -1e50;
            double ysf=0;
            for (size_t i=0; i < pts_row.size(); i++) {
                xmin = std::min(xmin, pts_row[i].x);
                xmax = std::max(xmax, pts_row[i].x);
                ysf = max(ysf, fabs(pts_row[i].y));
            }
            cf.xs_min = 0.5*(xmin + xmax);
            cf.xs_scale = 2.0/(xmax - xmin);
            cf.ysf = ysf = 1.0/ysf;
        }
        
        int tdim = cf.dimension();
        MatrixXd cov = MatrixXd::Zero(tdim, tdim);
        VectorXd b = VectorXd::Zero(tdim);
        VectorXd a = VectorXd::Zero(tdim);
        
        VectorXd sol;
        
        for (int iter=0; iter < 1; iter++) {
            cov.setZero();
            b.setZero();
            a.setZero();
            
            for (size_t i=0; i < pts_row.size(); i++) {
                const Sample& sp = pts_row[i];
                double w = sp.weight * sp.yweight;
                a[0] = 1*w;
                double prod = cf.scale(sp.x); // top poly
                for (int j=1; j <= cf.order_n; j++) {
                    a[j] = w*cf.cheb(j, prod);
                }
                // bottom poly
                for (int j=1; j <= cf.order_m; j++) {
                    a[j+cf.order_n] = cf.cheb(j, prod)*w*sp.y*cf.ysf;
                }
                
                for (int col=0; col < tdim; col++) { 
                    for (int icol=0; icol < tdim; icol++) { // build covariance of design matrix : A'*A
                        cov(col, icol) += a[col]*a[icol];
                    }
                    b[col] += cf.base_value*a[col]*sp.y*cf.ysf*w; // build rhs of system : A'*b
                }
            }
            
            for (int col=cf.order_n+1; col < cov.cols(); col++) {
                cov.col(col) = -cov.col(col);
            }
            
            sol = cov.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
        }
        
        // now perform non-linear optimization
        
        if (refine) {
            sol = cf.gauss_newton_armijo(sol);
        }
        
        return sol;
    }

    double evaluate(VectorXd&) {
        return 0;
    }
    
    int dimension(void) {
        return (order_n+1 + order_m);
    }
    
    double softclamp(double x, double lower, double upper, double p=0.98) {
        double s = (x - lower) / (upper - lower);
        if (s > p) {
            return 1.0/(1.0 + exp(-3.89182*s));
        }
        return s < 0 ? 0 : s;
    }
    
    vector< Mtf_profile_sample >& data;
    int order_n;
    int order_m;
    double maxy;
    double maxx;
    vector<Point2d> ridge;
    vector<Point2d> ridge_peaks;
    vector<Point2d> ridge_p05;
    vector<Point2d> ridge_p95;
    double focus_peak;
    double focus_peak_p05;
    double focus_peak_p95;
    
    double p2;
    double p98;
};

#endif
