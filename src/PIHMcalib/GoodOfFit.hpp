//
//  GoodOfFit.hpp
//  PIHMcalib
//
//  Created by Lele Shu on 8/18/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#ifndef GoodOfFit_hpp
#define GoodOfFit_hpp

#include <stdio.h>
#include <math.h>
#include "quicksort.hpp"
#include "FDC.hpp"
#include "functions.hpp"

void FDC(double *x, double *px, double *py, int n); /* Flow Duration Curve */

class GoodOfFit {
private:
    int initflag = 0;
    double max(double *x, int n); /**/
    double min(double *x, int n); /**/
    double mean(double *x, int n); /**/
    double sum(double *x, int n); /**/
    double sum_sq_diff(double *x, double *y, int n);
    double sum_sq_diff(double *x, double y, int n);
    void product(double *x, double *y, int n, double *ret);
    double sd(double *x, int n); /* Standard Deviation */
    double which_q(double *sortx, double *exceedance, double q);
    double *temp;
public:
    int nLen;
    double *sim;
    double *obs;
    double *ssim; // sorted
    double *sobs; // sorted
    double *psim; // Exceedence
    double *pobs; // Exceedence
    
    GoodOfFit();
    ~GoodOfFit();
    void init(double *xsim, double *xobs, int n);
    
    double gof_NSE(); /*Nash Efficiency*/
    double gof_RMSE(); /*Root Mean Sq Error*/
    double gof_nRMSE(int imean); /* Normalized Root Mean Sq Error*/
    double gof_MSE(); /* Mean Sq Error*/
    double gof_pBias(); /* Percent Bias*/
    double gof_pBiasFDC(); /* Percent Bias of the midsegment of the Flow Duration Curve */
    double gof_pBiasFDC(double lq, double hq);
    double gof_R2(); /* */
    double gof_ME(); /* Mean Error*/
    double gof_MAE(); /* Mean Absolute Error */
    double gof_ssq(); /* Sum of Squared Residuals between 'sim' and 'obs' */
    void callFDC();
    void print_gofname(FILE *fp);
    void print_gof(FILE *fp);
};

#endif /* GoodOfFit_hpp */

