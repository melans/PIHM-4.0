
#include "GoodOfFit.hpp"

GoodOfFit::GoodOfFit(){
}
GoodOfFit::~GoodOfFit(){
    if(initflag){
        delete[] sim;
        delete[] obs;
        delete[] ssim;
        delete[] sobs;
        delete[] psim;
        delete[] pobs;
        delete[] temp;
    }
}
void GoodOfFit::init(double *xsim, double *xobs, int n){
    nLen = n;
    sim = new double[nLen];
    obs = new double[nLen];
    ssim = new double[nLen];
    sobs = new double[nLen];
    psim = new double[nLen];
    pobs = new double[nLen];
    temp = new double[nLen];
    for(int i = 0; i < n; i++){
        sim[i] = xsim[i];
        obs[i] = xobs[i];
    }
    initflag = 1;
}
void GoodOfFit::callFDC(){
    FDC(sim, ssim, psim, nLen);
    FDC(obs, sobs, pobs, nLen);
}
double GoodOfFit::sum_sq_diff(double *x, double *y, int n){
    double ret = 0.;
    for (int i = 0; i < n; i++){
        ret += (x[i] - y[i]) * (x[i] - y[i]);
    }
    return ret;
}
double GoodOfFit::sum_sq_diff(double *x, double y, int n){
    double ret = 0.;
    for (int i = 0; i < n; i++){
        ret += (x[i] - y) * (x[i] - y);
    }
    return ret;
}
double GoodOfFit::which_q(double *sortx, double *exceedance, double q){
    int id = 0;
    for(int i = 0; i < nLen; i++){
        if(exceedance[i] >= q){
            id = i;
        }else{
            break;
        }
    }
    return sortx[id];
}
double GoodOfFit::max(double *x, int n){
    double ret = x[0];
    for (int i = 1; i < n; i++){
        ret = ret > x[i] ? ret : x[i];
    }
    return ret;
}
double GoodOfFit::min(double *x, int n){
    double ret = x[0];
    for (int i = 1; i < n; i++){
        ret = ret < x[i] ? ret : x[i];
    }
    return ret;
}
double GoodOfFit::mean(double *x, int n){
    /* Mean of X */
    double ret = 0.;
    ret = sum(x, n);
    return ret / n;
}
void GoodOfFit::product(double *x, double *y, int n, double *ret){
    for(int i = 0; i < n; i++){
        ret[i] = x[i] * y[i];
    }
}
double GoodOfFit::sum(double *x, int n){
    /* sum of X*/
    double ret = 0.;
    for (int i = 0; i < n; i++){
        ret += x[i];
    }
    return ret;
}
double GoodOfFit::sd(double *x, int n){
    /* sum of X*/
    double ret = 0.;
    double xmean = mean(x, n);
    ret = sum_sq_diff(x, xmean, n);
    ret = sqrt( ret / n );
    return ret;
}

/*********************GOOD OF FIT********************************/
double GoodOfFit::gof_NSE(){
    /***********************************************************************
     Nash-sutcliffe Efficiency
     Nash-Sutcliffe efficiencies (Nash and Sutcliffe, 1970) range from -Inf to 1.
     An efficiency of 1 (NSE = 1) corresponds to a perfect match of modeled to the observed data.
     An efficiency of 0 (NSE = 0) indicates that the model predictions are as accurate
     as the mean of the observed data, whereas
     an efficiency less than zero (-Inf < NSE < 0) occurs when the observed mean is a better predictor than the model.
     Essentially, the closer the model efficiency is to 1, the more accurate the model is.
     ********************************************************************************/
    double ret = 0.;
    double obsmean = mean(obs, nLen);
    double denominator = sum_sq_diff(obs, obsmean, nLen);
    double sum2 = sum_sq_diff(obs, sim, nLen);
    ret = 1. - ( sum2 / denominator );
    if (ret == NA_VALUE) {
        fprintf(stderr, "\nError of NSE: nLen = %d obsMean = %.2f, denominator = %.2f, sumqdiff(o,s)=sum2 \n", nLen, denominator, sum2);
    }
    return ret;
}
double GoodOfFit::gof_RMSE(){
    /* Root Mean Square Error */
    double ret = 0.;
    ret = gof_MSE();
    ret = sqrt(ret);
    return ret;
}
double GoodOfFit::gof_nRMSE(int imean){
    /* Normalized Root Mean Square Error */
    double ret = 0.;
    double cte = 0.;
    if(imean){
        cte = mean(obs, nLen);
    }else{
        cte = max(obs, nLen) - min(obs, nLen);
    }
    ret = gof_MSE();
    ret = sqrt(ret) / cte;
    return ret;
}
double GoodOfFit::gof_MSE(){
    /* Mean Square Error */
    double ret = sum_sq_diff(sim, obs, nLen);
    ret = ret / nLen;
    return ret;
}
double GoodOfFit::gof_ME(){
    /* Mean Error*/
    double ret = 0.;
    for (int i = 0; i < nLen; i++){
        ret += (sim[i] - obs[i]) ;
    }
    ret = ret / nLen;
    return ret;
}
double GoodOfFit::gof_MAE(){
    /* Mean Absolute Error*/
    double ret = 0.;
    for (int i = 0; i < nLen; i++){
        ret += fabs(sim[i] - obs[i]) ;
    }
    ret = ret / nLen;
    return ret;
}
double GoodOfFit::gof_pBias(){
    /* Percent Bias  */
    double ret = 0.;
    double denominator = sum(obs, nLen);
    for (int i = 0; i < nLen; i++){
        ret += (sim[i] - obs[i]) ;
    }
    ret = ret / denominator * 100.;
    return ret;
}
double GoodOfFit::gof_pBiasFDC(double lq, double hq){
    /* Percent Bias of the midsegment of the Flow Duration Curve  */
    double ret = 0.;
    double lq_obs, hq_obs;
    double lq_sim, hq_sim;
    lq_sim = which_q(ssim, psim, lq);
    hq_sim = which_q(ssim, psim, hq);
    lq_obs = which_q(sobs, pobs, lq);
    hq_obs = which_q(sobs, pobs, hq);
    double denominator = log(hq_obs) - log(lq_obs);
    ret = 100. * ( ( ( log(hq_sim) -  log(lq_sim) ) / denominator ) - 1. );
    return ret;
}

double GoodOfFit::gof_pBiasFDC(){
    /* Percent Bias of the midsegment of the Flow Duration Curve  */
    double ret = 0.;
    ret = gof_pBiasFDC(0.2, 0.8);
    return ret;
}
double GoodOfFit::gof_R2(){
    /* Pearson correlation */
    double ret = 0.;
    double obsmean = mean(obs, nLen);
    double simmean = mean(sim, nLen);
    product(obs, sim, nLen, temp);
    double osmean = mean(temp, nLen);
    double covariance = osmean - (obsmean * simmean);
    double sd1 = sd(obs, nLen);
    double sd2 = sd(sim, nLen);
    
    ret = covariance / (sd1 * sd2);
    ret = ret * ret;
    return ret;
}
double GoodOfFit::gof_ssq(){
    /* Percent Bias  */
    double ret = sum_sq_diff(sim, obs, nLen);
    double denominator = sum(obs, nLen);
    ret = ret / denominator * 100.;
    return ret;
}
void GoodOfFit::print_gofname(FILE *fp){
    fprintf(fp, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
            "NSE",
            "RMSE",
            "nRMSE",
            "MSE",
            "pBias",
            "pBiasFDC",
            "R2",
            "ME",
            "MAE",
            "ssq"
            );
}
void GoodOfFit::print_gof(FILE *fp){
    fprintf(fp, "%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
            gof_NSE(),
            gof_RMSE(),
            gof_nRMSE(1),
            gof_MSE(),
            gof_pBias(),
            gof_pBiasFDC(),
            gof_R2(),
            gof_ME(),
            gof_MAE(),
            gof_ssq()
            );
}
