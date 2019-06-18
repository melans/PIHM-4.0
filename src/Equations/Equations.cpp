//
//  Equations.cpp
//  Created by Lele Shu on 8/13/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#include "Equations.hpp"
double avgY(double z1, double y1, double z2, double y2, double threshold){
    double h1 = z1 + y1, h2 = z2 + y2;
    if (h1 > h2) {
        if (y1 > threshold) {
            //return 0.5 * (yi + yinabr);
            //return ((yinabr > yi) ? 0 : 1.0 * yi);
            /* Note the if-else TRUE case can be possible only for
             * Kinematic case */
            return y1;
        } else {
            return 0.;
        }
    } else {
        if (y2 > threshold) {
            //return 0.5 * (yi + yinabr);
            //return ((yi > yinabr) ? 0 : 1.0 * yinabr);
            /* Note the if-else TRUE case can be possible only for
             * Kinematic case */
            return y2;
        } else {
            return 0.;
        }
    }
}
double effKV(double ksatFunc, double gradY, double macKV, double KV, double areaF)
{
    if (ksatFunc >= 0.98) {
        return (macKV * areaF + KV * (1 - areaF) * ksatFunc);
    } else {
        if (fabs(gradY) * ksatFunc * KV <= 1 * KV * ksatFunc) {
            return KV * ksatFunc;
        } else {
            if (fabs(gradY) * ksatFunc * KV < (macKV * areaF + KV * (1 - areaF) * ksatFunc)) {
                return (macKV * areaF * ksatFunc + KV * (1 - areaF) * ksatFunc);
            } else {
                return (macKV * areaF + KV * (1 - areaF) * ksatFunc);
            }
        }
    }
}
double effKRech(double ksatFunc,  double macKV, double KV, double areaF)
{
    double keff = 0.0;
    if(areaF <= 0. ){
        keff = KV * ksatFunc;
    }else{
        keff = (KV * (1 - areaF) + macKV * areaF) * ksatFunc;
    }
    return keff;
}
double effKVnew(double ksatFunc,  double macKV, double KV, double areaF, double y0)
{
    double keff = 0.0;
    
    if(areaF <= 0. ){
        keff = KV * ksatFunc;
    }else{
        keff = KV * (1 - areaF) + macKV * areaF;
    }
//    if(areaF <= 0. ){
//        keff = KV * ksatFunc;
//    }else if (y0 > EPSILON){
//        /*  Ponding water exists */
//        keff = KV * (1 - areaF) + macKV * areaF;
//    }else{
//        keff = (KV * (1 - areaF) + macKV * areaF) * ksatFunc;
//    }
    return keff;
}
double effKH(double Ygw, double aqDepth, double MacD, double Kmac, double AF, double Kmx){
    double effk = 0;
    if (MacD <= EPS_DOUBLE || Ygw < aqDepth - MacD) {
        effk = Kmx;
    } else {
        if (Ygw > aqDepth) {
            effk = (Kmac * MacD * AF +
                    Kmx * (aqDepth - MacD * AF)) / aqDepth;
        } else {
            effk = (Kmac * (Ygw - (aqDepth - MacD)) * AF +
                    Kmx * (aqDepth - MacD + (Ygw - (aqDepth - MacD)) * (1 - AF))) / Ygw;
        }
    }
    if( effk < 0. || effk > 1e9 ){
        fprintf(stderr, "Wrong effKH for ground water: %f [ m/s ]\n", effk / 60.);
        myexit(ERRDATAIN);
    }
    return effk;
}

double satKfun(double elemSatn, double n){
    double temp = -1 + pow(1 - pow(elemSatn, n / (n - 1)), (n - 1) / n);
    double ret = sqrt(elemSatn) * temp * temp;
    return ret;
}
double ManningEquation(double Area, double rough, double R, double S){
    double Q = 0.;
    if (S > 0) {
        Q = sqrt(S) * Area * pow23(R) / rough;
        //return sqrt(S) * Area * pow(Area / Perem, 2. / 3.) / rough;
    } else {
        Q = -1.0 * sqrt(-S) * Area * pow23(R) / rough;
    }
    return Q;
}
double OverlandManning(double avg_y, double grad_y, double avg_sf, double A, double n){
    double Q;
    //flux[loci][locj] = crossA * pow23(avg_y) * grad_y / (sqrt(fabs(avg_sf)) * avg_rough);
    Q = A * pow23(avg_y) * grad_y / ( sqrt(fabs(avg_sf)) * n );
//    Q = A * pow23(avg_y) * grad_y / ( sqrt(fabs(avg_sf)) * n );
//        if (grad_y >= 0.) {
//            return A * pow23(avg_y) * sqrt(grad_y) / n;
//        } else {
//            return -1. * A * pow23(avg_y) * sqrt(-1. * grad_y) / n;
//        }
    return Q;
}
double sat2Y(double elemSatn, double alpha, double n, double mpsi){
    double temp = -(pow(pow(elemSatn, n / (1 - n)) - 1, 1 / n) / alpha);
    double ret = (temp < mpsi) ? mpsi : temp;
    return ret;
}
double fun_recharge(double effk_us, double kgw, double Deficit, double ygw, double hus)
{
    double effk;
    double q = 0.;
    if(effk_us < EPS_DOUBLE || kgw < EPS_DOUBLE){
        return 0.;
    }
    effk = meanHarmonic(kgw, effk_us, ygw, hus);
    q = effk * (1. + hus / (Deficit * 0.5));
    return q;
}
double FieldCapacity (double Alpha, double Beta, double deficit ){
    //    sfc=((alpha *deficit)^(-beta) +1) ^(-1/beta)/alpha
    double    sfc;
    sfc = pow( pow(Alpha * deficit, -Beta) +1, -1/Beta) / Alpha;
    return sfc;
}
double GreenAmpt(double k, double ti, double ts, double phi, double hf, double h0, double Sy){
    double dTheta = ts - ti;
//    double dh = h0 + hf - phi;
    double q = 0.;
//    double D0;
    if(h0 <= 0 ){
        /* No rainfall, no ponding */
        return 0.;
    }
    hf = max(EPSILON, hf);
    if(dTheta <= 0.){
        q = 0.;
    }else{
        q = k * ( (h0 + hf - phi) * dTheta / hf  );
    }
    if( q >= 0. ){
        if( q > h0 ){
            q = h0;
        }else{
            q = q;
        }
    }else{
        q = q;
    }
    return q;
}

