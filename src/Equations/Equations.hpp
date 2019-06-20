//  Equations.hpp
//  Created by Lele Shu (lele.shu@gmail.com) on 2018.
//  Copyright © 2018 Lele Shu. All rights reserved.
//
#ifndef Equations_hpp
#define Equations_hpp

#include <stdio.h>
#include <math.h>
#include "Macros.hpp"
#include "functions.hpp"


//double returnVal(double rArea, double rPerem, double eqWid, double ap_Bool);
//double CS_AreaOrPerem(int rivOrder, double rivDepth, double rivCoeff, double a_pBool);
double avgY(double z1, double y1, double z2, double y2, double threshold);
double effKV(double ksatFunc, double gradY, double macKV, double KV, double areaF);
double effKVnew(double ksatFunc, double macKV, double KV, double areaF, double y0);

double effKH(double tmpY, double aqDepth, double MacD, double MacKsatH, double areaF, double ksatH);
double satKfun(double elemSatn, double beta);
double sat2psi(double elemSatn, double alpha, double beta, double mpsi);
double fun_recharge(double effk_us, double kgw, double Deficit, double ygw, double yus);
double effKRech(double ksatFunc,  double macKV, double KV, double areaF);
double ManningEquation(double Area, double rough, double y, double S);
double OverlandManning(double avg_y, double grad_y, double avg_sf, double A, double n);
double FieldCapacity (double Alpha, double Beta, double deficit );
double GreenAmpt(double k, double ti, double ts, double phi, double hf, double prcp, double sy);


inline double pow23(double x){
    double t = cbrt(x);
    return t * t;
}
inline double sqpow2(double x, double y){
    return sqrt(x * x + y * y);
}
inline double meanHarmonic(double k1, double k2, double d1, double d2){
    return (d1 + d2) / ( d1 / k1 + d2 / k2);
}
inline double meanArithmetic(double k1, double k2, double d1, double d2){
    return (k1 * d1 + k2 * d2) / (d1 + d2);
}
#endif /* Equations_hpp */
