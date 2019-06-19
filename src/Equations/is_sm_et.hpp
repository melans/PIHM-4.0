//
//  is_sm_et.h
//  pihm2.2
//
//  Created by Lele Shu on 6/25/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#ifndef is_sm_et_h
#define is_sm_et_h
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "nvector_serial.h"
#include "sundials_types.h"
#include "ModelConfigure.hpp"
#include "Macros.hpp"
#include "functions.hpp"

//void is_sm_et(double t, double stepsize, void *DS, N_Vector VY);

double Penman_Monteith(double T, double RH, double Vel,
                       double P,double Rn, double rl,
                       double windH, double Gamma,
                       double R_ref, double lai);
double SoilMoistureStress(double ThetaS, double ThetaR, double SatRatio);
double Eact_Interception(double etp, double yis, double ymax);
double Eact_Vegetation(double etp, double yis, double ymax);
double ActualEvaporation(double etp, double ThetaS, double ThetaR, double SatRatio);

double VaporPressure_Sat(double T_in_C);
double VaporPressure_Act(double esat, double rh);
//double VaporPressure_Act(double rh, double T_C, double P);
double AerodynamicResistance(double v, double rl, double height);
double SlopeSatVaporPressure(double T);
double CanopyResistance(double rmin, double lai);
double PsychrometricConstant(double Pressure);
double PressureElevation(double z);
double BulkSurfaceResistance(double R_ref, double lai);



inline double Eact_Interception(double etp, double yis, double ymax){
    return etp * pow23(yis / ymax); // Deardorff (1978)
    /* Mukeksh(2009) eq 4.7c*/
}
inline double Eact_Vegetation(double etp, double yis, double ymax){
    double aet=0;
    return aet;
}
inline double ActualEvaporation(double etp, double ThetaS, double ThetaR, double SatRatio){
    return SoilMoistureStress(ThetaS, ThetaR, SatRatio) * etp;
}

inline double BulkSurfaceResistance(double R_ref, double lai){
    // R_ref     bulk stomatal resistance of the well-illuminated leaf [min m-1],
    return R_ref * 2. / lai; /* Allen(1998) */
}
inline double PressureElevation(double z){
    return 101.325 * pow((293. - 0.0065 * z) / 293, 5.26); /*Pressure based on Elevation*/
    //    P atmospheric pressure [kPa],
    //    z elevation above sea level [m],
    /* Allen(1998) Eq(7) */
}
inline double PsychrometricConstant(double Pressure){
    return 0.665e-3 * Pressure;
    //    γ psychrometric constant [kPa °C-1],
    /* Allen(1998) Eq(8) */
}
inline double VaporPressure_Sat(double T_in_C){
    //    e°(T) saturation vapour pressure at the air temperature T [kPa],
    //    T air temperature [°C],
    //    exp[..] 2.7183 (base of natural logarithm) raised to the power [..].
    //Zotarelli, L., & Dukes, M. (2010).
    // Allen(1998) Eq (11)
    return .61128 * exp( 17.27 * T_in_C / (T_in_C + 237.3));
    //    http://www.fao.org/docrep/X0490E/x0490e07.htm
}
inline double VaporPressure_Act(double esat, double rh){
    return esat * rh;
}
//double VaporPressure_Act(double RH, double TC, double P){
//    double VP = 611.2 * exp(17.67 * TC / (TC + 243.5)) * RH;
//    double qv_sat = 0.622 * (VP / RH) / P;
//    return qv_sat;
//}
inline double AerodynamicResistance(double Vel, double rl, double z){
    /* Allen, R. G., S, P. L., Raes, D., & Martin, S. (1998).
     Crop evapotranspiration : Guidelines for computing crop water requirements
     by Richard G. Allen ... [et al.].
     FAO irrigation and drainage paper: 56.
     https://doi.org/10.1016/j.eja.2010.12.001 */
    // Eq(4) in reference
    //    ra aerodynamic resistance [s m-1],
    //    zm height of wind measurements [m],
    //    zh height of humidity measurements [m],
    //    d zero plane displacement height [m],
    //    zom roughness length governing momentum transfer [m],
    //    zoh roughness length governing transfer of heat and vapour [m],
    //    k von Karman's constant, 0.41 [-],
    //    uz wind speed at height z [m s-1].
    //    r_a = 12 * 4.72 * log(Ele[i].windH / rl) / (0.54 * Vel / UNIT_C / 60 + 1) / UNIT_C / 60;    return r_a;
    double  r_a;
    double  zw, zp, rw, rp;
    zw = z;
    zp = 2.0;
    rw = rl;
    rp = rl * 0.1;
    //    r_a = 12. * 4.72 * log(z / rl) / (0.54 * Vel);
    r_a = log( zw / rw ) * log( zp / rp)
    / (VON_KARMAN * VON_KARMAN* Vel); /* bug?? */
    return r_a;
}

inline double SlopeSatVaporPressure(double T){
    //# Slope of saturation vapour pressure curve (D )
    //# http://www.fao.org/docrep/X0490E/x0490e07.htm
    //# D slope of saturation vapour pressure curve at air temperature T [kPa °C-1],
    //# T air temperature [°C],
    //# exp[..] 2.7183 (base of natural logarithm) raised to the power [..].
    double tt = (T + 237.3);
    double delta = 4098. * ( 0.6108* exp(17.27 * T / tt )  ) /  ( tt * tt);
    return delta;
}

inline double CanopyResistance(double rmin, double lai){
    double rc;
    rc = rmin * 2. / lai;
    return rc;
}


#endif /* is_sm_et_h */