#include "is_sm_et.hpp"

double Penman_Monteith(double T, double RH, double Vel, double P,
                       double Rn, double rl, double windH,
                       double Gamma,
                       double R_ref, double lai){
    /* http://www.fao.org/docrep/X0490E/x0490e06.htm#penman%20monteith%20equation
     Rn net radiation at the crop surface [MJ m-2 min-1],
     G soil heat flux density [MJ m-2 min-1],
     T mean daily air temperature at 2 m height [°C],
     u2 wind speed at 2 m height [m s-1],
     es saturation vapour pressure [kPa],
     ea actual vapour pressure [kPa],
     es - ea saturation vapour pressure deficit [kPa],
     ∆ slope vapour pressure curve [kPa °C-1],
     γ psychrometric constant [kPa °C-1].
    */
    double Delta;
    double  G = 0.;  // Soil heat flux. Assume G = 0;
    //    double    Gamma;  // psychrometric constant
    //    double    Delta;  // the slope of the saturation vapour pressure temperature relationship
    double    r_s;    // (bulk) surface resistances. [min m-1]
    double    r_a;    // aerodynamic resistances. [min m-1]
    double    ETp;    //
    double    ea, es, ed; //represents the vapour pressure deficit of the air
    double  rho = 1.225;    // Air density at sea level at 15C.
    Delta = SlopeSatVaporPressure(T);
    
    rho = 3.486 * P / (275 + T) ;  // David R Maidment, Hand of Hydrology Eq4.2.4.
    es = VaporPressure_Sat(T);
    ea = VaporPressure_Act(es, RH); //kpa
    ed = es - ea ; // [kpa]
    r_a = AerodynamicResistance(Vel, rl, windH); // [min m-1]
    r_s = BulkSurfaceResistance(R_ref, lai);
    ETp = (Delta * (Rn - G) + rho * Cp * ed / r_a)
    / (Delta + Gamma * (1 + r_s / r_a));
    ETp = ETp / LAMBDA * 0.001 ; // mm/min = m/min;   kg/m2/min = mm/min
    return ETp;
}
double SoilMoistureStress(double ThetaS, double ThetaR, double SatRatio){
    double fc, ThetaW, beta_s;
    fc = ThetaS * 0.75;
    ThetaW = ThetaR * 1.05;
    beta_s = (SatRatio * ThetaS + ThetaR - ThetaW) / (fc - ThetaW);
    if(beta_s > 1.){
        beta_s = 1.;
    }else if (beta_s < 0.){
        beta_s = 0.;
    }
    return beta_s;
}

