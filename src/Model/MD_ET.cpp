//
//  MD_ET.cpp
//  PIHM++
//
//  Created by Lele Shu on 10/26/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#include "Model_Data.hpp"
void Model_Data::updateforcing(double t){
    for (int i = 0; i < NumForc; i++){
        tsd_weather[i].movePointer(t);
    }
    tsd_MF.movePointer(t);
    tsd_LAI.movePointer(t);
    tsd_RL.movePointer(t);
    for(int i = 0; i < NumEle; i++){
        tReadForcing(t,i);
    }
}
void Model_Data::tReadForcing(double t, int i){
    t_prcp[i] = tsd_weather[Ele[i].iForc - 1].getX(t, i_prcp) * gc.cPrep;
    t_temp[i] = tsd_weather[Ele[i].iForc - 1].getX(t, i_temp) +  gc.cTemp;
    t_lai[i] = tsd_LAI.getX(t, Ele[i].iLC) * gc.cLAItsd ;
    t_mf[i] = tsd_MF.getX(t, Ele[i].iMF);
    t_rn[i] = tsd_weather[Ele[i].iForc - 1].getX(t, i_rn) * (1 - Ele[i].Albedo);
    t_wind[i] = (fabs(tsd_weather[Ele[i].iForc - 1].getX(t, i_wind) ) + 1); // +1 voids ZERO.
    t_rh[i] = tsd_weather[Ele[i].iForc - 1].getX(t, i_rh);
    t_rl[i] = tsd_RL.getX(t, Ele[i].iLC);
    
    t_rh[i] = min(max(t_rh[i], 0.01), 1.0);
    
    t_prcp[i]   = t_prcp[i] / 1440.; // debug
    t_rn[i]     = t_rn[i] * 1.0e-6 / 1440.;  // J/m2/day to [MJ m-2 min-1]
    t_wind[i]   = t_wind[i] / 1440. ; // Unit = m/s ; m/d to m/min
    
    qElePrep[i] = t_prcp[i];
    
    //    qEleETP[i] = Penman_Monteith(t_temp[i], t_rh[i], t_wind[i], Ele[i].FixPressure, t_rn[i], t_rl[i], Ele[i].windH);
    qEleETP[i] = Penman_Monteith(t_temp[i], t_rh[i], t_wind[i],
                                 Ele[i].FixPressure, t_rn[i], t_rl[i],
                                 Ele[i].windH, Ele[i].FixGamma,
                                 Ele[i].Rmin, t_lai[i]);
    //    if(i==1){
    //            printf("%.1f@%5d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", t, i+1, t_prcp[i], t_temp[i], t_rh[i], t_rn[i] , t_wind[i], t_lai[i] , t_rl[i], t_mf[i], qElePrep[i], qEleETP[i]);
    //        printf("etp=%E\n", qEleETP[i]);
    //    }
}

void Model_Data::EvapoTranspiration(N_Vector uY, double t, double dt){
//    double  Rn=NA_VALUE, Vel=NA_VALUE, RH=NA_VALUE, rl=NA_VALUE;
    double  T=NA_VALUE,  LAI=NA_VALUE;
    double  isval = 0.;// etval = 0;
    double  DT_min = NA_VALUE;
    double  r_ISMax = 0.0;
    double  fracSnow = NA_VALUE,
    snowRate=NA_VALUE,
    MeltRateGrnd=NA_VALUE,
    MeltRateCanopy=NA_VALUE,
    MF=NA_VALUE,
    ret=NA_VALUE;
    DT_min = dt ; /* dt [min] */
    
    for(int i = 0; i < NumEle; i++) {
        /* Note the dependence on physical units */
        T = t_temp[i];
        LAI = t_lai[i];
        MF = t_mf[i];
//        Rn = t_rn[i];
//        Vel = t_wind[i];
//        RH = t_rh[i];
//        rl = t_rl[i];
        /* Snow Accumulation/Melt Calculation*/
        if( T < Ts){
            fracSnow = 1.;
        }else if(T > Tr){
            fracSnow = 0.;
        }else{
            fracSnow = (Tr - T) / (Tr - Ts);
        }
        snowRate = fracSnow * qElePrep[i];
        /* EleSnowGrnd, yEleSnowCanopy, yEleISsnowmax,
         * MeltRateGrnd,MeltRateCanopy are the average value prorated
         * over the whole elemental area */
        yEleSnowGrnd[i] += (1 - Ele[i].VegFrac) * snowRate * DT_min;
        yEleSnowCanopy[i] += Ele[i].VegFrac * snowRate * DT_min;
        yEleISsnowmax[i] = yEleSnowCanopy[i] > 0 ? 0.003 * LAI * Ele[i].VegFrac : 0;
        yEleISsnowmax[i] = yEleISsnowmax[i];
        if (yEleSnowCanopy[i] > yEleISsnowmax[i]) {
            yEleSnowCanopy[i] = yEleISsnowmax[i];
            yEleSnowGrnd[i] = yEleSnowGrnd[i] + yEleSnowCanopy[i] - yEleISsnowmax[i];
        }
        MeltRateGrnd = MeltRateCanopy = (T > To ? (T - To) * MF : 0.);    /* Note the units for
                                                                           * MF. */
        if (yEleSnowGrnd[i] > MeltRateGrnd * DT_min) {
            yEleSnowGrnd[i] = yEleSnowGrnd[i] - MeltRateGrnd * DT_min;
        } else {
            MeltRateGrnd = yEleSnowGrnd[i] / DT_min;
            yEleSnowGrnd[i] = 0;
        }
        if (yEleSnowCanopy[i] > MeltRateCanopy * DT_min) {
            yEleSnowCanopy[i] = yEleSnowCanopy[i] - MeltRateCanopy * DT_min;
        } else {
            MeltRateCanopy = yEleSnowCanopy[i] / DT_min;
            yEleSnowCanopy[i] = 0;
        }
        yEleSnow[i] = yEleSnowCanopy[i] + yEleSnowGrnd[i];
        /************************************************************************/
        /* ThroughFall and Evaporation from canopy             */
        /************************************************************************/
        /*
         * EleIS, EleET[0] and ret are prorated for the whole
         * element. Logistics are simpler if assumed in volumetric
         * form by multiplication of Area on either side of equation
         */
        yEleISmax[i] = gc.cISmax * 0.0002 * LAI * Ele[i].VegFrac;
        /* Note the dependence on physical units */
        yEleIS[i] = min(yEleIS[i], yEleISmax[i]);
        r_ISMax = yEleIS[i] / yEleISmax[i];
        if (LAI > 0.0 && Ele[i].VegFrac > 0.0) {
            if(yEleIS[i] < 0){
                qEleET[i][0] = 0.;
            }else{
                qEleET[i][0] = gc.cEt0 * Ele[i].VegFrac * sqrt(r_ISMax) * qEleETP[i];
            }
            qEleET[i][0] = qEleET[i][0] < 0 ? 0 : qEleET[i][0];
            
            if(yEleIS[i] <= 0 ){
                qEleTF[i] = 0.;
            }else{
                qEleTF[i] = 0.0565 * yEleISmax[i] * exp(3.89 * r_ISMax);
            }
        } else {
            qEleET[i][0] = 0.0;
            qEleTF[i] = 0.0;
        }
        if (yEleIS[i] >= yEleISmax[i]) {
            if (((1 - fracSnow) * qElePrep[i] * Ele[i].VegFrac + MeltRateCanopy) >= qEleET[i][0] + qEleTF[i]) {
                qEleETloss[i] = qEleET[i][0];
                ret = qEleTF[i] + yEleIS[i] - yEleISmax[i] + (((1 - fracSnow) * qElePrep[i] * Ele[i].VegFrac + MeltRateCanopy) - (qEleET[i][0] + qEleTF[i]));
                isval = yEleISmax[i];
                //EleIS[i] = yEleISmax[i];
            } else if ((((1 - fracSnow) * qElePrep[i] * Ele[i].VegFrac + MeltRateCanopy) < qEleET[i][0] + qEleTF[i]) && (yEleIS[i] + DT_min * ((1 - fracSnow) * qElePrep[i] * Ele[i].VegFrac + MeltRateCanopy - qEleET[i][0] - qEleTF[i]) <= 0)) {
                qEleET[i][0] = (qEleET[i][0] / (qEleET[i][0] + qEleTF[i])) * (yEleIS[i] / DT_min + ((1 - fracSnow) * qElePrep[i] * Ele[i].VegFrac + MeltRateCanopy));
                ret = (qEleTF[i] / (qEleET[i][0] + qEleTF[i])) * (yEleIS[i] / DT_min + ((1 - fracSnow) * qElePrep[i] * Ele[i].VegFrac + MeltRateCanopy));
                //EleIS[i] = 0;
                isval = 0;
                qEleETloss[i] = qEleET[i][0];
            } else {
                isval = yEleIS[i] + DT_min * (((1 - fracSnow) * qElePrep[i] * Ele[i].VegFrac + MeltRateCanopy) - qEleET[i][0] - qEleTF[i]);
                //EleIS[i] = EleIS[i] + dt * (((1 - fracSnow) * ElePrep[i] * Ele[i].VegFrac + MeltRateCanopy) - EleET[i][0] - EleTF[i]);
                ret = qEleTF[i];
                qEleETloss[i] = qEleET[i][0];
            }
        } else if ((yEleIS[i] < yEleISmax[i]) && ((yEleIS[i] + (((1 - fracSnow) * qElePrep[i] * Ele[i].VegFrac + MeltRateCanopy) - qEleET[i][0] - qEleTF[i]) * DT_min) >= yEleISmax[i])) {
            qEleETloss[i] = qEleET[i][0];
            isval = yEleISmax[i];
            ret = qEleTF[i] + ((yEleIS[i] + (((1 - fracSnow) * qElePrep[i] * Ele[i].VegFrac + MeltRateCanopy) - qEleET[i][0] - qEleTF[i]) * DT_min) - yEleISmax[i]);
        } else if ((yEleIS[i] < yEleISmax[i]) && ((yEleIS[i] + (((1 - fracSnow) * qElePrep[i] * Ele[i].VegFrac + MeltRateCanopy) - qEleET[i][0] - qEleTF[i]) * DT_min) <= 0)) {
            if ((qEleET[i][0] > 0) || (qEleTF[i] > 0)) {
                qEleET[i][0] = (qEleET[i][0] / (qEleET[i][0] + qEleTF[i])) * (yEleIS[i] / DT_min + ((1 - fracSnow) * qElePrep[i] * Ele[i].VegFrac + MeltRateCanopy));
                ret = (qEleTF[i] / (qEleET[i][0] + qEleTF[i])) * (yEleIS[i] / DT_min + ((1 - fracSnow) * qElePrep[i] * Ele[i].VegFrac + MeltRateCanopy));
            } else {
                qEleET[i][0] = 0;
                ret = 0;
            }
            qEleETloss[i] = qEleET[i][0];
            isval = 0;
        } else {
            isval = yEleIS[i] + (((1 - fracSnow) * qElePrep[i] * Ele[i].VegFrac + MeltRateCanopy) - qEleET[i][0] - qEleTF[i]) * DT_min;
            qEleETloss[i] = qEleET[i][0];
            ret = qEleTF[i];
        }
        qEleNetPrep[i] = (1 - Ele[i].VegFrac) * (1 - fracSnow) * qElePrep[i] + ret + MeltRateGrnd;
        qEleTF[i] = ret;
        yEleIS[i] = isval;
#ifdef _DEBUG
        CheckNANi(qEleNetPrep[i], i, "qEleNetPrep[i]");
        CheckNANi(qElePrep[i], i, "qElePrep[i]");
        CheckNANi(qEleET[i][0], i, "qEleET[i][0]");
        CheckNANi(qEleET[i][1], i, "qEleET[i][1]");
        CheckNANi(qEleET[i][2], i, "qEleET[i][2]");
        CheckNANi(qEleTF[i], i, "qEleTF");
        CheckNANi(yEleIS[i], i, "yEleIS");
        CheckNANi(qEleETloss[i], i, "qEleETloss");
#endif
        f_etFlux(i, t);
    } // end of for i=1:NumEle
}

void Model_Data::f_etFlux(int i, double t){
    double  elemSatn, beta_s, LAI, ETp ;
    double  Et, Ev;
    ETp = qEleETP[i];
    LAI = t_lai[i];
    elemSatn = Ele[i].u_satn;
    beta_s = SoilMoistureStress(Soil[(Ele[i].iSoil - 1)].ThetaS, Soil[(Ele[i].iSoil - 1)].ThetaR, elemSatn);
    Ev = gc.cEt2 * (1 - Ele[i].VegFrac) * (1 - Ele[i].Landcover::ImpAF) * ETp * beta_s;
    Ev = max(0.0,  Ev);
    
    if(LAI > 0.){
        if(uYgw[i] > Ele[i].RootReachLevel){
            beta_s = 1.;
        }
        Et = gc.cEt2 * Ele[i].VegFrac * (1 - Ele[i].Landcover::ImpAF) * ETp * beta_s;
    }else{
        Et = 0.;
    }
    qEleET[i][1] = Et;
    qEleET[i][2] = Ev;
    qEleETA[i] = qEleET[i][0] + qEleET[i][1] + qEleET[i][2];
#ifdef _DEBUG
    CheckNANi(ETp, i, "Potential ET (Model_Data::EvapoTranspiration)");
    CheckNANi(qEleET[i][1], i, "Transpiration (Model_Data::EvapoTranspiration)");
    CheckNANi(qEleET[i][2], i, "Soil Evaporation (Model_Data::EvapoTranspiration)");
#endif
}

