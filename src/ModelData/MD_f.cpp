//
//  MD_f.cpp
//  PIHM++ (v 4.0)
//
//  Created by Lele Shu on 1/27/19.
//  Copyright © 2019 Lele Shu. All rights reserved.
//

#include "MD_f.hpp"
void Model_Data:: f_loop( double  *Y, double  *DY, double t){
    int i;
    tnow = t;
    for (i = 0; i < NumEle; i++) {
        /*DO INFILTRATION FRIST, then do LATERAL FLOW.*/
        /*========ET Function==============*/
        //        f_etFlux(i, t); // et is moved out of f_loop()
        /*========infiltration/Recharge Function==============*/
        Ele[i].updateElement(uYsf[i] , uYus[i] , uYgw[i] ); // step 1 update the kinf, kh, etc. for elements.
        f_InfilRecharge(i, t); // step 2 calculate the infiltration and recharge.
    }
    for (i = 0; i < NumEle; i++) {
        /*DO INFILTRATION FRIST, then do LATERAL FLOW.*/
        /*========ET Function==============*/
        //        f_etFlux(i, t); // et is moved out of f_loop()
        /*========infiltration/Recharge Function==============*/
        Ele[i].updateElement(uYsf[i] , uYus[i] , uYgw[i] ); // step 1 update the kinf, kh, etc. for elements.
        f_InfilRecharge(i, t); // step 2 calculate the infiltration and recharge.
        /*========surf/gw flow Function==============*/
        f_lateralFlux(i, t); // AFTER infiltration, do the lateral flux. ESP for overland flow.
//        CheckNANi(qEleInfil[i], i, "qEleInfil[i]::f_recharge()");
//        CheckNANi(qEleRecharge[i], i, "qEleRecharge[i]::f_recharge()");
    } //end of for loop.
    for (i = 0; i < NumSegmt; i++) {
        f_Segement_surface(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
        f_Segement_sub(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
    }
    for (i = 0; i < NumRiv; i++) {
        Flux_RiverDown(t, i);
    }
    /* Shared for both OpenMP and Serial, to update */
    PassValue();
}
void Model_Data::f_applyDY(double *DY, double t){
    double area;
    for (int i = 0; i < NumEle; i++) {
        area = Ele[i].area;
        QeleSurfTot[i] = Qe2r_Surf[i];
        QeleSubTot[i] = Qe2r_Sub[i];
        for (int j = 0; j < 3; j++) {
            QeleSurfTot[i] += QeleSurf[i][j];
            QeleSubTot[i] += QeleSub[i][j];
//            CheckNANi(QeleSubTot[i], 1, "QeleSubTot[i]");
//            CheckNANi(QeleSurfTot[i], 1, "QeleSurfTot[i]");
        }
        DY[i] = qEleNetPrep[i] - qEleInfil[i] + qEleExfil[i] - QeleSurfTot[i] / area;
        DY[iUS] = qEleInfil[i] - qEleRecharge[i];
        DY[iGW] = qEleRecharge[i] - qEleExfil[i] - QeleSubTot[i] / area;
        
        if(uYsf[i] < EPSILON){ /* NO ponding water*/
            if (uYgw[i] < Ele[i].WetlandLevel){
                /*Evaporate from unsat soil*/
                DY[iUS] += -qEleET[i][2];
            }else{
                /*Evaporate from Ground water*/
                DY[iGW] += -qEleET[i][2];
            }
        }else{ /*Ponding water*/
            DY[i] +=  - qEleET[i][2];
        }
        if (uYgw[i] > Ele[i].RootReachLevel) {
            /*Vegetation sucks water from Ground water*/
            DY[iGW] += - qEleET[i][1];
        } else {
            DY[iUS] += - qEleET[i][1];
        }
        
        /* Boundary condition and Source/Sink */
        if(Ele[i].iBC > 0){ // Fix head of GW.
            DY[iGW] = 0;
        }else if(Ele[i].iBC < 0){ // Fix flux in GW
            DY[iGW] += Ele[i].QBC / area;
        }else{ /* Void */}
        
        if(Ele[i].iSS > 0){ // SS in Landusrface
            DY[iSF] += Ele[i].QSS / area;
        }else if(Ele[i].iSS < 0){ // SS in GW
            DY[iGW] += Ele[i].QSS / area;
        }else{}
        
        /* Convert with specific yield */
        DY[iUS] /= Ele[i].Sy;
        DY[iGW] /= Ele[i].Sy;
//                DY[iSF] =0.0;  // debug only.
//                DY[iUS] =0.0;
//                DY[iGW] =0.0;
//        printf("%d:%e\t%e\t%e\t%e\t%e\t%e\t%e\n", i+1, qEleNetPrep[i], qEleInfil[i], qEleExfil[i], QeleSurfTot[i] / area);
#ifdef _DEBUG
        CheckNANi(DY[i], i, "DY[i] (Model_Data::f_applyDY)");
        CheckNANi(DY[iUS], i, "DY[iUS] (Model_Data::f_applyDY)");
        CheckNANi(DY[iGW], i, "DY[iGW] (Model_Data::f_applyDY)");
#endif
    }
    for (int i = 0; i < NumRiv; i++) {
        if(Riv[i].BC > 0){
//            Newmann condition.
            DY[iRIV] = 0.;
        }else{
            DY[iRIV] = (- QrivUp[i] - QrivSurf[i] - QrivSub[i] - QrivDown[i] + Riv[i].qBC) / Riv[i].u_TopArea;
        }
//        DY[iRIV] = 0.0;
#ifdef _DEBUG
        CheckNANi(DY[i + 3 * NumEle], i, "DY[i] of river (Model_Data::f_applyDY)");
#endif
    }
}

void Model_Data::f_Segement_update(int iEle, int iRiv, int i){
    QrivSurf[iRiv] += QsegSurf[i]; // Positive from River to Element
    QrivSub[iRiv] += QsegSub[i];
    Qe2r_Surf[iEle] += -QsegSurf[i]; // Positive from Element to River
    Qe2r_Sub[iEle] += -QsegSub[i];
    CheckNANi(QrivSub[i], i, "xxx");
    CheckNANi(QrivSurf[i], i, "xxx");
}
void Model_Data::PassValue(){
    int i, j, inabr, jnabr, ie, ir;
    for (i = 0; i < NumEle; i++) { /*Check flux A->B  = Flux B->A*/
        for (j = 0; j < 3; j++) {
            inabr = Ele[i].nabr[j] - 1;
            if (inabr >= 0) {
                jnabr = Ele[i].nabrToMe[j];
                if(Ele[inabr].iupdSF[jnabr]){
                    //void
                }else{
                    QeleSurf[inabr][jnabr] = - QeleSurf[i][j];
                    Ele[inabr].iupdSF[jnabr] = 1;
                    QeleSub[inabr][jnabr] = - QeleSub[i][j];
                    Ele[inabr].iupdGW[jnabr] = 1;
                }
            }
        }
    }
    for (i = 0; i < NumSegmt; i++) {
        ie = RivSeg[i].iEle-1;
        ir = RivSeg[i].iRiv-1;
//        f_Segement_update(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
        QrivSurf[ir] += QsegSurf[i]; // Positive from River to Element
        QrivSub[ir] += QsegSub[i];
        Qe2r_Surf[ie] += -QsegSurf[i]; // Positive from Element to River
        Qe2r_Sub[ie] += -QsegSub[i];
    }
    for (i = 0; i < NumRiv; i++) {
        if(iDownStrm >= 0){
            QrivUp[iDownStrm] += - QrivDown[i];
        }
    }
}

void Model_Data::applyBCSS(double *DY, int i){
    if(Ele[i].iBC > 0){ // Fix head of GW.
        DY[iGW] = 0;
    }else if(Ele[i].iBC < 0){ // Fix flux in GW
        DY[iGW] += Ele[i].QBC / Ele[i].area;
    }else{}
    
    if(Ele[i].iSS > 0){ // SS in Landusrface
        DY[iSF] += Ele[i].QSS / Ele[i].area;
    }else if(Ele[i].iSS < 0){ // SS in GW
        DY[iGW] += Ele[i].QSS / Ele[i].area;
    }else{}
}
