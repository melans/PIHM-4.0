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
#ifdef _PIHMOMP
#pragma omp parallel  default(shared) private(i) num_threads(CS.num_threads)
    {
#pragma omp for
        for (i = 0; i < NumEle; i++) {
            /*DO INFILTRATION FRIST, then do LATERAL FLOW.*/
            /*========ET Function==============*/
            //        f_etFlux(i, t); // et is moved out of f_loop()
            /*========infiltration/Recharge Function==============*/
            Ele[i].updateElement(uYsf[i] , uYus[i] , uYgw[i] ); // step 1 update the kinf, kh, etc. for elements.
            f_InfilRecharge(i, t); // step 2 calculate the infiltration and recharge.
            /*========surf/gw flow Function==============*/
            f_lateralFlux(i, t); // AFTER infiltration, do the lateral flux. ESP for overland flow.
        } //end of for loop.
#pragma omp for
        for (i = 0; i < NumSegmt; i++) {
            f_Segement_surface( RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
            f_Segement_sub( RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
        }
#pragma omp for
        for (i = 0; i < NumRiv; i++) {
            Flux_RiverDown(t, i);
        }
    }
#else
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
#endif
    
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
        }
        DY[i] = qEleNetPrep[i] - qEleInfil[i] - QeleSurfTot[i] / area;
        DY[iUS] = qEleInfil[i] - qEleRecharge[i];
        DY[iGW] = qEleRecharge[i] - QeleSubTot[i] / area;
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
        DY[iUS] /= Ele[i].Sy;
        DY[iGW] /= Ele[i].Sy;
//                DY[iSF] =0.0;  // debug only.
//                DY[iUS] =0.0;
//                DY[iGW] =0.0;
#ifdef _DEBUG
        CheckNANi(DY[i], i, "DY[i] (Model_Data::f_applyDY)");
        CheckNANi(DY[iUS], i, "DY[iUS] (Model_Data::f_applyDY)");
        CheckNANi(DY[iGW], i, "DY[iGW] (Model_Data::f_applyDY)");
#endif
    }
    for (int i = 0; i < NumRiv; i++) {
        DY[iRIV] = (- QrivUp[i] - QrivSurf[i] - QrivSub[i] - QrivDown[i]) / Riv[i].u_TopArea;
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
