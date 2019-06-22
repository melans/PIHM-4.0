//
//  MD_f.cpp
//  PIHM++
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
#pragma omp for schedule(static)
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
    } //end of for loop.
    for (i = 0; i < NumSegmt; i++) {
        f_Segement_surface(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
        f_Segement_sub(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
    }
    for (i = 0; i < NumRiv; i++) {
        Flux_RiverDown(t, i);
    }
//    for (i = 0; i < NumEle; i++) {
//        f_waterbalance(i);
//    }
#endif
}
void Model_Data::f_applyDY(double *DY, double t){
    double area;
    for (int i = 0; i < NumEle; i++) {
        area = Ele[i].area;
        QeleSurfTot[i] = 0.;
        QeleSubTot[i] = 0.;
        for (int j = 0; j < 3; j++) {
            QeleSurfTot[i] += QeleSurf[i][j];
            QeleSubTot[i] += QeleSub[i][j];
        }
        QeleSurfTot[i] += Qe2r_Surf[i];
        QeleSubTot[i] += Qe2r_Sub[i];
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
        DY[iGW] /= Ele[i].Sy ;
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
//    FILE *fid = fopen("DY_debug.dat", "ab+");
//    fwrite(&t, sizeof(t), 1, fid);
//    fwrite(DY, sizeof(DY), NumY, fid);
//    fclose(fid);
//    printf("");
}
