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
#pragma omp for
        for (i = 0; i < NumEle; i++) {
            f_waterbalance(i);
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
void Model_Data::f_wb_updateQsurf(int i, int j){
    int inab;
    inab = Ele[i].nabr[j] - 1;
    for(int k = 0; k < 3; k++){
        if ( Ele[inab].nabr[k] == i + 1){
            QeleSurf[inab][k] = - QeleSurf[i][j];
        }
    }
}
void Model_Data::f_waterbalance(int i){
    double r = 1;
    QoutSurf[i] = 0.;
    if(qEleInfil[i] > 0.){
        QoutSurf[i] += qEleInfil[i] * Ele[i].area;
    }
    for(int j = 0; j < 3; j++){
        if (QeleSurf[i][j] > 0.){
            QoutSurf[i] += QeleSurf[i][j];
        }
    }
    if (qEleInfil[i] > 0.){
        if( uYsf[i] + qEleNetPrep[i] < qEleInfil[i] ){
            qEleInfil[i] = uYsf[i] + qEleNetPrep[i];
        }
    }
    if (QoutSurf[i] > 0.){
        r = (uYsf[i] + qEleNetPrep[i] - qEleInfil[i])* Ele[i].area / QoutSurf[i];
        if(r < 0.){
            r = r;
        }else if (r < 1. && r >= 0.){
            //        qEleInfil[i] *= r;
            for(int j = 0; j < 3; j++){
                if (QeleSurf[i][j] > 0.){
                    QeleSurf[i][j] *= r;
                    f_wb_updateQsurf(i, j);
                }
            }
        }
    }
    //    Ele[i].wf += qEleInfil[i] - (qEleRecharge[i] > 0. ? qEleRecharge[i] : 0.0) - qEleET[i][2];
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
        
        DY[iUS] /= Ele[i].Sy * UNIT_C;
        DY[iGW] /= Ele[i].Sy * UNIT_C;
        DY[i] /= UNIT_C;
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
        DY[iRIV] = (- QrivUp[i] - QrivSurf[i] - QrivSub[i] - QrivDown[i]) / Riv[i].u_TopArea / UNIT_C;
//        DY[iRIV] = 0.0;
//        if(DY[iRIV] + uYriv[i] < 0){
//            printf("debug: riv%d\t %f\t %f \t qup=%f\t qsf=%f\t qsub=%f\t qdown= %f\n",
//                   i+1, uYriv[i], Riv[i].depth,
//                   - QrivUp[i]/ Riv[i].u_TopArea,
//                   - QrivSurf[i]/ Riv[i].u_TopArea,
//                   - QrivSub[i]/ Riv[i].u_TopArea,
//                   - QrivDown[i] / Riv[i].u_TopArea);
//            i=i;
//        }
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
