#include "Model_Data.hpp"
void Model_Data::f_update(double  *Y, double *DY, double t){
    /* Initialization of temporary state variables */
    for (int i = 0; i < NumY; i++) {
        DY[i] = 0.;
    }
    for (int i = 0; i < NumEle; i++) {
        uYsf[i] = (Y[iSF] >= 0.) ? Y[iSF] : 0.;
        uYus[i] = (Y[iUS] >= 0.) ? Y[iUS] : 0.;
//        uYgw[i] = (Y[iGW] >= 0.) ? Y[iGW] : 0.;
//        code above will cause balance issue, even though it is not very significant mostly.
//        uYsf[i] = max(Y[iSF], 0.);
//        uYus[i] = Y[iUS];
//        uYgw[i] = Y[iGW];
        if(Ele[i].iBC == 0){ // NO BC
            uYgw[i] = max(0.0, Y[iGW]);
            Ele[i].QBC = 0.;
        }else if(Ele[i].iBC > 0){ // BC fix head
            Ele[i].yBC = tsd_eyBC.getX(t, Ele[i].iBC);
            uYgw[i] = Ele[i].yBC;
            Ele[i].QBC = 0.;
        }else{ // BC fix flux to GW
            Ele[i].QBC = tsd_eqBC.getX(t, -Ele[i].iBC);
        }
        /***** SS and BC *****/
        for(int j = 0; j<3;j++){
            QeleSub[i][j] = 0.;
            QeleSurf[i][j] = 0.;
            Ele[i].iupdSF[j] = 0;
            Ele[i].iupdGW[j] = 0;
        }
        Qe2r_Surf[i] = 0.;
        Qe2r_Sub[i] = 0.;
        qEleExfil[i] = 0.;
        qEleInfil[i] = 0.;
/********* Below are remove because the bass-balance issue. **********/
//        for (int j = 0; j < 3; j++) {
//            if(Ele[i].nabr[j] > 0){
//                Ele[i].surfH[j] = (Ele[Ele[i].nabr[j] - 1].zmax + uYsf[Ele[i].nabr[j] - 1]);
//            }else{
//                Ele[i].surfH[j] = (Ele[i].zmax + uYsf[i]);
//            }
//        }
//        Ele[i].dhBYdx = dhdx(Ele[i].surfX, Ele[i].surfY, Ele[i].surfH);
//        Ele[i].dhBYdy = dhdy(Ele[i].surfX, Ele[i].surfY, Ele[i].surfH);
//        Ele[i].Avg_Sf = sqpow2(Ele[i].dhBYdx, Ele[i].dhBYdy);
    }//end of for j=1:NumEle
    
    for (int i = 0; i < NumRiv; i++ ){
//        uYriv[i] = (Y[iRIV] >= 0.) ? Y[iRIV] : 0.;
        /* qrivsurf and qrivsub are calculated in Element fluxes.
         qrivDown and qrivUp are calculated in River fluxes. */
        QrivSurf[i] = 0.;
        QrivSub[i] = 0.;
        QrivUp[i] = 0.;
        QrivDown[i] = 0.;
        Riv[i].updateRiver(uYriv[i]);
        /***** SS and BC *****/
        Riv[i].qBC = 0.0;
        uYriv[i] = max(0.0,  Y[iRIV]);
        if(Riv[i].BC == 0){
            /* Void */
        }else if(Riv[i].BC < 0){ // Fixed Flux INTO river Reaches.
            Riv[i].qBC = tsd_rqBC.getX(t, -Riv[i].BC);
        }else if (Riv[i].BC > 0){ // Fixed Stage of river reach.
            Riv[i].yBC = tsd_ryBC.getX(t, Riv[i].BC);
            uYriv[i] = Riv[i].yBC;
        }
    }
    for (int i = 0; i < NumSegmt; i++ ){
        QsegSurf[i] = 0.;
        QsegSub[i] = 0.;
    }
}

//void Model_Data::updateWF(double dt){
//    for(int i = 0; i < NumEle; i++){
//        if( dt > 0.){
//            Ele[i].updateWF(max(0.0, qEleInfil[i]), dt);
//        }
//        yEleWetFront[i] = Ele[i].u_wf;
//    }
//}


void Model_Data::summary (N_Vector udata){
    double  *Y;
#ifdef _PIHMOMP
    Y = NV_DATA_OMP(udata);
#else
    Y = NV_DATA_S(udata);
#endif
    for (int i = 0; i < NumEle; i++){
        yEleSurf[i] = Y[iSF];
        yEleUnsat[i] = Y[iUS];
        
        if(Ele[i].iBC > 0){
            yEleGW[i] = Ele[i].yBC;
        }else{
            yEleGW[i] = Y[iGW];
        }
    }
    for (int i = 0; i < NumRiv; i++){
        yRivStg[i] = Y[iRIV];
        //        uYriv[i] = Y[iRIV];
        if(Riv[i].BC > 0){
            yRivStg[i] = Riv[i].yBC;
        }else{
            yRivStg[i] = Y[iRIV];
        }
    }
}
void Model_Data::summary (N_Vector udata1, N_Vector udata2){
    double  *Y1, *Y2;
#ifdef _PIHMOMP
    Y1 = NV_DATA_OMP(udata1);
    Y2 = NV_DATA_OMP(udata2);
#else
    Y1 = NV_DATA_S(udata1);
    Y2 = NV_DATA_S(udata2);
#endif
    for (int i = 0; i < NumEle; i++){
        yEleSurf[i] = Y1[iSF];
        yEleUnsat[i] = Y1[iUS];
        yEleGW[i] = Y1[iGW];
    }
    for (int i = 0; i < NumRiv; i++){
        yRivStg[i] = Y2[i];
        if(Y2[i] < 0){
            printf("Yriv[%d]=%f\n", i+1, Y2[i]);
        }
    }
}

