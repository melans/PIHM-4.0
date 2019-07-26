#include "Model_Data.hpp"
void Model_Data::f_update(double  *Y, double *DY, double t){
    /* Initialization of temporary state variables */
    for (int i = 0; i < NumY; i++) {
        DY[i] = 0.;
    }
    for (int i = 0; i < NumEle; i++) {
        //        uYsf[i] = (Y[iSF] >= 0.) ? Y[iSF] : 0.;
        //        uYus[i] = (Y[iUS] >= 0.) ? Y[iUS] : 0.;
        //        uYgw[i] = (Y[iGW] >= 0.) ? Y[iGW] : 0.;
        // code above will cause balance issue, even though it is not very significant mostly.
        uYsf[i] = Y[iSF];
        uYus[i] = Y[iUS];
        uYgw[i] = Y[iGW];
    }
    for (int i = 0; i < NumEle; i++) {
        for(int j = 0; j<3;j++){
            QeleSub[i][j] = 0.;
            QeleSurf[i][j] = 0.;
            Ele[i].iupdSF[j] = 0;
            Ele[i].iupdGW[j] = 0;
        }
        Qe2r_Surf[i] = 0.;
        Qe2r_Sub[i] = 0.;
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
    }//end of for j=1:3
    
    for (int i = 0; i < NumRiv; i++ ){
        uYriv[i] = (Y[iRIV] >= 0.) ? Y[iRIV] : 0.;
        /* qrivsurf and qrivsub are calculated in Element fluxes.
         qrivDown and qrivUp are calculated in River fluxes. */
        QrivSurf[i] = 0.;
        QrivSub[i] = 0.;
        QrivUp[i] = 0.;
        QrivDown[i] = 0.;
        Riv[i].updateRiver(uYriv[i]);
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

