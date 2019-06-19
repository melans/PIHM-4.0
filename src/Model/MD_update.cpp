#include "Model_Data.hpp"
void Model_Data::f_update(double  *Y, double *DY, double t){
    //    for (int i = 0; i<NumForc; i++){
    //        tsd_weather[i].movePointer(t);
    //    }
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
        }
//        qEleET[i][1] = 0.; // cannot put them zero, otherwise water balance issue.
//        qEleET[i][2] = 0.;
        QeleSurfTot[i] = 0.;
        QeleSubTot[i] = 0.;
        Qe2r_Surf[i] = 0.;
        Qe2r_Sub[i] = 0.;
        for (int j = 0; j < 3; j++) {
            if(Ele[i].nabr[j] > 0){
                Ele[i].surfH[j] = (Ele[Ele[i].nabr[j] - 1].zmax + uYsf[Ele[i].nabr[j] - 1]);
            }else{
                Ele[i].surfH[j] = (Ele[i].zmax + uYsf[i]);
            }
            
        }
        Ele[i].dhBYdx = -1 * (Ele[i].surfY[2] * (Ele[i].surfH[1] - Ele[i].surfH[0]) +
                              Ele[i].surfY[1] * (Ele[i].surfH[0] - Ele[i].surfH[2]) +
                              Ele[i].surfY[0] * (Ele[i].surfH[2] - Ele[i].surfH[1])) /
        (Ele[i].surfX[2] * (Ele[i].surfY[1] - Ele[i].surfY[0]) +
         Ele[i].surfX[1] * (Ele[i].surfY[0] - Ele[i].surfY[2]) +
         Ele[i].surfX[0] * (Ele[i].surfY[2] - Ele[i].surfY[1]));
        Ele[i].dhBYdy = -1 * (Ele[i].surfX[2] * (Ele[i].surfH[1] - Ele[i].surfH[0]) +
                              Ele[i].surfX[1] * (Ele[i].surfH[0] - Ele[i].surfH[2]) +
                              Ele[i].surfX[0] * (Ele[i].surfH[2] - Ele[i].surfH[1])) /
        (Ele[i].surfY[2] * (Ele[i].surfX[1] - Ele[i].surfX[0]) +
         Ele[i].surfY[1] * (Ele[i].surfX[0] - Ele[i].surfX[2]) +
         Ele[i].surfY[0] * (Ele[i].surfX[2] - Ele[i].surfX[1]));
    }//end of for i=1:3
    
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
}

void Model_Data::updateWF(double dt){
    for(int i = 0; i < NumEle; i++){
        if( dt > 0.){
            Ele[i].updateWF(max(0.0, qEleInfil[i]), dt);
        }
        yEleWetFront[i] = Ele[i].u_wf;
    }
}

