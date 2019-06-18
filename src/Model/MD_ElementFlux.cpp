#include "Model_Data.hpp"
void Model_Data::f_InfilRecharge(int i, double t){
    Ele[i].Flux_InfiRech(uYsf[i] , uYus[i] , uYgw[i], qEleNetPrep[i]  );
    qEleInfil[i] = Ele[i].u_qi;
    if(qEleInfil[i] > 0. && qEleInfil[i] > uYsf[i]){
        if(uYsf[i] > 0){
            qEleInfil[i] = uYsf[i];
        }else{
            qEleInfil[i] = 0.;
        }
    }
    qEleRecharge[i] = Ele[i].u_qr;
    if(uYsf[i] < EPSILON){ // those should out of this function. Next.
        qEleET[i][2] = Ele[i].u_satn * qEleET[i][2];
        if(Ele[i].u_satn <= EPS){
            qEleET[i][2] = 0.;
        }
    }
}


void Model_Data::f_lateralFlux(int i, double t){
    int j, inabr;
    double  Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, Avg_Sf, CrossA;
    double  Avg_Y_Sub, dy_sub, Avg_Ksat, Grad_Y_Sub, effK, effKnabr;
    double isf, nsf; // Available Y in Surface of this/nabor element
    isf = uYsf[i] - qEleInfil[i];
    isf = isf < 0. ? 0. : isf;
//    if(isf < 0. && uYsf[i] > 0.){
//        isf = 0.;
//    }
    for (j = 0; j < 3; j++) {
        inabr = Ele[i].nabr[j] - 1;
        if (inabr >= 0) {
            nsf = uYsf[inabr] - qEleInfil[inabr];
            nsf = nsf < 0. ? 0. : nsf;
            /***************************************************************************/
            /* Subsurface Lateral Flux Calculation between Triangular elements Follows */
            /***************************************************************************/
            dy_sub = (uYgw[i] + Ele[i].zmin) - (uYgw[inabr] + Ele[inabr].zmin);
            Avg_Y_Sub = avgY(uYgw[i], Ele[i].zmin, uYgw[inabr], Ele[inabr].zmin, 0.002);
            Grad_Y_Sub = dy_sub / Ele[i].Dist2Nabor[j];
            /* take care of macropore effect */
            effK = effKH( uYgw[i], Ele[i].AquiferDepth, Ele[i].macD, Ele[i].macKsatH, Ele[i].geo_vAreaF, Ele[i].KsatH);
            effKnabr = effKH(uYgw[inabr], Ele[inabr].AquiferDepth, Ele[inabr].macD, Ele[inabr].macKsatH, Ele[inabr].geo_vAreaF, Ele[inabr].KsatH);
            /* It should be weighted average. However, there is an ambiguity about distance used */
            Avg_Ksat = 0.5 * (effK + effKnabr);
            QeleSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * Ele[i].edge[j];
            /***************************************************************************/
            /* Surface Lateral Flux Calculation between Triangular elements Follows */
            /***************************************************************************/
            Dif_Y_Surf = (isf + Ele[i].zmax) - (nsf + Ele[inabr].zmax);
            Avg_Y_Surf = avgY(Ele[i].zmax, isf, Ele[inabr].zmax, nsf, Ele[i].depression);
            if(Avg_Y_Surf < 0){
                QeleSurf[i][j] = 0.;
            }else{
                Grad_Y_Surf = Dif_Y_Surf / Ele[i].Dist2Nabor[j];
                Avg_Sf = sqpow2(Ele[i].dhBYdx, Ele[i].dhBYdy);
                if(Avg_Sf > EPS_SLOPE){
                    /*void*/
                }else{
                    Avg_Sf = EPS_SLOPE;
                }
                CrossA = Avg_Y_Surf * Ele[i].edge[j];
                QeleSurf[i][j] = OverlandManning(Avg_Y_Surf, Grad_Y_Surf, Avg_Sf, CrossA, Ele[i].avgRough[j]);
                CheckNANij(QeleSurf[i][j], i, "QeleSurf[i][j] (Model_Data::f_applyDY)");
            } //end of ifelse Avg_Y_Surf < EPSilon
        } else {
            QeleSurf[i][j] = 0;
            QeleSub[i][j] = 0;
        } // end of if
    } // end of for loop
}// end of functions
