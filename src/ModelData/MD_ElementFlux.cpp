#include "Model_Data.hpp"
void Model_Data::fun_Ele_Recharge(int i, double t){
    qEleRecharge[i] = Ele[i].Flux_Recharge(uYus[i] , uYgw[i]);
//    qEleRecharge[i] = Ele[i].u_qr;
}

void Model_Data::fun_Ele_Infiltraion(int i, double t){
    Ele[i].Flux_Infiltration(uYsf[i] , uYus[i] , uYgw[i], qEleNetPrep[i]  );
    qEleInfil[i] = Ele[i].u_qi;
    qEleExfil[i] = Ele[i].u_qex;
}
void Model_Data::fun_Ele_surface(int i, double t){
    int j, inabr;
    double  Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, CrossA;
    double isf, nsf; // Available Y in Surface of this/nabor element    
    isf = uYsf[i] - qEleInfil[i] + qEleExfil[i];
    isf = isf < 0. ? 0. : isf;
    for (j = 0; j < 3; j++) {
        inabr = Ele[i].nabr[j] - 1;
        if (inabr >= 0) {
            /***************************************************************************/
            /* Surface Lateral Flux Calculation between Triangular elements Follows */
            /***************************************************************************/
            nsf = uYsf[inabr] - qEleInfil[inabr] + qEleExfil[inabr];
            nsf = nsf < 0. ? 0. : nsf;
            Dif_Y_Surf = (isf + Ele[i].zmax) - (nsf + Ele[inabr].zmax);
            Avg_Y_Surf = avgY_sf(Ele[i].zmax, isf, Ele[inabr].zmax, nsf, Ele[i].depression);
            if(Avg_Y_Surf < 0.){
                QeleSurf[i][j] = 0.;
            }else{
                Grad_Y_Surf = Dif_Y_Surf / Ele[i].Dist2Nabor[j];
                CrossA = Avg_Y_Surf * Ele[i].edge[j];
                QeleSurf[i][j] = ManningEquation(CrossA, Ele[i].avgRough[j], Avg_Y_Surf, Grad_Y_Surf);
//                CheckNANi(QeleSurf[i][j], i, "QeleSurf[i][j]");
            } //end of ifelse Avg_Y_Surf < 0.
        } else {
            QeleSurf[i][j] = 0;
        } // end of if
    } // end of for loop
}// end of functions


void Model_Data::fun_Ele_sub(int i, double t){
    int j, inabr;
    double  Avg_Y_Sub, dy_sub, Avg_Ksat, Grad_Y_Sub, effK, effKnabr;
    
    for (j = 0; j < 3; j++) {
        inabr = Ele[i].nabr[j] - 1;
        if (inabr >= 0) {
            /***************************************************************************/
            /* Subsurface Lateral Flux Calculation between Triangular elements Follows */
            /***************************************************************************/
            dy_sub = (uYgw[i] + Ele[i].zmin) - (uYgw[inabr] + Ele[inabr].zmin);
            if(dy_sub > 0. && uYgw[i] <=EPSILON){
                dy_sub = 0;
                QeleSub[i][j] = 0.;
            }else if(dy_sub < 0. && uYgw[inabr]<= EPSILON){
                dy_sub = 0.;
                QeleSub[i][j] = 0.;
            }else{
                Avg_Y_Sub = avgY_gw(Ele[i].zmin, uYgw[i], Ele[inabr].zmin, uYgw[inabr], 0.002);  
                Grad_Y_Sub = dy_sub / Ele[i].Dist2Nabor[j];
                /* Take care of macropore effect */
                effK = effKH( uYgw[i], Ele[i].AquiferDepth, Ele[i].macD, Ele[i].macKsatH, Ele[i].geo_vAreaF, Ele[i].KsatH);
                effKnabr = effKH(uYgw[inabr], Ele[inabr].AquiferDepth, Ele[inabr].macD, Ele[inabr].macKsatH, Ele[inabr].geo_vAreaF, Ele[inabr].KsatH);
                /* It should be weighted average. However, there is an ambiguity about distance used */
                Avg_Ksat = 0.5 * (effK + effKnabr);
                QeleSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * Ele[i].edge[j];
            }
        } else {
            QeleSub[i][j] = 0;
        } // end of if
    } // end of for loop
}// end of functions
