#include "Model_Data.hpp"
void Model_Data::f_channel_surface(int iEle, int iRiv, int i){
    //Surface Flux from River Channel to Element;
    double Q = 0.;
    double zbed, zbank;
    zbed = Ele[iEle].zmax - Riv[iRiv].depth;
    zbank = Ele[iEle].zmax;
//    Q = WeirFlow(Ele[iEle].zmax, uYsf[iEle],
//                 Ele[iEle].zmax - Riv[iRiv].depth, uYriv[iRiv],
//                 Ele[iEle].zmax, channel[i].Cwr, channel[i].length, Ele[iEle].depression);
    Q = flux_R2E_SF(uYriv[iRiv], Ele[iEle].zmax - Riv[iRiv].depth,
                    uYsf[iEle], Ele[iEle].zmax, Ele[iEle].Rough,
                    channel[i].length,
                    Ele[iEle].area / channel[i].length * 0.25, Ele[iEle].depression);
//    if( Q < 0. ){
//        if( -Q > uYsf[iEle] * Ele[iEle].area * UNIT_C){
//            Q = max(Q , -1.0 * uYsf[iEle] * Ele[iEle].area * UNIT_C);
//        }
//    }else{
//        if ( uYsf[iEle] < Ele[iEle].depression){
//            Q = 0.;
//        }
//    }
    if(iEle == ID){
        i=i;
    }
    QrivSurf[iRiv] += Q; // Positive from River to Element
    Qe2r_Surf[iEle] += -Q; // Positive from Element to River
    
//    if(uYsf[iEle] > 0.1 && iEle == ID){
//        printf("%d %.2e -> %.2e \n",
//               iEle+1, uYsf[iEle], -Qe2r_Surf[iEle] / Ele[iEle].area);
//        i=i;
//    }
//    if(Qe2r_Surf[iEle] / Ele[iEle].area>1.){
//        printf("%d %.2f -> %.2f %.2f %.2f\n",
//               iEle+1, uYsf[iEle], qEleNetPrep[iEle], qEleInfil[iEle], Q / Ele[iEle].area);
//        i=i;
//    }
    
#ifdef _DEBUG
    CheckNANi(Q, i, "River Flux Surface (Functopm:f_channel_surface)");
#endif
}

void Model_Data::f_channel_sub( int iEle, int iRiv, int i){
    //Subsurface Flux from River Channel to Element;
    double  Q = 0.;
    Q = flux_R2E_GW(uYriv[iRiv], Ele[iEle].zmax - Riv[iRiv].depth,
                    uYgw[iEle], Ele[iEle].zmin,
                    Ele[iEle].u_effKH, Riv[iRiv].KsatH,
                    channel[i].length,Riv[iRiv].BedThick);
    QrivSub[iRiv] += Q;
    Qe2r_Sub[iEle] += -Q;
#ifdef _DEBUG
    CheckNANi(Q, i, "River Flux Sub(Functopm:f_channel_sub)");
#endif
}
