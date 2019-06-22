#include "Model_Data.hpp"
void Model_Data::f_Segement_surface(int iEle, int iRiv, int i){
    //Surface Flux from River Segment to Element;
    double Q = 0.;
    Q = WeirFlow(Ele[iEle].zmax, uYsf[iEle],
                 Ele[iEle].zmax - Riv[iRiv].depth, uYriv[iRiv],
                 Ele[iEle].zmax, RivSeg[i].Cwr, RivSeg[i].length, Ele[iEle].depression);
    if(Q > 0 || Q < 0){
        i=i;
    }
    QrivSurf[iRiv] += Q; // Positive from River to Element
    Qe2r_Surf[iEle] += -Q; // Positive from Element to River
#ifdef _DEBUG
    CheckNANi(Q, i, "River Flux Surface (Functopm:f_Segement_surface)");
#endif
}
void Model_Data::f_Segement_sub( int iEle, int iRiv, int i){
    //Subsurface Flux from River Segment to Element;
    double  Q = 0.;
    Q = flux_R2E_GW(uYriv[iRiv], Ele[iEle].zmax - Riv[iRiv].depth,
                    uYgw[iEle], Ele[iEle].zmin,
                    Ele[iEle].u_effKH, Riv[iRiv].KsatH,
                    RivSeg[i].length,Riv[iRiv].BedThick);
    QrivSub[iRiv] += Q;
    Qe2r_Sub[iEle] += -Q;
#ifdef _DEBUG
    CheckNANi(Q, i, "River Flux Sub(Functopm:f_Segement_sub)");
#endif
}
