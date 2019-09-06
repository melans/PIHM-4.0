//
//  Element.cpp
//  PIHM++ (v 4.0)
//
//  Created by Lele Shu on 7/17/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#include "Element.hpp"
void Triangle::printHeader(FILE *fp){
    for(int i = 0; i < 3; i++){
        fprintf(fp, "node%d\t", i);
    }
    for(int i = 0; i < 3; i++){
        fprintf(fp, "nabr%d\t", i);
    }
    for(int i = 0; i < 3; i++){
        fprintf(fp, "edge%d\t", i);
    }
    fprintf(fp, "%s\t", "area");
    fprintf(fp, "%s\t", "x");
    fprintf(fp, "%s\t", "y");
    fprintf(fp, "%s\t", "zmin");
    fprintf(fp, "%s\t", "zmax");
}
void Triangle::printInfo(FILE *fp){
    for(int i = 0; i < 3; i++){
        fprintf(fp, "%d\t", node[i]);
    }
    for(int i = 0; i < 3; i++){
        fprintf(fp, "%d\t", nabr[i]);
    }
    for(int i = 0; i < 3; i++){
        fprintf(fp, "%g\t", edge[i]);
    }
    fprintf(fp, "%g\t", area);
    fprintf(fp, "%g\t", x);
    fprintf(fp, "%g\t", y);
    fprintf(fp, "%g\t", zmin);
    fprintf(fp, "%g\t", zmax);
}
void AttriuteIndex::printHeader(FILE *fp){
    fprintf(fp, "%s\t", "iSoil");
    fprintf(fp, "%s\t", "iGeol");
    fprintf(fp, "%s\t", "iLC");
    fprintf(fp, "%s\t", "IC");
    fprintf(fp, "%s\t", "iForc");
    fprintf(fp, "%s\t", "iMF");
    fprintf(fp, "%s\t", "iBC");
}
void AttriuteIndex::printInfo(FILE *fp){
    fprintf(fp, "%d\t", iSoil);
    fprintf(fp, "%d\t", iGeol);
    fprintf(fp, "%d\t", iLC);
    fprintf(fp, "%d\t", IC);
    fprintf(fp, "%d\t", iForc);
    fprintf(fp, "%d\t", iMF);
    fprintf(fp, "%d\t", iBC);
    fprintf(fp, "%d\t", iSS);
}

void _Element::applyGeometry(_Node *Node){
    double  a_x, a_y, b_x, b_y, c_x, c_y;
    double  a_zmin, a_zmax, b_zmin, b_zmax, c_zmin, c_zmax;
    double  aqd, z0;
    a_x = Node[node[0] - 1].x;
    b_x = Node[node[1] - 1].x;
    c_x = Node[node[2] - 1].x;
    a_y = Node[node[0] - 1].y;
    b_y = Node[node[1] - 1].y;
    c_y = Node[node[2] - 1].y;
    
    a_zmin = Node[node[0] - 1].zmin;
    b_zmin = Node[node[1] - 1].zmin;
    c_zmin = Node[node[2] - 1].zmin;
    a_zmax = Node[node[0] - 1].zmax;
    b_zmax = Node[node[1] - 1].zmax;
    c_zmax = Node[node[2] - 1].zmax;
    
    area = 0.5 * ((b_x - a_x) * (c_y - a_y) - (b_y - a_y) * (c_x - a_x));
    
    zmax = (a_zmax + b_zmax + c_zmax) / 3.0;
    zmin = (a_zmin + b_zmin + c_zmin) / 3.0;
    if(zcentroid != NA_VALUE){
        z0 = zmax;
        aqd = zmax - zmin;
//        zmax = (zmax+ zcentroid) * 0.5;
        zmax = (a_zmax + b_zmax + c_zmax + zcentroid) / 4.0;
//        zmax = (a_zmax + b_zmax + c_zmax) / 3.0;  // debug
        zmin = zmax - aqd;
        if(fabs(z0 - zmax) > 10.){
//            printf("DZ(%d) = %f\n", index, z0 - zmax);
        }
    }
    /* calculate centroid of triangle */
    x = (a_x + b_x + c_x) / 3.0;
    y = (a_y + b_y + c_y) / 3.0;
    edge[0] = Eudist(b_x, b_y, c_x, c_y);
    edge[1] = Eudist(c_x, c_y, a_x, a_y);
    edge[2] = Eudist(a_x, a_y, b_x, b_y);
}
void _Element::InitElement(){
    AquiferDepth = zmax - zmin;
    WetlandLevel = AquiferDepth - infD;
    RootReachLevel = AquiferDepth - RzD;
    FixPressure = PressureElevation(zmax); // P atmospheric pressure [kPa]
//    FixGamma = PsychrometricConstant(FixPressure);
    MacporeLevel = AquiferDepth - macD;
    
    if (AquiferDepth < macD)
        macD = AquiferDepth;
    
    CheckNonZero(RootReachLevel, index-1, "RootReachLevel");
    CheckNonZero(FixPressure, index-1, "FixPressure");
    CheckNonZero(area, index-1, "area");
    CheckNA(zmax, "zmax");
    CheckNA(zmin, "zmin");
    for(int i = 0; i < 3; i++){
        CheckNonZero(edge[i], index-1, "edge");
    }
    
}
void _Element::applyNabor(_Node *Node, _Element *Ele){
    int eNabor;
    for (int j = 0; j < 3; j++) {
        if(nabr[j]>0){
            for(int k = 0; k < 3; k++){
                if(Ele[nabr[j] - 1].nabr[k] == this->index){
                    this->nabrToMe[j] = k;
                    /* Neighbor's K direction and My J direction shared a edge*/
                }
            }
        }else{
            /* VOID */
        }
    }

    for(int j = 0; j < 3; j++){
        eNabor = nabr[j] - 1;
        Dist2Edge[j] = sqpow2(edge[0] * edge[1] * edge[2] / (4 * area), edge[j] / 2 );
        //          Dist2Edge[j] =  sqrt(pow(edge[0] * edge[1] * edge[2] / (4 * area), 2) - pow(edge[j] / 2, 2));
        if(eNabor >= 0){
            Dist2Nabor[j] = Eudist(x, y,
                                   Ele[eNabor].x, Ele[eNabor].y);
            avgRough[j] = 0.5 * (Rough + Ele[eNabor].Rough);
        }else{
            Dist2Nabor[j] = 0.;
            avgRough[j] = Rough;
        }
        //            Dist2Nabor[j] = sqrt(pow((x - Ele[eNabor].x), 2) + pow((y - Ele[eNabor].y), 2));
    }
    
}
void _Element::Flux_InfiRech(double Ysurf, double Yunsat, double Ygw, double netprcp){
    double ke=0.;
    double grad;
    if(Ygw > AquiferDepth){
        /* GW reach the surface */
        u_qex = (Ygw - AquiferDepth) / 1. * Kmax;
        u_qi = 0. ;
    }else{
        u_qex = 0.;
        Ysurf = max(Ysurf, 0.);
        u_effkInfi = infKsatV * (1 - hAreaF) + hAreaF * macKsatV * u_satn;
        /* NORMAL GW level */
        u_qi = (Ysurf + infD) / infD * u_effkInfi;
        if(u_qi > 0.){
            if(Ysurf > 0.){
                u_qi = min(Ysurf, u_qi);
            }else{
                u_qi = 0.;
            }
        }else{
            /* void */
        }
    }
    /* Recharge */
//    grad = (1. + u_phius * 2. / AquiferDepth);
    grad = Yunsat / AquiferDepth;
    if(grad >= 0.){
//        ke = meanHarmonic(infKsatV * u_satKr, KsatV, u_deficit, Ygw);
        ke = meanArithmetic(infKsatV * u_satKr, KsatV, u_deficit, Ygw);
        u_qr = grad * ke;
        if(u_deficit > 0.){
            u_qr *= Yunsat / u_deficit;
        }
        u_qr = min(Yunsat, u_qr);
    }else{
        u_qr = 0.;
//        u_qr = max(-1.*Ygw, u_qr);
    }
    CheckNANi(u_qi, 0, "u_qi");
}
void _Element::updateElement(double Ysurf, double Yunsat, double Ygw){
    u_effKH = effKH(Ygw,  AquiferDepth,  macD,  macKsatH,  geo_vAreaF,  KsatH);
    u_deficit = AquiferDepth - Ygw;
    Kmax = infKsatV * (1. - hAreaF) + macKsatV * hAreaF ;
    u_qex = 0.;
    u_qi = 0.;
    u_qr = 0.;
    if(u_deficit <= 0. ){
        u_deficit = 0.;
        u_satn = 1.;
    }else{
        u_satn = (Yunsat / u_deficit * ThetaS - ThetaR) / (ThetaS - ThetaR) ;
    }
    if(u_satn > 0.99 ){
        u_satn = 1.0;
        u_satKr = 1.0;
        u_phius = 0.;
    }else if(u_satn <= EPS_DOUBLE){
        u_satn = 0.;
        u_satKr = 0.;
        u_phius = MINPSI;
    }else{
        u_satKr = satKfun(u_satn, Beta);
        u_phius = sat2psi(u_satn, Alpha, Beta);
        u_phius = max(MINPSI, u_phius);
    }
    u_effkInfi = infKsatV * (1 - hAreaF) + u_satn * macKsatV * hAreaF ;
    u_effkInfi *=  1. - u_satn;
#ifdef _DEBUG
    if (u_effkInfi < 0.){
        printf("WARNING: Negative effective conductivity for infiltration.\n");
    }
    if (u_effKH < 0.000000001){
        printf("WARNING: Negative effective conductivity for infiltration.\n");
    }
#endif
}

void _Element::copyGeol(Geol_Layer *g){
    KsatH        = g[iGeol - 1].KsatH;
    KsatV        = g[iGeol - 1].KsatV;
    geo_ThetaS    = g[iGeol - 1].geo_ThetaS;
    geo_ThetaR    = g[iGeol - 1].geo_ThetaR;
    geo_vAreaF    = g[iGeol - 1].geo_vAreaF;
    macKsatH    = g[iGeol - 1].macKsatH;
    macD        = g[iGeol - 1].macD;
    Sy    = g[iGeol - 1].Sy;
}
void _Element::copySoil(Soil_Layer *g){
    infKsatV = g[iSoil - 1].infKsatV;
    ThetaS = g[iSoil - 1].ThetaS;
    ThetaR = g[iSoil - 1].ThetaR;
    Alpha = g[iSoil - 1].Alpha;
    Beta = g[iSoil - 1].Beta;
    hAreaF = g[iSoil - 1].hAreaF;
    macKsatV = g[iSoil - 1].macKsatV;
    infD = g[iSoil - 1].infD;
    
    CheckNonZero(ThetaS, index-1, "ThetaS");
    CheckNonZero(ThetaR, index-1, "ThetaR");
    CheckNonZero(infD, index-1, "infD");
}
void _Element::copyLandc(Landcover *g){
//    LAImax   = g[iLC - 1].LAImax;
    VegFrac  = g[iLC - 1].VegFrac;
    Albedo   = g[iLC - 1].Albedo;
    Rs_ref   = g[iLC - 1].Rs_ref;
    Rmin     = g[iLC - 1].Rmin;
    Rough    = g[iLC - 1].Rough;
    RzD      = g[iLC - 1].RzD;
    SoilDgrd = g[iLC - 1].SoilDgrd;
    ImpAF    = g[iLC - 1].ImpAF;
}
//void _Element::updateWF(double qi, double dt){
//    double dh = 0.;
//    double grad;
//    double effk;
//    if( qi > 0. ){
//        /* Infiltration occur
//         wf increases, bottom loss
//         */
//        dh += qi;
//        if( u_wf > 0.){
//            /* available water in wf */
//            effk = u_satKr * infKsatV;
//            if(u_deficit - u_wf <= 0){
//                grad = 0;
//            }else{
//                grad = (u_deficit * 0.5 - u_phius) / ( 0.5 * (u_deficit - u_wf) );
//            }
//            dh += -effk * grad;
//        }
//        u_wf +=  dh * dt;
//        if(u_wf > 0.01){
//            dh = dh;
//        }
//    }else{
//        /* No infiltration: sat layer moves forward */
//        u_wf = infD;
//    }
//    CheckNA(u_wf, "Weting Front");
//}
void _Element::printHeader(FILE *fp){
    fprintf(fp, "%s\t", "index");
    AttriuteIndex::printHeader(fp);
    Triangle::printHeader(fp);
    Soil_Layer::printHeader(fp);
    Geol_Layer::printHeader(fp);
    Landcover::printHeader(fp);
    for(int i = 0; i < 3; i++){
        fprintf(fp, "Dist2Edge%d\t", i);
    }
    for(int i = 0; i < 3; i++){
        fprintf(fp, "Dist2Nabor%d\t", i);
    }
    for(int i = 0; i < 3; i++){
        fprintf(fp, "avgRough%d\t", i);
    }
    fprintf(fp, "%s\t", "windH");
    fprintf(fp, "%s\t", "FixPressure");
    fprintf(fp, "%s\t", "AquiferDepth");
    fprintf(fp, "%s\t", "WetlandLevel");
    fprintf(fp, "%s\t", "RootReachLevel");
    fprintf(fp, "%s\t", "MacporeLevel");
    fprintf(fp, "\n");
}
void _Element::printInfo(FILE *fp){
    fprintf(fp, "%d\t", index);
    AttriuteIndex::printInfo(fp);
    Triangle::printInfo(fp);
    Soil_Layer::printInfo(fp);
    Geol_Layer::printInfo(fp);
    Landcover::printInfo(fp);
    for(int i = 0; i < 3; i++){
        fprintf(fp, "%g\t", Dist2Edge[i]);
    }
    for(int i = 0; i < 3; i++){
        fprintf(fp, "%g\t", Dist2Nabor[i]);
    }
    for(int i = 0; i < 3; i++){
        fprintf(fp, "%g\t", avgRough[i]);
    }
    fprintf(fp, "%g\t", windH);
    fprintf(fp, "%g\t", FixPressure);
    fprintf(fp, "%g\t", AquiferDepth);
    fprintf(fp, "%g\t", WetlandLevel);
    fprintf(fp, "%g\t", RootReachLevel);
    fprintf(fp, "%g\t", MacporeLevel);
    fprintf(fp, "\n");
}
