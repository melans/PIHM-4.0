//
//  Element.cpp
//  PIHM++
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
    FixGamma = PsychrometricConstant(FixPressure);
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
    double        a_x, a_y, b_x, b_y, c_x, c_y, distX, distY;
    int eNabor ;
    a_x = Node[node[0] - 1].x;
    b_x = Node[node[1] - 1].x;
    c_x = Node[node[2] - 1].x;
    a_y = Node[node[0] - 1].y;
    b_y = Node[node[1] - 1].y;
    c_y = Node[node[2] - 1].y;
    for (int j = 0; j < 3; j++) {
        /*
         * Note: Assumption here is that the forumulation is
         * circumcenter based
         */
        switch (j) {
            case 0:
                distX = (x - 0.5 * (b_x + c_x));
                distY = (y - 0.5 * (b_y + c_y));
                break;
            case 1:
                distX = (x - 0.5 * (c_x + a_x));
                distY = (y - 0.5 * (c_y + a_y));
                break;
            case 2:
                distX = (x - 0.5 * (a_x + b_x));
                distY = (y - 0.5 * (a_y + b_y));
                break;
            default:
                distX = 0.;
                distY = 0.;
                
        }
        surfH[j] = (nabr[j] > 0) ? (Ele[nabr[j] - 1].zmax) : (zmax);
        surfX[j] = (nabr[j] > 0)
        ? Ele[nabr[j] - 1].x
        : (x - 2 * distX);
        surfY[j] = nabr[j] > 0
        ? Ele[nabr[j] - 1].y
        : (y - 2 * distY);
    }
    dhBYdx = -(surfY[2] * (surfH[1] - surfH[0]) + surfY[1] * (surfH[0] - surfH[2]) + surfY[0] * (surfH[2] - surfH[1])) / (surfX[2] * (surfY[1] - surfY[0]) + surfX[1] * (surfY[0] - surfY[2]) + surfX[0] * (surfY[2] - surfY[1]));
    dhBYdy = -(surfX[2] * (surfH[1] - surfH[0]) + surfX[1] * (surfH[0] - surfH[2]) + surfX[0] * (surfH[2] - surfH[1])) / (surfY[2] * (surfX[1] - surfX[0]) + surfY[1] * (surfX[0] - surfX[2]) + surfY[0] * (surfX[2] - surfX[1]));
    
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
    double H_us =  u_phius + Yunsat * .5;
    if(u_deficit <= infD){
        /* gw reach surface -> Rechage is negative -- from gw to surface */
//        u_qr = AquiferDepth - Ygw;
//        u_qi = u_qr < 0. ? u_qr : u_qi;
        u_qr = 0.;
        u_qi = 0.;
        if(u_qi < 0.){
            u_qi = 0.;
        }
    }else{
        if(ImpAF < 1.){
            u_qi = GreenAmpt(u_effkInfi, u_Theta, ThetaS, u_phius, u_wf, Ysurf + netprcp, Sy);
        }else{
            u_qi = 0.0;
        }
        u_qr = fun_recharge(u_satKr * infKsatV, KsatV, u_deficit, Ygw, H_us);
    }
    if(u_qr > 0. && Yunsat <= EPS){
        u_qr = 0.;
    }
    if(u_qr < 0. && Ygw <= EPS){
        u_qr = 0.;
    }
}
void _Element::updateElement(double Ysurf, double Yunsat, double Ygw){
    u_effKH = effKH(Ygw,  AquiferDepth,  macD,  macKsatH,  geo_vAreaF,  KsatH);
    u_deficit = AquiferDepth - Ygw;
    if(u_deficit <= 0. ){
        u_deficit = 0.;
        u_satn = 1.;
    }else{
        u_satn = Yunsat / u_deficit;
    }
    if(u_satn > 0.99 ){
        u_satn = 1.0;
        u_satKr = 1.0;
        u_phius = 0.;
    }else{
        u_satKr = satKfun(u_satn, Beta);
        u_phius = sat2Y(u_satn, Alpha, Beta, -1 * u_deficit);
    }
    //    u_ThetaFC = FieldCapacity(Alpha, Beta, u_deficit);
//    u_Theta = u_satn * ThetaS;
    u_effkInfi = effKVnew(u_satKr, macKsatV, infKsatV, hAreaF, Ysurf);
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
    LAImax   = g[iLC - 1].LAImax;
    VegFrac  = g[iLC - 1].VegFrac;
    Albedo   = g[iLC - 1].Albedo;
    Rs_ref   = g[iLC - 1].Rs_ref;
    Rmin     = g[iLC - 1].Rmin;
    Rough    = g[iLC - 1].Rough;
    RzD      = g[iLC - 1].RzD;
    SoilDgrd = g[iLC - 1].SoilDgrd;
    ImpAF    = g[iLC - 1].ImpAF;
}
void _Element::updateWF(double qi, double dt){
    double dh = 0.;
    double grad;
    double effk;
    if( qi > 0. ){
        /* Infiltration occur
         wf increases, bottom loss
         */
        dh += qi;
        if( u_wf > 0.){
            /* available water in wf */
            effk = u_satKr * infKsatV;
            if(u_deficit - u_wf <= 0){
                grad = 0;
            }else{
                grad = (u_deficit * 0.5 - u_phius) / ( 0.5 * (u_deficit - u_wf) );
            }
            dh += -effk * grad;
        }
        u_wf +=  dh * dt / UNIT_C;
        if(u_wf > 0.01){
            dh = dh;
        }
    }else{
        /* No infiltration: sat layer moves forward */
        u_wf = infD;
    }
    CheckNA(u_wf, "Weting Front");
}
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
