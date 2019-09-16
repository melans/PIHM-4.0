#include "ModelConfigure.hpp"
#include "Macros.hpp"
#include "functions.hpp"
#include "Model_Data.hpp"
void Model_Data::Flux_RiverDown(double t, int i){
    double  Distance, CSarea, Perem, R, s, n;
    double  sMean = 0.;
    int id = Riv[i].down - 1;
    n = Riv[i].avgRough;
    if (id >= 0) {
        sMean = (Riv[i].BedSlope + Riv[id].BedSlope) * 0.5 ;
        Distance =  Riv[i].Dist2DownStream;
        s = (uYriv[i] - uYriv[id]) / Distance + sMean;
//        s = sMean;
//        CSarea = 0.5 * (Riv[i].u_CSarea + Riv[iDownStrm].u_CSarea);
//        Perem = 0.5 * (Riv[i].u_CSperem + Riv[iDownStrm].u_CSperem);
        CSarea = Riv[i].u_CSarea;  // if avg, this occurs 0.5 * (10w * 0.1h + 100w * 0.1) = 55wh;
        Perem = Riv[i].u_CSperem;
        R = (Perem <= 0.) ? 0. : (CSarea / Perem);
        if(s > 0.){
            QrivDown[i] = ManningEquation(CSarea, n, R, s);
        }else{
            QrivDown[i] = ManningEquation(CSarea, n, R, s);
        }
    } else {
        switch (Riv[i].down) {
            case -1:
                /* Newman Condition, not ready yet*/
                //                break;
            case -2:
                //                break;
            case -3:
                /* zero-depth-gradient boundary conditions */
                Perem = Riv[i].u_CSperem;
                s = Riv[i].BedSlope + uYriv[i] * 2. / Riv[i].Length;
//                s = Riv[i].BedSlope;
                CSarea = Riv[i].u_CSarea;
                R = (Perem <= 0.) ? 0. : (CSarea / Perem);
                QrivDown[i] = ManningEquation(CSarea, n, R, s);
                break;
            case -4:
                /* Critical Depth boundary conditions */
                QrivDown[i] = Riv[i].u_CSarea * sqrt(GRAV * uYriv[i]) * 60.;    /* Note the dependence on physical units */
                break;
            default:
                printf("Fatal Error: River Routing Boundary Condition Type Is Wrong!");
                exit(1);
        }//end of switch
    } //end of if
#ifdef _DEBUG
    CheckNANi(QrivDown[i], i, "RiverFlux Down");
#endif
}

double Model_Data::WeirFlow(double ze, double ye,
                            double zr, double yr,
                            double zbank,
                            double cwr, double rivLen,
                            double threshold){
    double he, hr, Q=0.;
    double dh, y;
    he = ye + ze;
    hr = yr + zr;
    dh = hr - he;
    if(dh > 0.){ // from River to Element. Q is Positive.
        y = dh;
        Q = 2. / 3. * cwr * sqrt(2. * GRAV * y) * rivLen * y* 60.; /* 60 is for m3/s  to m3/min, because GRAV is m/s2*/
    }else{ // from Element to River. Q is Negative.
        if( ye > threshold){
            if( hr > ze){
                y = -dh;
            }else{
                y = ye;
            }
            Q = -2./ 3. * cwr * sqrt(2. * GRAV * y) * rivLen * y* 60.;
        }else{
            Q = 0.;
        }
    }
    return Q;
}

void Model_Data::fun_Seg_surface(int iEle, int iRiv, int i){
    //Surface Flux from River Segment to Element;
    double isf  = uYsf[iEle] - qEleInfil[iEle] + qEleExfil[iEle];
    isf = max(0., isf);
    QsegSurf[i] = WeirFlow(Ele[iEle].zmax, isf,
                           Ele[iEle].zmax - Riv[iRiv].depth, uYriv[iRiv],
                           Ele[iEle].zmax, RivSeg[i].Cwr, RivSeg[i].length, Ele[iEle].depression);
    
#ifdef _DEBUG
    CheckNANi(QsegSurf[i], i, "River Flux Surface (Functopm:f_Segement_surface)");
#endif
}
void Model_Data::fun_Seg_sub( int iEle, int iRiv, int i){
    //Subsurface Flux from River Segment to Element;
    QsegSub[i] = flux_R2E_GW(uYriv[iRiv], Ele[iEle].zmax - Riv[iRiv].depth,
                             uYgw[iEle], Ele[iEle].zmin,
                             Ele[iEle].u_effKH, Riv[iRiv].KsatH,  
                             RivSeg[i].length,Riv[iRiv].BedThick);
#ifdef _DEBUG
    CheckNANi(Q, i, "River Flux Sub(Functopm:fun_Seg_sub)");
#endif
}
