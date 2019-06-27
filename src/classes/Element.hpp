//
//  Element.hpp
//  PIHM++ (v 4.0)
//
//  Created by Lele Shu on 7/17/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#ifndef Element_hpp
#define Element_hpp

#include <stdio.h>
#include "Macros.hpp"
#include "Node.hpp"
#include "ModelConfigure.hpp"
#include "Equations.hpp"
#include "is_sm_et.hpp"

//class Element_Calib{
//    calib_soil csoil;
//    calib_geol cgeol;
//    calib_landcover clandc;
//};
class Triangle{
public:
    int node[3];/* anti-clock-wise */
    int nabr[3];/* neighbor i shares edge i (0: on boundary) */
    int nabrToMe[3] = {-1, -1, -1}; /* The index j of my nabor */
    double edge[3];/* edge i is from node i to node i+1 */
    double area = NA_VALUE;    /* area of element */
    
    double x = NA_VALUE;    /* x of centroid */
    double y = NA_VALUE;    /* y of centroid */
    double zmin = NA_VALUE;    /* z_min of centroid */
    double zmax = NA_VALUE;    /* z_max of centroid */
    double  zcentroid = NA_VALUE; 
    void printHeader(FILE *fp);
    void printInfo(FILE *fp);
};
/* Definition of Global Variable Types */
class AttriuteIndex{
public:
    int iSoil = NA_VALUE;   /* soil type */
    int iGeol = NA_VALUE;   /* geology type */
    int iLC = NA_VALUE;     /* Land Cover type */
    int IC  = NA_VALUE;     /* initial condition type */
    int iForc = NA_VALUE;   /* precipitation (forcing) type */
    int iMF = NA_VALUE;     /* meltFactor */
    int iBC = NA_VALUE;     /* boundary type. 0:natural bc (no flow); 1:Dirichlet BC; 2:Neumann BC */
    void printHeader(FILE *fp);
    void printInfo(FILE *fp);
};

class _Element : public Triangle,
public Soil_Layer,
public Geol_Layer,
public Landcover ,
public AttriuteIndex
{    /* Data model for a triangular element */
public:
    int index;    /* Element No. */
    int RivID = 0;
    double windH = NA_VALUE;    /* wind measurement height */
    /* for calculation of dh/ds */
//    double surfH[3];    /* Total head in neighboring cells */
//    double surfX[3];    /* Center X location of neighboring cells */
//    double surfY[3];   /* Center Y location of neighboring cells */
//    double dhBYdx = NA_VALUE;    /* Head gradient in x dirn. */
//    double dhBYdy = NA_VALUE;    /* Head gradient in y dirn. */
//    double Avg_Sf = NA_VALUE;    /* Head gradient in normal */
    double Dist2Edge[3];
    double Dist2Nabor[3];
    double FixPressure = NA_VALUE;  /* Pressure [Pa]*/
    double FixGamma = NA_VALUE;     /* Psychrometric Constant [kPa °C-1] */
    double AquiferDepth = NA_VALUE; // Zmax - Zmin
    double WetlandLevel = NA_VALUE; // Aquiferdepth - infD
    double RootReachLevel = NA_VALUE; //Aquiferdepth - RzD
    double MacporeLevel = NA_VALUE; //Aquiferdepth - macD
    double avgRough[3];
    double depression = 0.0002; //Depression value. No overland flow before filling the depression. Default = 0.2 mm.
    
    int iupdGW[3] = {0,0,0}; /* whether the Groundwater Flux value on this J is update */
    int iupdSF[3] = {0,0,0};  /* whether the Surface Flux value on this J is update */
    /* Value must be updated each loop */
    double u_qi; /* infiltration from surface to unsat zone */
    double u_qr; /* recharge into ground water */
    /* GW flow*/
    double u_effKH; /* Horizontal effective flow, for gw*/
    double u_satn; /* Saturation ratio */
    double u_wf = 0.;
    double u_deficit; /* deficit. aquiferdepth - Ygw */
private:
    /* Infiltration */
    //    double u_Yus; /* Water hea\d of the unsat zone from zmin*/
    double u_phius; /* pressure head of the unsat zone from zmin*/
    //    double u_Hus; /* Water head of the unsat zone from z = 0*/
    double u_Ginfi; /* Gradient to infiltration*/
    double u_satKr; /* ratio to effective Infiltration K */
    double u_effkInfi; /* Effective K for infiltration */
    /* Recharge */
//    double u_effkRech; /* Effective K for recharge */
//    double u_Grech; /* Gradient for recharge */
//    double u_ThetaFC; /* Field capacity */
    double u_Theta; /* Soil Moisture Content */
    
    /*==== Methods ===================*/
public:
    void InitElement();
    void copyGeol(Geol_Layer *);
    void copySoil(Soil_Layer *);
    void copyLandc(Landcover *);
    void applyGeometry(_Node *Node);
    void applyNabor(_Node *Node, _Element *Ele);
    void updateElement(double Ysurf, double Yunsat, double Ygw);
    void updateWF(double dh, double dt);
    void Flux_InfiRech(double Ysurf, double Yunsat, double Ygw, double netprcp);
    void printHeader(FILE *fp);
    void printInfo(FILE *fp);
};
#endif /* Element_hpp */

