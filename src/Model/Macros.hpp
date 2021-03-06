//  Macros.h
//
//  Created by Lele Shu on 2/16/17.
//  Copyright (c) 2017 Lele Shu. All rights reserved.

#ifndef PIHM3_0_Macros_h
#define PIHM3_0_Macros_h
#include <math.h>
#include <vector>

#ifdef _PIHMOMP
#include "omp.h"
#include "nvector_openmp.h" /* contains the definition of type N_Vector for openmp */
#include "nvector/nvector_openmp.h" /* serial N_Vector types, fcts., macros */
#define SET_VALUE(v, i) NV_Ith_OMP(v,i)
#else
#include "nvector_serial.h" /* contains the definition of type N_Vector */
#define SET_VALUE(v, i) NV_Ith_S(v,i)
#endif

///*========index===============*/
#define iSF     i
#define iUS     i + NumEle
#define iGW     i + 2 * NumEle
#define iRIV    i + 3 * NumEle
#define iLAKE    i + 3 * NumEle + NumRiv
#define iDownStrm Riv[i].down - 1


/*========Misc constant===============*/
#define MAXLEN 2048  /*Max Str Length*/
#define EPSILON 0.005
#define EPS_DOUBLE 1.0e-10 // precision of double type in 1.e-14.
#define EPS_SLOPE   0.05e-6
#define MINPSI -1000000
#define FieldCapacityRatio 0.75
#define MAXQUE 10000
#define Nforc 5
#define i_prcp 1
#define i_temp 2
#define i_rh 3
#define i_wind 4
#define i_rn 5

/*========Physical Constant value===============*/
#define PI 3.1415926
#define MINRIVSLOPE 1e-4
#define C_air 1004.0
#define THRESH 0.0

#define GRAV 9.8		/* m/s^2 Note the dependence on physical units */

#define C_air 1004.0
#define Lv 2503000.0 //Volumetric latent heat of vaporization. Energy required per water volume vaporized. (Lv = 2453 MJ m−3 on wiki)
#define SIGMA 3.402e-6
#define R_dry 287.04
#define R_v 461.5
#define Ts  -3.0  //Threshold for Snow
#define Tr  1.0
#define To  0.0

#define IC_MAX 0.0002  // Maximum Interception on caonpy.
#define IC_MAX_SNOW  0.003

#define CKconst  273.15 /* Kelvin Constant */
#define VON_KARMAN     0.4        /* Von Karman's constant */
#define HeightWindMeasure   10  /* Height Wind/Relative Humidity Measure */
#define Cp  1.013e-3    /* cp specific heat at constant pressure, 1.013E-3 [MJ kg-1 °C-1] Allen(1998) eq(8) */
//#define LAMBDA 2.45 /* λ latent heat of vaporization, 2.45 [MJ kg-1] Allen(1998) eq(8) */
//#define VAPRATIO 0.622 /* ε ratio molecular weight of water vapour/dry air = 0.622. Allen(1998) eq(8) */

/*========ERROR CODE===============*/
#define ERRSUCCESS  0
#define ERRNAN      10
#define ERRFileIO   12
#define ERRDATAIN   13
#define ERRCVODE    19
#define ERRCONSIS   20
#define NA_VALUE -9999
/*=======================*/
#define ID 94 // debug only.
extern int debug_mode;
extern int dummy_mode;
extern int verbose_mode;
extern int sinks_remove;
extern int smooth_river;
extern int quiet_mode;
extern int ilog;

extern double *uYsf;
extern double *uYus;
extern double *uYgw;
extern double *uYriv;
extern double *uYlake;
//extern double *DY;
extern double *globalY;

extern double timeNow;

#endif
