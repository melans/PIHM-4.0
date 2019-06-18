//  Model_Control.hpp
//  Created by Lele Shu (lele.shu@gmail.com) on 2018.
//  Copyright © 2018 Lele Shu. All rights reserved.
//
#ifndef Model_Control_hpp
#define Model_Control_hpp

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "Macros.hpp"
#include "ModelConfigure.hpp"

#ifdef _CALIBMODE
#include "ObsSim.hpp"
#endif

class Print_Ctrl{
public:
    char    filename[MAXLEN];
    int     Interval = NA_VALUE;
    int     NumVar = NA_VALUE;
    double  **PrintVar = NULL;
    double  *buffer = NULL;
    int     Binary = 1;
    int     Ascii = 0;
    double  tau = 1440. / UNIT_C;    // time unit in calculation. [min]
    FILE    *fid_bin = NULL;
    FILE    *fid_asc = NULL;
    char    filea[MAXLEN];
    char    fileb[MAXLEN];
    
    Print_Ctrl();
    ~Print_Ctrl();
    void    open_file(int a, int b);
    void    PrintData (double dt, double t);
    void    Init(long st, int n, const char *s, int dt, double *x, int iFlux);
    void    InitIJ(long st, int n, const char *s, int dt, double **x, int j, int iFlux);
private:
    void    fun_printASCII(double t, double dt);
    void    fun_printBINARY(double t, double dt);
    void    close_file();
private:
    long   StartTime;
};
class PrintOutDt {
public:
    /* Element storage */
    int dt_ye_gw = 0;
    int dt_ye_surf = 0;
    int dt_ye_snow = 0;
    int dt_ye_ic = 0;
    int dt_ye_unsat = 0;
    
    /* Element Fluxes */
    int dt_qe_prcp = 1440; // default output PRCP.
    int dt_qe_infil = 0;
    int dt_qe_et[3] = {0,0,0};
    int dt_qe_rech = 0;
    int dt_qe_etp = 0;
    
    /* Element volume Fluxes */
    int dt_Qe_sub = 0;
    int dt_Qe_surf = 0;
    
    /* River Stage */
    int dt_yr_stage = 0;
    /* River volume Fluxes */
    int dt_Qr_up = 0;
    int dt_Qr_down = 0;
    int dt_Qr_sub = 0;
    int dt_Qr_surf = 0;
    
    /* Lake  stage */
    int dt_yl_stage     = 0;
    /* Lake  Fluxes */
    int dt_Ql_surf      = 0;
    int dt_Ql_sub       = 0;
    int dt_Ql_chn       = 0;

    int dt_ql_et        = 0;
    int dt_ql_prcp      = 0;
    int dt_Ql_rivin     = 0;
    int dt_Ql_rivout    = 0;
    
    void calibmode(int dt);
    void defaultmode();
};
class Control_Data : public PrintOutDt{
private:
    double DayStart = 0;
    double DayEnd = 10;
    //    int outtype;
    double a =1;    /* External time stepping controls */
    double b = 1;
public:
    int Verbose = 0;
    int Debug = 0;
    int Ascii = 0;
    int Binary = 1;
    
    //    int Solver;    /* Solver type */
    int NumSteps;    /* Number of external time steps
                      * (when results can be printed) for
                      * the whole simulation */
    int num_threads;    /* Number of Threads in OPENmp only*/
    int init_type = 3;    /* initialization mode */
    
    double abstol = 1.0e-4;    /* absolute tolerance */
    double reltol = 1.0e-3;    /* relative tolerance */
    double InitStep = 1e-5;    /* initial step size */
    double MaxStep = 300;/* Maximum step size */
    double ETStep = 30;    /* Step for et from interception */
    double StartTime = 0.;    /* Start time of simulation */
    double EndTime = 14400;/* End time of simulation */
    double dt = 1;
    double *Tout;
    
    //    globalCal Cal;    /* Convert this to pointer for localized  calibration */
    
    int NumPrint = 0;;
    
    Print_Ctrl PCtrl[100];
#ifdef _CALIBMODE
    ObsnSim CV;
#endif
    
    /* Methods */
    Control_Data();
    ~Control_Data();
    void ExportResults(double t);
    void updateSimPeriod(double day0, double day1);
    void read(const char *fn);
    void write(const char *fn);
    void getValue(const char *varname);
private:
    void updateSimPeriod();
} ;

#endif /* Model_Control_hpp */
