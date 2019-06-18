#include <stdio.h>
#include <iostream>
//#include "f_element.hpp"
//#include "f_River.hpp"
#include "f.hpp"
#include "IO.hpp"
#include "ModelConfigure.hpp"
#include "print.hpp"
#include "Macros.hpp"
#include "functions.hpp"
//#include "is_sm_et.hpp"
#include "cvode_config.hpp"
#include "Model_Data.hpp"
#include "ForcingData.hpp"
#include "FloodAlert.hpp"
#include "CommandIn.hpp"

double *uYsf;
double *uYus;
double *uYgw;
double *uYriv;
double *dy;
using namespace std;
double PIHM(FileIn *fin, FileOut *fout){
    double ret = 0.;
    Model_Data  *MD;        /* Model Data                */
    N_Vector    udata;
    N_Vector    du;
    
    void    *mem = NULL;
    SUNLinearSolver LS = NULL;
    int     flag;            /* flag to test return value */
    double  t, dt, tout;    /* stress period & step size */
    int NY = 0;
    int ierr = 0;
    /* allocate memory for model data structure */
    MD = new Model_Data();
    MD->loadinput(fin);
    MD->initialize();
    MD->CheckInputData();
    fout->updateFilePath();
    NY = MD->NumY;
    dy = new double[NY];
#ifdef _PIHMOMP
    omp_set_num_threads(MD->CS.num_threads);
    screeninfo("\nopenMP enabled. No of Threads = %d\n", MD->CS.num_threads);
    udata = N_VNew_OpenMP(NY, MD->CS.num_threads);
#else
    screeninfo("\nopenMP disabled\n");
    udata = N_VNew_Serial(NY);
    du = N_VNew_Serial(NY);
#endif
    MD->InitialCondition(fin, udata);
    
#ifdef _CALIBMODE
    MD->CS.CV.readobs(fin->file_obs);
    MD->CS.CV.Pointer2Sim(MD->QrivDown, 1, 0);
    MD->CS.updateSimPeriod(MD->CS.CV.DayStart, MD->CS.CV.DayEnd);
    MD->CS.calibmode(1440);
#else
    MD->initialize_output(fout);
    MD->PrintInit(fout->Init_bak);
    MD->InitFloodAlert(fout->floodout);
#endif
    SetCVODE(mem, f, MD, udata, LS);
    
    /* set start time */
    t = MD->CS.StartTime;
    double tnext = t;
    //CheckInput(MD, &CS);
    /* start solver in loops */
    getSecond();
    MD->modelSummary(fin, 0);
    MD->debugData(fout->outpath);
    MD->gc.write(fout->Calib_bak);

    FILE *fid = fopen("DY_debug.dat", "wb");
    fclose(fid);
    
    for (int i = 0; i < MD->CS.NumSteps && !ierr; i++) {        /* inner loops to next output points with ET step size control */
        while (t < tnext ) {
            if (t + MD->CS.ETStep >=tnext) {
                tout = tnext;
            } else {
                tout = t + MD->CS.ETStep;
            }
            dt = tout - t;
            MD->updateforcing(t);
            /* calculate Interception Storage */
            MD->EvapoTranspiration(udata, t, dt);
            MD->updateWF(dt);
//            fSolver(MD, udata, t, tout);
            flag = CVode(mem, tout, udata, &t, CV_NORMAL);
            check_flag(&flag, "CVode", 1);
        }
        tnext += MD->CS.MaxStep;
        MD->summary(udata);
#ifdef _CALIBMODE
        MD->CS.CV.pushsim(t);
        if (atInterval(t, 1440 * 30)) {
            MD->PrintInit(fout->Init_update);
        }
#else
//        f(t, udata, du, MD);
        MD->CS.ExportResults(tnext);
//        flag = ScreenPrint(t, MD->CS.ETStep, i, MD->CS.NumSteps);
        flag = ScreenPrint(t, MD->CS.ETStep, i, MD->CS.NumSteps, MD->nFCall);
        if (flag) {
            if (atInterval(t, 1440)) {
                MD->PrintInit(fout->Init_update);
            }
        }
        MD->flood->FloodWarning(t);
#endif
    }
#ifndef _CALIBMODE
    PrintFinalStats(mem);
#endif
    MD->modelSummary(fin, 1);
    /* Free memory */
    N_VDestroy_Serial(udata);
    N_VDestroy_Serial(du);
    /* Free integrator memory */
    CVodeFree(&mem);
    
#ifdef _CALIBMODE
    if(ierr){
        ret = 9999;
    }else{
        MD->CS.CV.printData(fout->obs_sim);
        MD->CS.CV.callGOF();
        MD->CS.CV.printGOF(stdout);
        //    printf("OBJvalue = %f\n", MD->CS.CV.gof->gof_nRMSE(1) );
        ret = MD->CS.CV.gof->gof_nRMSE(1);
        if(ret < 0.01){
            /* if nRMSE is less than 1%, it is good enough */
            ret = 0.0;
        }
    }
#endif
    
    delete MD;
    delete dy;
    return ret;
}


int PIHM(int argc, char *argv[]){
    CommandIn CLI;
    FileIn *fin = new FileIn;
    FileOut *fout = new FileOut;
    CLI.parse(argc, argv);
    CLI.setFileIO(fin, fout);
    
    PIHM(fin, fout);
    
    delete fin;
    delete fout;
    return 0;
}

