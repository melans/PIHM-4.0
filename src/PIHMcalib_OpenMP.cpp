/*******************************************************************************
 * File        : PIHMcalib_OpenMP.cpp                                          *
 * Version     : Oct, 2018 (PIHM++ 4.0)                                        *
 * Function    : PIHM (Penn State Integrated Hydrologic Model                  *
 * Developer   :     Lele Shu (lele.shu@gmail.com)                             *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *..............Functionality of the Calibration tool..........................*
 * 0) CMA-ES (Covariance Matrix Adaption Evolution Strategy
 * 1) Support OpenMP parrellel. Parrellelism is only applied to CMA-ES section.
    e.g. not inside of PIHM simulation. If you have many cores/nodes resources
    for calibration, you can try PIHMcalib_MPI which support both OpenMP and OpenMPI.
 * 2) There are several objective functions are available, such as NSE, RMSE, R2
    nRMSE, Bias, etc.
 * 3) All the samples of parameters are saved for deep analysis.
 * 4) lambda > 10 is recommended, to increase the odd of global optimization.
 * 5) Final calibration file is exported as Final.calib.
 * 6) Observation vs Simulation data are saved.
 *******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <math.h>
#include "cmaes_interface.hpp"
#include "omp.h"
#include "ModelConfigure.hpp"
#include "Model_Data.hpp"
#include "ModelCalibration.hpp"
#include "pihm.hpp"
#include "Function_misc.hpp"
#include "print.hpp"
#include "CommandIn.hpp"
int main(int argc, char **argv) {
    CommandIn CLI;
    FileIn *fin = new FileIn();
    FileOut *fout = new FileOut();
    CLI.parse(argc, argv);
    CLI.setFileIO(fin, fout);
    
    ModelCalibration CB;
    
    CB.gc0.read(fin->file_calib); /* The calibration data */
    CB.init(CLI.n_lambda, fin->outpath, fin->file_calib, CLI.dir_cmaes);
    int     NDimPara = CB.N; /* Number of parameters */
    int     num_threads = CB.lamda;
    int     i;
    int     rank;
    
//    omp_set_num_threads(num_threads);
//    num_threads = omp_get_max_threads();
    FileIn *ifin = new FileIn[num_threads];
    FileOut *ifout = new FileOut[num_threads];

    int         nbetter = 0, Nrepeat = 0;
    int         Niter = 0;
    int         buf_stopsign = 1;       /* LOOP stops when  buf_stopsign = 0; */
    int         iBest = 0;     /* Best dataset, range in (1 ~ lamda) */
    double      BestObjVal = 1.0E16;    /* Minimum Object Value */
    double *    arFunvals = NULL;
    double *const * pop;        /* Population of calibration dataset */
    double          objValue[num_threads];
    cmaes_t         *evo = new cmaes_t();/* an CMA-ES type struct or "object" */
    string  str;
    
    for(int i = 0; i < num_threads; i++){
        objValue[i] = 0.;
        ifin[i].setInFilePath(fin->inpath, fin->projectname);
        str = "_" + to_string(i + 1);
        ifout[i].setsuffix(str.c_str());
        ifout[i].setOutFilePath(fout->outpath, fout->projectname);
    }
    /* Initialize everything into the struct evo, 0 means default */
    BestObjVal = 1e16;
    iBest = 0;
    arFunvals = cmaes_init(evo, NDimPara, NULL, NULL, 0, num_threads, CB.fn_cinit);
    cmaes_ReadSignals(evo, CB.fn_csignal);  /* write header and initial values */
    buf_stopsign=(cmaes_TestForTermination(evo)==NULL) ? 1 : 0;
    
    /* Iterate until stop criterion holds */
    
    while(buf_stopsign){
        /* generate lambda new search points, sample population */
        printf("\n\n***************************\n");
        pop = cmaes_SamplePopulation(evo);            /* do not change content of pop */
        for (i = 0; i < cmaes_Get(evo, "popsize"); ++i) {
            while ( !is_feasible(pop[i], NDimPara, -1., 1.) ){
                cmaes_ReSampleSingle(evo, i);
            }
        }
        printf("\t Step = %d \t sigma = %f\n", ++Niter, cmaes_Get(evo, "sigma"));
        CB.setCalibration(pop); /* Convert pop to parameter ranges and to .calib files.*/
        
        for(int i = 0; i < num_threads; i++){
            CB.setIO(&ifin[i], &ifout[i], i);
        }
#pragma omp parallel  default(shared)  shared(ifin, ifout, objValue) private(i, rank) num_threads(num_threads)
        {
#pragma omp for
            for (int i = 0; i < num_threads; i++) {
                /* Let all threads wait until pop ready in Rank0. */
                //        MPI_Barrier(MPI_COMM_WORLD);
                /******************************************/
                rank = omp_get_thread_num();
                /* Core PIHM simulations */
                printf("Process %d:\t%s\t%s\n", rank, ifin[rank].file_calib, ifout[rank].outpath);
                objValue[rank] = PIHM(&ifin[rank], &ifout[rank]);
                printf("Process %d: gof = %f\n", rank, objValue[rank]);
                CB.printObjValue(objValue[rank], rank);
                arFunvals[rank] = objValue[rank];
                /* Done PIHM simulations */
                /******************************************/
            } //end of for loop.
        }
        printf("\tAll process Finished. \n");
        nbetter = 0;
        for(i = 0; i < num_threads; i++){
            printf("Process %d's value is %f\n",i,arFunvals[i]);
            if(BestObjVal>arFunvals[i]){
                iBest = i;
                BestObjVal=arFunvals[iBest];
                CB.outBestData(iBest);
                nbetter++;
            }
        }
        if(nbetter > 0){
            Nrepeat = 0;
            printf("Minimum value of objective function %lf,which is found in the foler %d\n",BestObjVal,iBest);
        }else{
            Nrepeat++;
            printf("No Better simulation found. Best ObjValue = %f (Repeat %d Time(s))\n", BestObjVal, Nrepeat);
            if(Nrepeat>num_threads){
                printf("WARNING: Model may stucks on the setting, check the data and try again\n");
            }
            
            cmaes_UpdateDistribution(evo, arFunvals);
            /* read instructions for printing output or changing termination conditions */
            cmaes_ReadSignals(evo, CB.fn_csignal);
            buf_stopsign = (cmaes_TestForTermination(evo)==NULL) ? 1 : 0;
            if(Nrepeat > num_threads){
                buf_stopsign = 0;
            }
            printf("stopsign=%d\n",buf_stopsign);
        } // end of nbetter;
        
    } //end of while(buf_stopsign)
    printf("Stop:\n%s\n",  cmaes_TestForTermination(evo)); /* print termination reason */
    cmaes_WriteToFile(evo, "all", CB.fn_result);         /* write final results */
    
    /* get best estimator for the optimum, xmean */
//    FinalCalib = cmaes_GetNew(evo, "xmean"); /* "xbestever" might be used as well */
    cmaes_exit(evo); /* release memory */
    /* do something with final solution and finally release memory */
    delete[] ifin;
    delete[] ifout;
    delete fin;
    delete fout;
    return 0;
}


