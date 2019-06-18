/*******************************************************************************
 * File        : PIHMcalib_MPI.cpp                                             *
 * Version     : Oct, 2018 (PIHM++ 4.0)                                        *
 * Function    : PIHM (Penn State Integrated Hydrologic Model                  *
 * Developer   :     Lele Shu (lele.shu@gmail.com)                             *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *..............Functionality of the Calibration tool..........................*
 * 0) CMA-ES (Covariance Matrix Adaption Evolution Strategy
 * 1) Support OpenMPI parrellel. Parrellelism is NOT ONLY applied to CMA-ES section,
    BUT ALSO inside of PIHM simulation with (OpenMP). So you need install both
    OpenMPI and OpenMP to accelarate the calibration. If you have limited
    resources, particularly processors on your machine, calibration with OpenMP
    is your better option.
 * 2) There are several objective functions are available, such as NSE, RMSE, R2
    nRMSE, Bias, etc.
 * 3) All the samples of parameters are saved for deep analysis.
 * 4) lambda > 10 is recommended, to increase the odd of global optimization.
 * 5) Final calibration file is exported as Final.calib.
 * 6) Observation vs Simulation data are saved.
 *******************************************************************************/
#include <stdio.h>
#include <stdlib.h> /* free() */
#include "cmaes_interface.hpp"
#include <unistd.h>
#include <sys/types.h>
//#include <errno.h>
#include <math.h>
#include "mpi.h"

#include "ModelConfigure.hpp"
#include "Model_Data.hpp"
#include "ModelCalibration.hpp"
#include "pihm.hpp"
#include "Function_misc.hpp"
#include "print.hpp"
#include "CommandIn.hpp"
int main(int argc, char **argv) {
    CommandIn CLI;
    int         i, flag, namelen;
    int         numprocs, rank;
    char        processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   /* Get the Rank in MPI */
    MPI_Get_processor_name(processor_name, &namelen);
    //    MPI_Status status;
    if(rank == 0){
        PIHMlogo();
        fprintf(stdout, "MPI mode enable\n");
    }
    FileIn *fin = new FileIn();
    FileOut *fout = new FileOut();
    CLI.parse(argc, argv);
    CLI.setFileIO(fin, fout);
    
    ModelCalibration CB;
    
    CB.gc0.read(fin->file_calib); /* The calibration data */
    CB.init(numprocs, fout->outpath, fin->file_calib, CLI.dir_cmaes);
    int         NDimPara = CB.N; /* Number of parameters */
    
    int             nbetter = 0, Nrepeat = 0;
    int             buf_stopsign = 1;       /* LOOP stops when  buf_stopsign = 0; */
    int             Niter = 0;
    int             iBest = 0;     /* Best dataset, range in (1 ~ lamda) */
    double          BestObjVal = 1.0E16;    /* Minimum Object Value */
    
    double *        arFunvals = NULL;
    double *const * pop;        /* Population of calibration dataset */
    double *        FinalCalib;     /* Final (best) calibration dataset */
    double          objValue[numprocs];
    cmaes_t         *evo = new cmaes_t();/* an CMA-ES type struct or "object" */
    
    for(int i = 0; i < numprocs; i++){
        objValue[i] = 0.;
    }
    /* Initialize everything into the struct evo, 0 means default */
    if( ifRoot(rank) ){
        BestObjVal = 1e16;
        iBest = 0;
        arFunvals = cmaes_init(evo, 0, NULL, NULL, 0, 0, CB.fn_cinit);
        cmaes_ReadSignals(evo, CB.fn_csignal);  /* write header and initial values */
        buf_stopsign=(cmaes_TestForTermination(evo)==NULL) ? 1 : 0;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    flag=MPI_Bcast(&buf_stopsign, 1 , MPI_INT, 0, MPI_COMM_WORLD);
    /* Iterate until stop criterion holds */
    while(buf_stopsign){
        /* Here you may resample each solution point pop[i] until it
         becomes feasible, e.g. for box constraints (variable
         boundaries). function is_feasible(...) needs to be
         user-defined.
         Assumptions: the feasible domain is convex, the optimum is
         not on (or very close to) the domain boundary, initialX is
         feasible and initialStandardDeviations are sufficiently small
         to prevent quasi-infinite looping.
         */
        if( ifRoot(rank) ){
            printf("\n\n***************************\n");
            /* generate lambda new search points, sample population */
            pop = cmaes_SamplePopulation(evo);            /* do not change content of pop */
            for (i = 0; i < cmaes_Get(evo, "popsize"); ++i) {
                while ( !is_feasible(pop[i], NDimPara, 0, 1) ){
                    cmaes_ReSampleSingle(evo, i);
                }
            }
            printf("\t Step = %d \t sigma = %f\n", ++Niter, cmaes_Get(evo, "sigma"));
            CB.setCalibration(pop); /* Convert pop to parameter ranges and to .calib files.*/
        }
        /* Let all threads wait until pop ready in Rank0. */
        MPI_Barrier(MPI_COMM_WORLD);
        /******************************************/
        /* Core PIHM simulations */
        CB.setIO(fin, fout, rank); /* Specify the input .calib files and output folder*/
        printf("Process %d:\t%s\t%s\n", rank, fin->file_calib, fout->outpath);
        objValue[rank] = PIHM(fin, fout);
        CB.printObjValue(objValue[rank], rank);
        printf("Process %d: gof = %f\n", rank, objValue[rank]);
        for(i = 0; i < numprocs; i++){
            MPI_Reduce(&objValue[i], &arFunvals[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        /* Done PIHM simulations */
        /******************************************/

        
        flag=MPI_Barrier(MPI_COMM_WORLD);
        if( ifRoot(rank) ) {
            nbetter = 0;
            for(i = 0; i < numprocs; i++){
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
                if(Nrepeat>numprocs){
                    printf("WARNING: Model may stucks on the setting, check the data and try again\n");
                }
            }
            
            cmaes_UpdateDistribution(evo, arFunvals);
            /* read instructions for printing output or changing termination conditions */
            cmaes_ReadSignals(evo, CB.fn_csignal);
            buf_stopsign = (cmaes_TestForTermination(evo)==NULL) ? 1 : 0;
            if(Nrepeat > numprocs){
                buf_stopsign = 0;
            }
            printf("stopsign=%d\n",buf_stopsign);
        } // end of ifRoot(rank)
        flag=MPI_Bcast(&buf_stopsign,1 , MPI_INT, 0, MPI_COMM_WORLD);
        flag=MPI_Barrier(MPI_COMM_WORLD);
    } //end of while(buf_stopsign)
    if( ifRoot(rank) )  printf("Stop:\n%s\n",  cmaes_TestForTermination(evo)); /* print termination reason */
    if( ifRoot(rank) )  cmaes_WriteToFile(evo, "all", CB.fn_result);         /* write final results */
    
    /* get best estimator for the optimum, xmean */
    if( ifRoot(rank) )   FinalCalib = cmaes_GetNew(evo, "xmean"); /* "xbestever" might be used as well */
    if( ifRoot(rank) )   cmaes_exit(evo); /* release memory */
    /* do something with final solution and finally release memory */
    //    if( ifRoot(rank) )  free(xfinal);
    MPI_Finalize();
    
    delete fin;
    delete fout;
    return 0;
}


