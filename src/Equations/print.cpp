/*******************************************************************************
 * File        : print.c	                                                   *
 * Version     : Oct, 2018 (PIHM++ 4.0)                                        *
 * Function    : PIHM (Penn State Integrated Hydrologic Model                  *
 * Developer of PIHM++ 4.0:      Lele Shu (lele.shu@gmail.com)                 *
 * Developer of PIHM 3.0:        Gopal Bhatt (gopal.bhatt@psu.edu)             *
 * Developer of PIHM 2.0:        Mukesh Kumar (muk139@psu.edu)                 *
 * Developer of PIHM 1.0:        Yizhong Qu   (quyizhong@gmail.com)            *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *                                                                             *
 *..............MODIFICATIONS/ADDITIONS in PIHM++ 4.0...........................*
 * a) This file is downgraded from Version 1.0, as no ancillary results are    *
 *    being output			                                                   *
 * b) Only state variables and flux to/in/accross river and its bed are being  *
 *    output							                                       *
 * c) Addition of Average Function to output average variables at regular time *
 *    intervals								                                   *
 *******************************************************************************/

#include "print.hpp"
void Model_Data::summary (N_Vector udata){
    double  *Y;
#ifdef _PIHMOMP
    Y = NV_DATA_OMP(udata);
#else
    Y = NV_DATA_S(udata);
#endif
    for (int i = 0; i < NumEle; i++){
        yEleSurf[i] = Y[iSF];
        yEleUnsat[i] = Y[iUS];
        yEleGW[i] = Y[iGW];
        uYsf[i] = Y[iSF];
        uYus[i] = Y[iUS];
        uYgw[i] = Y[iGW];
    }
    for (int i = 0; i < NumRiv; i++){
        yRivStg[i] = Y[iRIV];
        uYriv[i] = Y[iRIV];
    }
}
void Model_Data::summary (N_Vector udata1, N_Vector udata2){
    double  *Y1, *Y2;
#ifdef _PIHMOMP
    Y1 = NV_DATA_OMP(udata1);
    Y2 = NV_DATA_OMP(udata2);
#else
    Y1 = NV_DATA_S(udata1);
    Y2 = NV_DATA_S(udata2);
#endif
    for (int i = 0; i < NumEle; i++){
        yEleSurf[i] = Y1[iSF];
        yEleUnsat[i] = Y1[iUS];
        yEleGW[i] = Y1[iGW];
        
        uYsf[i] = Y1[iSF];
        uYus[i] = Y1[iUS];
        uYgw[i] = Y1[iGW];
    }
    for (int i = 0; i < NumRiv; i++){
        yRivStg[i] = Y2[i];
        if(Y2[i] < 0){
            printf("Yriv[%d]=%f\n", i+1, Y2[i]);
        }
        uYriv[i] = Y2[i];
    }
}

void PIHMlogo(void){
    printf ("\n");
    printf ("\t########  ####  ##     ##  ##     ##             #    \n");
    printf ("\t##     ##  ##   ##     ##  ###   ###             #    \n");
    printf ("\t##     ##  ##   ##     ##  #### ####     #    ####### \n");
    printf ("\t########   ##   #########  ## ### ##     #       #    \n");
    printf ("\t##         ##   ##     ##  ##     ##  #######    #    \n");
    printf ("\t##         ##   ##     ##  ##     ##     #            \n");
    printf ("\t##        ####  ##     ##  ##     ##     # Verion 4.0 \n");
    
    printf ("\n\t\tThe Penn State Integrated Hydrologic Model v4.0\n");
    
#ifdef _CALIBMODE
    printf("\nCalibration mode enable.\n");
#endif
#ifdef _DEBUG
    printf("\nModel-debug mode enable.\n");
#endif
    
#ifdef _PIHMOMP
    printf("\nopenMP enabled. Maximum Threads = %d\n", omp_get_max_threads());
#else
    printf("\nopenMP disabled.\n");
#endif
}


void Model_Data::PrintInit (const char *fn)
{
    FILE           *fp;
    fp = fopen (fn, "w");
    CheckFile(fp, fn);
    
    fprintf (fp, "%d\t%d\n", NumEle, 6);
    fprintf (fp, "%s\t%s\t%s\t%s\t%s\t%s\n","Index",
             "Canopy", "Snow", "Surface", "Unsat", "GW");
    for (int i = 0; i < NumEle; i++){
        fprintf (fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\n", i+1, yEleIS[i], yEleSnow[i], yEleSurf[i], yEleUnsat[i], yEleGW[i]);
    }
    fprintf (fp, "%d\t%d\n", NumRiv, 2);
    fprintf (fp, "%s\t%s\n", "Index", "Stage");
    for (int i = 0; i < NumRiv; i++){
        fprintf (fp, "%d\t%lf\n", i+1, yRivStg[i]);
    }
    if(NumLake > 0){
        fprintf (fp, "%d\t%d\n", NumLake, 2);
        fprintf (fp, "%s\t%s\n", "Index", "LakeStage");
        for (int i = 0; i < NumLake; i++){
            fprintf (fp, "%d\t%lf\n", i+1,yLakeStg[i]);
        }
    }
    fclose (fp);
}

