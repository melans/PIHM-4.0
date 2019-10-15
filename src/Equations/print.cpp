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


