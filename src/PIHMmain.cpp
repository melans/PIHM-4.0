/*******************************************************************************
 * File        : PIHMmain.cpp                                                  *
 * Version     : Oct, 2018 (PIHM++ 4.0)                                        *
 * Function    : PIHM (Penn State Integrated Hydrologic Model                  *
 * Developer of PIHM++ 4.0:      Lele Shu (lele.shu@gmail.com)                 *
 * Developer of PIHM 3.0:        Gopal Bhatt (gopal.bhatt@psu.edu)             *
 * Developer of PIHM 2.0:        Mukesh Kumar (muk139@psu.edu)                 *
 * Developer of PIHM 1.0:        Yizhong Qu   (quyizhong@gmail.com)            *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *..............MODIFICATIONS/ADDITIONS in PIHM++ 4.0..........................*
 * 0) Change the language and structure of code from C to C++;
 * 1) Update the CVODE from v2.2 to v3.x
 * 2) Support OpenMP Parrallel computing
 * 3) Change the input/output format. Check the Manual of PIHM++ v4.0;
 * 4) Change the structure of River
 * 5) The functions to handle the time-series data, including forcing, LAI,
 *    Roughness Length, Boundary Condition, Melt factor
 * 6) Add Lakes into the hydrologic process
 * 7) CMA-ES calibration, with either OpenMP or OpenMPI.
 *******************************************************************************/
#include "pihm.hpp"
#include "print.hpp"
int main(int argc, char *argv[]){
    PIHMlogo();
    PIHM(argc, argv);
    return 0;
}

