/*******************************************************************************
 * File        : pihm.hpp                                                      *
 * Version     : Oct, 2018 (PIHM++ 4.0)                                        *
 * Function    : PIHM (Penn State Integrated Hydrologic Model                  *
 * Developer of PIHM++ 4.0:      Lele Shu (lele.shu@gmail.com)                 *
 * Developer of PIHM 3.0:        Gopal Bhatt (gopal.bhatt@psu.edu)             *
 * Developer of PIHM 2.0:        Mukesh Kumar (muk139@psu.edu)                 *
 * Developer of PIHM 1.0:        Yizhong Qu   (quyizhong@gmail.com)            *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *.............................................................................*
 *******************************************************************************/

#ifndef pihm_hpp
#define pihm_hpp

#include "IO.hpp"

int PIHM(int argc, char *argv[]);
double PIHM(FileIn *fin, FileOut *fout);

#endif /* pihm.hpp */
