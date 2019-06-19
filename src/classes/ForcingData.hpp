//  ForcingData.hpp
//  Created by Lele Shu (lele.shu@gmail.com) on 2018.
//  Copyright © 2018 Lele Shu. All rights reserved.
//
#ifndef ForcingData_hpp
#define ForcingData_hpp
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include "Macros.hpp"
#include "IO.hpp"
#include "functions.hpp"

class _DataQueue {
public:
    std::string fn;
    
    _DataQueue();
    ~_DataQueue();
    void    readDimensions();
    void    read_csv();
    void    tsd_interpolation(double t);
    double  getX(double t, int column);
    void    applyCalib(double prcp, double temp);
    void    movePointer(double t);
    void    initialize(int n);
    void    checkValue(int icol, double xmin, double xmax, const char *varname);
    int     get_Ncol();
private:
    long StartTime;
    int ncol = 0;
    int Length;
    int eof;
    int iNow, iNext;
    int nQue = 0;
    double *ts[MAXQUE + 1];
    int pRing[MAXQUE + 1];
    
    //    void    buildfn(std::string fforc);
    double  interpolation(double t);
};

void CheckFile(std::ifstream * fp, const char *s);
#endif                /* ForcingData_hpp */


