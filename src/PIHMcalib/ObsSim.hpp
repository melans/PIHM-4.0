//  ObsSim.hpp
//  Created by Lele Shu (lele.shu@gmail.com) on 2018.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#ifndef ObsSim_hpp
#define ObsSim_hpp

#include <stdio.h>
#include "Macros.hpp"
#include "functions.hpp"
#include "GoodOfFit.hpp"

class ObsnSim{
public:
    int nTotal, nSpin, nCalib;
    int targetID = 0;
    double *time;
    double *sim;
    double *obs;
    double *ptr;
    double ctm, cta; /* Conversion factor multiple and plus. y = x * ctm + cta */
    double SimMean = 0.;
    double DayStart, DayEnd;
    GoodOfFit *gof;
    
    ObsnSim();
    ~ObsnSim();
    
    void Pointer2Sim(double *x, double ctMultiple, double ctPlus);
    void readobs(const char *fn);
    void pushsim(double t);
    void callGOF(void);
    void printGOF(FILE *fp);
    double getNSE(void);
    void printData(const char *fn);
private:
    int iNow = 0;
    void init(int n);
};
#endif /* ObsSim_hpp */
