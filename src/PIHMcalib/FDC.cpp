//
//  FDC.cpp
//  PIHMcalib
//
//  Created by Lele Shu on 8/18/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#include "FDC.hpp"

/*********************FLOW DURATION CURVE********************************/

void FDC(double *data, double *px, double *py, int n){
    /* Flow Duration Curve */
    // copy x to px
    for(int i = 0; i < n; i++){
        py[i] = data[i] ;
    }
    quickSort(py, 0, n);
    double key;
    int nlt;
    for(int i = 0; i < n; i++){
        key = py[i];
        nlt = 0;
        for(int j = 0; j < n; j++){
            if( key >= py[i] ){
                nlt++;
                continue;
            }else{
                break;
            }
        }
        px[i] = (double) nlt / n;
    }
}
