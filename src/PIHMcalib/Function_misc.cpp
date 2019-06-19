//
//  Function_misc.cpp
//  PIHMcalib
//
//  Created by Lele Shu on 9/25/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#include "Function_misc.hpp"
int ifRoot(int rank){
    if(rank == RootRank){
        return 1;
    }else{
        return 0;
    }
}

/* the objective (fitness) function to be minized */
/* the optimization loop */
int is_feasible(double const *x, int N, double xmin, double xmax){
    for(int i = 0; i < N; i++){
        if( x[i] < xmin || x[i] > xmax ){
            return 0;
        }
    }
    return 1;
}

