//
//  Function_misc.hpp
//  PIHMcalib
//
//  Created by Lele Shu on 9/25/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#ifndef Function_misc_hpp
#define Function_misc_hpp

#include <stdio.h>

//#define SLEEP_TIME 2
//#define SLEEP_TIME1 3
//extern int errno;


#define RootRank 0
int ifRoot(int rank);
int is_feasible(double const *x, int N, double xmin, double xmax);


#endif /* Function_misc_hpp */
