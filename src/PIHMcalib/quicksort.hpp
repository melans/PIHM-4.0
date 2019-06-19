//  quicksort.hpp
//  Created by Lele Shu (lele.shu@gmail.com) on 2018.
//  Copyright © 2018 Lele Shu. All rights reserved.
//
#ifndef quicksort_hpp
#define quicksort_hpp

#include <stdio.h>
void swap(double * a, double * b);
int partition(double * fptr, int left, int right);
void quickSort(double x[], int low, int high);

#endif /* quicksort_hpp */
