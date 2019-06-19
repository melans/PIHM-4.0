//
//  quicksort.cpp
//  PIHMcalib
//
//  Created by Lele Shu on 8/18/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#include "quicksort.hpp"

void swap(double * a, double * b) {
    double tmp = * a;
    * a = * b;
    * b = tmp;
}
int partition(double* fptr, int left, int right)
{
    int iwall = left;     //The location of Wall, Initialization.
    double pivot = fptr[right];   //Use last var as the intial guess
    for (int i = left; i < right; i++){
        
        if(pivot > fptr[i])
        {
            swap(fptr+i,fptr+iwall);
            iwall++;
        }
    }
    swap(fptr+iwall,fptr+right);
    return iwall;    //return the location of wall.
}
void quickSort(double x[], int left, int right){
    if (left < right){
        /* pi is partitioning index, x[p] is nowat right place */
        int pi = partition(x, left, right);
        quickSort(x, left, pi - 1);  // Before pi
        quickSort(x, pi + 1, right); // After pi
    }
}
