//
//functions.cpp
//
//Created by Lele Shu on 6 / 23 / 18.
// Copyright � �2018 Lele Shu.All rights reserved.
//

#include "functions.hpp"

void myexit(int flag){
    switch (flag) {
        case ERRFileIO:
            fprintf(stderr, "\nEXIT with error code %d(FILEIO)\n", flag);
            break;
        case ERRNAN:
            fprintf(stderr, "\nEXIT with error code %d(NAN/INF VALUE)\n", flag);
            break;
        case ERRCVODE:
            fprintf(stderr, "\nEXIT with error code %d(CVODE)\n", flag);
            break;
        case ERRCONSIS:
            fprintf(stderr, "\nEXIT with error code %d(Data Consistency)\n", flag);
            break;
        case ERRDATAIN:
            fprintf(stderr, "\nEXIT with error code %d(Data validation)\n", flag);
            break;
        case ERRSUCCESS:
            fprintf(stderr, "\nEXIT without any error (%d)\n", flag);
            break;
        default:
            fprintf(stderr, "\nEXIT with error code %d(Undefined error)\n", flag);
            break;
    }
    fprintf(stderr, "\n\n\n");
    exit(flag);
}
int checkRange(int x, int xmin, int xmax, int i, const char *s){
    if(x < xmin || x > xmax){
        fprintf(stderr, "\nInvalid value %d for %s(%d) is out of range [%d, %d]\n", x, s, i+1,  xmin, xmax);
        return 1;
    }else{
        return 0;
    }
}
int checkRange(double x, double xmin, double xmax, int i, const char *s){
    if(x < xmin || x > xmax){
        fprintf(stderr, "\nInvalid value %f for %s(%d) is out of range [%g, %g]\n", x, s, i+1,  xmin, xmax);
        return 1;
    }else{
        return 0;
    }
}
double timeInterp(double t, double t0, double t1, double x0, double x1){
    double dt = t1 - t0;
    if (fabs(dt) < EPS_DOUBLE) {
        return x1;
    } else {
        return ((t1 - t) * x0 + (t - t0) * x1) / dt;
    }
}

double min(double *x, int n){
    double  ret;
    ret = x[0];
    for(int i = 0; i < n; i++){
        if(ret > x[i]){
            ret = x[i];
        }
    }
    return ret;
}
double max(double *x, int n){
    double  ret;
    ret = x[0];
    for(int i = 0; i < n; i++){
        if(ret < x[i]){
            ret = x[i];
        }
    }
    return ret;
}
void CheckNANi(double x, int i, const char *s)
{
    if (isnan(x) || isinf(x)) {
        printf("\nERROR: NAN error for %s %d\n", s, i + 1);
        myexit(ERRNAN);
    }
}
void CheckNANij(double x, int i, const char *s)
{
    if (isnan(x) || isinf(x)) {
        printf("\nERROR: NAN error for %s %d\n", s, i + 1);
        myexit(ERRNAN);
    }
}
double getSecond(void)
{
#ifdef _PIHMOMP
    static double t0 = 0.;
    double t1;
    double sec;
    t1 = omp_get_wtime();
    sec = (double)(t1 - t0);
    t0 = t1;
#else
    static clock_t t0 = 0;
    clock_t t1;
    double sec;
    t1 = clock();
    sec = (double)(t1 - t0) / CLOCKS_PER_SEC;
    t0 = t1;
#endif
    return sec;
}


void CheckNonZero(double x, int i, const char *s)
{
    if (x <= 0.0 || isnan(x) || isinf(x) || fabs(x - NA_VALUE) < EPS_DOUBLE) {
        printf("ERROR: Value %e for %s of Element %d is not allowed. Please check again.\n", x, s, i + 1);
        myexit(ERRNAN);
    }
}
void CheckNonZero(int x, int i, const char *s)
{
    if (x <= 0) {
        printf("ERROR: Value %d for %s of Element %d is not allowed. Please check again.\n", x, s, i + 1);
        myexit(ERRNAN);
    }
}
void compareVal(double x, double y)
{
    double d = fabs(x - y);
    if (d != 0.) {
        if (d > EPS_DOUBLE) {
            //if (d / fabs(x) > 1.e-4 || d / fabs(y) > 1.e-4) {
            printf("\nValues are different\n\n");
            printf("%f \t %f \t %e\n\n", x, y, d);
            printf("\n");
        } else {
            d = d;
            //printf("\nValues are different, due to truncation error.\n\n");
            //printf("%f \t %f, error=%e\n\n", x, y, d);
        }
    }
}

double mean(double x[], int n){
    double y = x[0];
    for (int i = 1; i < n; i++) {
        y += x[i];
    }
    return y / n;
}

double stddeviation(double x[], int n){
    double ret = 0.;
    double xmean;
    xmean = mean(x, n);
    for (int i = 0; i < n; ++i) {
        ret += (x[i] - xmean) * (x[i] - xmean);
    }
    return sqrt(ret / n);
}
void CheckFile(FILE * fp, const char *s){
    if (fp == NULL) {
        fprintf(stderr, "\n  Fatal Error: \n %s is in use or does not exist!\n", s);
        myexit(ERRFileIO);
    }
}

void CheckNA(double x, const char *s){
    if (fabs(x - NA_VALUE) < EPS_DOUBLE) {
        fprintf(stderr, "\n  Fatal Error: %s is NA_value \n", s);
        myexit(ERRDATAIN);
    }
}
void CheckNA(int x, const char *s){
    if (x == NA_VALUE || fabs(x) > 1e100 ) {
        fprintf(stderr, "\n  Fatal Error: %s is NA_value \n", s);
        myexit(ERRDATAIN);
    }
}
void creatFile(const char *fn){
    creatFile(fn, "w");
}
void creatFile(const char *fn, const char *mode){
    FILE *fp = fopen(fn, mode);
    CheckFile(fp, fn);
    fclose(fp);
}
void screeninfo(const char *s){
    screeninfo(s,"");
}
void printVectorBin(FILE *fid, double * x, int xstart, int n, double t){
    fwrite (&t, sizeof (double), 1, fid);
    fwrite (x, sizeof (double), n, fid);
    fflush (fid);
}
void printVector(FILE *fid, double * x, int xstart, int n, double t){
    fprintf(fid, "%f\t", t);;
    for(int i = 0; i < n; i++){
        fprintf(fid, "%d:%.3e\t", i+1, x[i + xstart]);
    }
    fprintf(fid, "\n");;
}
