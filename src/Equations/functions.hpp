#ifndef functions_h
#define functions_h
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "Macros.hpp"

/* Time Interpolation */
double timeInterp(double t, double t0, double t1, double x0, double x1);

/*==========Screen print function===============*/
int    ScreenPrint(double t, int dt, unsigned long it, unsigned long nt, unsigned long nfcall);
//int ScreenPrint(double t, int dt, unsigned long it, unsigned long nt, unsigned long  nf1, unsigned long  nf2);
int    atInterval(double x, int intv);
double getSecond(void);

/*==========Math function===============*/
double Eudist(double x1, double y1, double x2, double y2);
double min(double a, double b);
double max(double a, double b);
double min(double x[], int n);
double max(double x[], int n);

double mean(double *x, int n);
double stddeviation(double x[], int n);
//double polygonArea(double X[], double Y[], int n);/* Not ready yet.*/

/*==========misc function===============*/
void CheckNonZero(double x, int i, const char *s);
void CheckNonZero(int x, int i, const char *s);
void CheckNonZero(double x, int i, int j, const char *s);
void CheckNANi(double x, int i, const char *s);
void CheckNANij(double x, int i, const char *s);
void CheckNAvalue(double x, const char *s);
void CheckNA(int x, const char *s);
void creatFile(const char *fn);
void creatFile(const char *fn, const char *mode);
void compareVal(double x, double y);
void CheckFile(FILE * fp, const char *s);

int checkRange(int x, int xmin, int xmax, int i, const char *s);
int checkRange(double x, double xmin, double xmax, int i, const char *s);

void myexit(int flag);
void screeninfo(const char *s);

template<typename  T>
void screeninfo(const char *s, T x){
#ifndef _CALIBMODE
    char str[MAXLEN];
    sprintf(str, s, x);
    fprintf(stdout, "%s", str);
#endif
}

double Quadratic(double a, double b, double c);
double fun_dAtodY(double dA, double w0, double s);


inline double Eudist(double x1, double y1, double x2, double y2){
    //Euclidean distance calculator
    double dx = x2 - x1;
    double dy = y2 - y1;
    double d = sqrt(dx * dx + dy * dy);
    return d;
}
inline double min(double a, double b){
    return (a > b ? b : a);
}
inline double max(double a, double b){
    return (a < b ? b : a);
}

inline double Quadratic(double a, double b, double c){
    double ret = 0.;
    double cc = b * b - 4 * a * c;
    if (cc > 0.) {
        ret = (-b + sqrt(cc)) / (2*a);
    }else if (cc == 0) {
        ret = -1. * b / ( 2. * a );
    }else{
        fprintf(stderr, "Error in Quadratic\n");
    }
    return ret;
}
inline double fun_dAtodY(double dA, double w0, double s){
    /* Delta_A_cs of flux to Dy in river stage @t=t+1*/
    /* w0 = topwidth at t */
    double dy = 0.;
    if(dA == 0.) return 0.;
    
    if( fabs(s) < EPS_SLOPE ){
        dy = dA / w0;
    }else{
        if(dA > 0){
            dy = Quadratic(s, w0, -1. * dA);
        }else{
            dy = Quadratic(-s, w0, -1. * dA);
        }
    }
    return dy;
}

inline
double dhdx(double *x, double *y, double *h){
    return -1. * (y[2] * (h[1] - h[0]) +
                  y[1] * (h[0] - h[2]) +
                  y[0] * (h[2] - h[1])) /
    (x[2] * (y[1] - y[0]) +
     x[1] * (y[0] - y[2]) +
     x[0] * (y[2] - y[1]));
}
inline
double dhdy(double *x, double *y, double *h){
    return -1. * (x[2] * (h[1] - h[0]) +
                  x[1] * (h[0] - h[2]) +
                  x[0] * (h[2] - h[1])) /
    (y[2] * (x[1] - x[0]) +
     y[1] * (x[0] - x[2]) +
     y[0] * (x[2] - x[1]));
}
#endif /* functions_h */



