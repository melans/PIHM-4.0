#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "ForcingData.hpp"

void CheckFile(std::ifstream * fp, const char *s)
{
    if (fp == NULL) {
        fprintf(stderr, "\n  Fatal Error: \n %s is in use or does not exist!\n", s);
        myexit(ERRFileIO);
    }
}

_DataQueue::_DataQueue(){
}

_DataQueue::~_DataQueue(){
    for (int i = 0; i < MAXQUE; i++) {
        delete[] ts[i];
    }
#ifdef _DEBUG
    printf("TSD to %s Destructed.\n", fn.c_str());
#endif
}

void _DataQueue::initialize(int n)
{
    ncol = n; /* Include the Time and Data columns */
    
    iNow = 0;
    nQue = 0;
    eof = 0;
    for (int i = 0; i < MAXQUE + 1; i++) {
        pRing[i] = i + 1;
    }
    pRing[MAXQUE] = 0;
    //for repeating last value;
    iNow = 0;
    iNext = pRing[iNow];
    
    for (int i = 0; i < MAXQUE + 1; i++) {
        ts[i] = new double[n];
    }
}
void _DataQueue::read_csv()
{
    if (!eof) {
        std::ifstream file(fn);
        CheckFile(&file, fn.c_str());
        std::string str;
        Length = 0;
#ifdef _DEBUG
        if (nQue > 0) {
            std::cout << "No of Queue = " << nQue << std::endl;;
        }
#endif
        for (int i = 0; i < MAXQUE * nQue + 2; i++) {
            /* Line 1= size of table Line 2= Head of table */
            getline(file, str);
        }
        for (int i = 0; i < MAXQUE && getline(file, str); i++) {
            //std::cout << str << endl;
            std::istringstream iss(str);
            Length++;
            for (int j = 0; j < ncol; j++) {
                iss >> ts[i][j];
                if (j == 0) {
                    ts[i][j] *= 1440.;    /* Day to minute */
                }
                //cout << ts[i][j] << "\t";
            }
            //cout << endl;
        }
        nQue++;            /* Number of the Queue was reload */
        if (!file.eof()) {
            eof = 0;
        } else {
            eof = 1;
        }
#ifdef _DEBUG
        std::cout << fn << "\n\tUpdate queue. Length = " << Length << std::endl;
#endif
        if(Length <= 0){
            fprintf(stderr, "Reading fail, file = %s\n", fn.c_str());
            myexit(ERRFileIO);
        }
        file.close();
    }
}
void _DataQueue::readDimensions()
{
    int tmp, nc;
    FILE *fp = fopen(fn.c_str(), "r");
    CheckFile(fp, fn.c_str());
    fscanf(fp, "%d %d %ld", &tmp, &nc, &StartTime);
#ifdef _DEBUG
    fprintf(stdout, "Header of %s : %d\t%d\t%ld ",fn.c_str(),  tmp, nc, StartTime);
#endif
    fclose(fp);
    initialize(nc);
}
double _DataQueue::getX(double t, int col)
{
    return ts[iNow][col];
}
int _DataQueue::get_Ncol(){
    return ncol;
}
void _DataQueue::applyCalib(double prcp, double temp)
{
    for (int i = 0; i < Length; i++) {
        ts[i][0] *= prcp;    /* Calibration of prcp */
        ts[i][1] += temp;    /* Calibration of temp */
    }
}
void _DataQueue::movePointer(double t)
{
    
    while (t >= ts[iNext][0] && ts[iNext][0] >= ts[iNow][0]) {
        if (iNow == 0) {
            for (int i = 0; i < ncol; i++) {
                ts[Length][i] = ts[Length - 1][i];
                //repeat last item;
            }
        }
        iNow = iNext;
        iNext = pRing[iNow];
        if (iNext == 0) {
            read_csv();
        }
        //printf("%s,  %.1f \t [%.4f] \t %.1f\n", fn.c_str(), ts[iNow][0], t, ts[iNext][0]);
    }
    if (ts[iNext][0] < ts[iNow][0] && t - ts[iNow][0] > 1) {
        fprintf(stderr, "\n\nError: missing forcing data after t=%.3lf\n\n", ts[iNow][0] / 1440.);
        myexit(ERRFileIO);
    }
}

void _DataQueue::checkValue(int icol, double xmin, double xmax, const char *varname){
    for(int i = 0; i < MAXQUE & i < Length; i++){
        if( ts[i][icol] < xmin || ts[i][icol] > xmax){
            fprintf(stderr, "Warning: value of %s(t=%f day) = %f is out of range (%f, %f).\n", varname, ts[i][0], ts[i][icol], xmin, xmax);
        }
    }
}
