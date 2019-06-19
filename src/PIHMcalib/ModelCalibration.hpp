//
//  CalibDataset.hpp
//  PIHMcalib
//
//  Created by Lele Shu on 9/19/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#ifndef ModelCalibration_hpp
#define ModelCalibration_hpp

#include <stdio.h>
#include <string>
#include "ModelConfigure.hpp"
#include "IO.hpp"

using namespace std;
class ModelCalibration{
public:
    int N = 0; /* Number of calib parameters*/
    int lamda = 0; /* Number of samples */
    int nBest = 0;
    globalCal cMax;
    globalCal cMin;
    globalCal *gc;
    globalCal gc0;
    globalCal gcbest;
    string  *fns_calib;      /* Calibfiles of each PIHM simulations*/
    string  *outdirs; /* outdirs of each PIHM simulations*/
    char  cmaes_dir[MAXLEN];
    char  fn_config[MAXLEN];
    char  fn_range[MAXLEN];
    char  fn_cinit[MAXLEN];
    char  fn_csignal[MAXLEN];
    char  fn_result[MAXLEN];
    string vn[100];
    int logscale[100];
    int nRound = 0;
    
    ModelCalibration();
    ~ModelCalibration();
    void init(int nlambda, const char *outdir, const char *cbfile, const char *cmaesdir);
    void setCalibration(double * const *pop);
    void setCalibration(const char *var,  double sig, int igc);
    void backupCalibdata();
    void writefiles();
    void outBestData(int igc);
    void setIO(FileIn *fin, FileOut *fout, int index);
    void printObjValue(double ov, int index);
    
    void readRange();
    void readConfig();
private:
    string outpath;
    string file_calibdata;
    string file_popdata;
    string file_final;
    string file_steps;
    string file_ov;
    double **xrange;
    double ConvertRange(const char *varname, double sig, int ilog);
    void Range(double *xmin, double *xmax);
    void setFiles(const char *str,const char *dir);
};
#endif /* CalibDataset_hpp */
