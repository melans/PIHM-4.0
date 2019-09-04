//  Model_Data.hpp
//  Created by Lele Shu (lele.shu@gmail.com) on 2018.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#ifndef Model_Data_hpp
#define Model_Data_hpp

#include <stdio.h>
#include "TimeSeriesData.hpp"
#include "ModelConfigure.hpp"
#include "IO.hpp"
#include "River.hpp"
#include "Element.hpp"
#include "Model_Control.hpp"
#include "TabularData.hpp"
#include "FloodAlert.hpp"
#include "Lake.hpp"
#include "is_sm_et.hpp"
#include "Flux_RiverElement.hpp"
#include "Macros.hpp"
using namespace std;
class Model_Data {        /* Model_data definition */
public:
    double tnow;       /* Current time tag*/
    //int UnsatMode;        /* Unsat Mode */
    //int SurfMode;        /* Surface Overland Flow Mode */
    //int RivMode;        /* River Routing Mode */
//    unsigned long nFCall1;
//    unsigned long nFCall2;
    unsigned long nFCall = 0;
    double tic;
    int NumForc = 0;
    int NumEle = 0;            /* Number of Elements */
    int NumNode = 0;        /* Number of Nodes */
//    int NumY1;
//    int NumY2;
    int NumY = 0;
    int NumSSEle = 0;       /* Number of Souce/Sink for Elements */
    int NumBCEle1 = 0;      /* Number of Boundary Condition for Elements */
    int NumBCEle2 = 0;      /* Number of Boundary Condition for Elements */
    int NumBCRiv1 = 0;      /* Number of Boundary Condition for Rivers */
    int NumBCRiv2 = 0;      /* Number of Boundary Condition for Rivers */
    int NumBCLake1 = 0;      /* Number of Boundary Condition for Lakes */
    int NumBCLake2 = 0;      /* Number of Boundary Condition for Lakes */
    int NumRiv = 0;            /* Number of Rivere Reaches */
    
    
    int NumSoil = 0;        /* Number of Soils */
    int NumGeol = 0;        /* Number of Geologies */
    int NumLC = 0;            /* Number of Land Cover Index Data */
    int NumMeltF = 0;        /* Number of Melt Factor Time series */
    int NumRivType = 0;        /* Number of River Shape */
    int NumRivNode = 0;
    
    _TimeSeriesData *tsd_weather;
    _TimeSeriesData tsd_LAI;
    _TimeSeriesData tsd_RL;
    _TimeSeriesData tsd_MF;
    _TimeSeriesData tsd_eleSS; /* Element Source/Sink Term [L3/T] */
    _TimeSeriesData tsd_eyBC; /* Element Y BC */
    _TimeSeriesData tsd_eqBC; /* Element Q BC */
    _TimeSeriesData tsd_ryBC;
    _TimeSeriesData tsd_rqBC;
    _TimeSeriesData tsd_lyBC;
    _TimeSeriesData tsd_lqBC;
    int ieBC1 = 0;
    int ieBC2 = 0;
    int irBC1 = 0;
    int irBC2 = 0;
    int ilBC1 = 0;
    int ilBC2 = 0;
    int ieSS = 0;
    
    globalCal gc;
    
    _Element *Ele;        /* Store Element Information */
    _Node *Node;        /* Store Node Information */
    //element_IC * Ele_IC;    /* Store Element Initial Condtion */
    Soil_Layer *Soil;        /* Store Soil Information */
    Geol_Layer *Geol;            /* Store Geology Information */
    Landcover *LandC;        /* Store Land Cover Information */
    
    _River *Riv;        /* Store River Reach Information */
    river_para *Riv_Type;    /* Store River Shape Information */
    _Node *rivNode;
    FloodAlert *flood;
    
    double WatershedArea = 0.;
    double *ISFactor;        /* ISFactor is used to calculate ISMax from LAI */
    double *windH;        /* Height at which wind velocity is measured */
    _Lake *lake;
    int NumLake = 0;
    double *QoutSurf;
    
    double **QeleSurf;    /* Overland Flux */
    double **QeleSub;        /* Subsurface Flux */
    //double ** FluxRiv;    /* River Segment Flux */
    double *QrivSurf;        /* surface Flux between river and element */
    double *QrivSub;        /* gw Flux between river and element */
    double *QrivDown;
    double *QrivUp;
    double *QsegSurf;
    double *QsegSub;    
    
    double *QeleSurfTot;
    double *QeleSubTot;
    
    double *Qe2r_Surf;
    double *Qe2r_Sub;
    
    double *yEleWetFront;        /* Weting Front */
    
    double *qElePrep;        /* Precep. on each element */
    double *qEleETloss;
    double *qEleNetPrep;    /* Net precep. on each elment */
    double *qEleInfil;    /* Variable infiltration rate */
    double *qEleExfil;    /* Variable exfiltration rate */
    double *qEleRecharge;    /* Recharge rate to GW */
    double *yEleSnowGrnd;    /* Snow depth on ground element */
    double *yEleSnowCanopy;    /* Snow depth on canopy element */
    double *yEleISmax;    /* Maximum interception storage (liquid
                           * precep) */
    double *yEleISsnowmax;    /* Maximum interception storage (snow) */
    double *qEleTF;        /* Through Fall */
    double **qEleET;        /* Evapo-transpiration (from canopy,
                             * ground, subsurface, transpiration) */
    
    double *yEleIS;        /* Interception storage */
    double *yEleSnow;        /* Snow depth on each element */
    double *yEleGW;   // debug may not necessary
    double *yEleSurf;   // debug may not necessary
    double *yEleUnsat;   // debug may not necessary
    double *qEleETP;
    double *iBeta;
    double *iPC;
    double *qEleETA;
    double *yRivStg;   // debug may not necessary
    /* Lake variables */
    double *yLakeStg;
    double *QLakeSurf;
    double *QLakeSub;
    double *QLakeRiv;
    double *qLakeEvap;
    double *qLakePrcp;
    
//    double *uYele;
//    double *uYsf;
//    double *uYus;
//    double *uYgw;
//    double *uYriv;
    
    int NumSegmt;
    RiverSegement *RivSeg;
    
    long ForcStartTime;
    
    Control_Data CS;
private:
    double *t_prcp;
    double *t_temp;
    double *t_rh;
    double *t_wind;
    double *t_rn;
    double *t_vp;
    double *t_lai;
    double *t_mf;
    double *t_rl;
    
public:
    /* Methods: */
    Model_Data();
    ~Model_Data();
    /* Model input/output */
    void loadinput(FileIn * fin);
    void initialize();
    void initialize_output(FileOut * fout);
    void InitialCondition(FileIn *fin, N_Vector CV_Y);
    void InitialCondition(FileIn *fin, N_Vector udata1, N_Vector udata2);
    /* screen print */
    void modelSummary(FileIn * fileIn, int end);
    void PrintInit(const char *fn);
    
    void summary(N_Vector u1, N_Vector u2);
    void summary(N_Vector u);
    
    /* methods in f function */
    void f_loop(double * Y, double * DY, double t);
    void f_applyDY(double * DY, double t);
    void f_update(double * Y, double * DY, double t);
    
//    void updateWF(double dt);
    void CheckInputData();
    void InitFloodAlert(const char *fn);
    void updateRiverStage(N_Vector uY);
    void debugData();
    void debugData(const char *fn);
    void f_etFlux(int i, double t);
    void EvapoTranspiration(N_Vector uY, double t, double dt);
    void updateforcing(double t);
    double getArea();
private:
    void fillpits(int i);
    
    void tReadForcing(double t, int i);
    
    void ElementTable(const char *fn);
    void RiverTable(const char *fn);
    void LakeTable(const char *fn);
    
    /* Memory management */
    void allocateMemory();
    void FreeData();
    
    /* put calibration file into the parameters */
    void copyCalib();
    void calibSoil();
    void calibGeol();
    void calibLandc();
    
    /* Check the input data */
    void CheckInput_forc();
    void CheckInput_mesh();
    void CheckInput_att();
    void CheckInput_soil();
    void CheckInput_geol();
    void CheckInput_landcover();
    
    /* Read input data */
    void read_calib(const char *fn);
    void read_para(const char *fn);
    void read_riv(const char *fn);
    void read_rivchn(const char *fn);
    void read_mesh(const char *fn);
    
    void read_lake(const char *fn);
    void read_lakeBathy(const char *fn);
    
    void read_att(const char *fn);
    void read_soil(const char *fn);
    void read_geol(const char *fn);
    void read_lc(const char *fn);
    void read_forc_csv(const char *fn);
    void read_rl(const char *fn);
    void read_lai(const char *fn);
    void read_mf(const char *fn);
    
    void read_ssEle(const char *fn);
    void read_bcEle1(const char *fn);
    void read_bcEle2(const char *fn);
    void read_bcRiv1(const char *fn);
    void read_bcRiv2(const char *fn);
    void read_bcLake1(const char *fn);
    void read_bcLake2(const char *fn);
//    void CorrectRiver(double eps);
    void rmSinks();
    
    /* Physical processes */
    void Flux_RiverDown(double t, int i);
    void f_Segement_surface(int iEle, int iRiv, int i);
    void f_Segement_sub(int iEle, int iRiv, int i);
    void f_Segement_update(int iEle, int iRiv, int i);
    void PassValue();
    void applyBCSS(double *DY, int i);
    
    /* Methods for element calculation */
    void f_lateralFlux(int i, double t);
    void f_InfilRecharge(int i, double t);
    /* Functions */
    void TimeSpent();
    double WeirFlow(double ze, double ye, double zr, double yr,
                    double zbank, double cwr, double rivLen, double threshold);
    double updateArea();
};
#endif                /* Model_Data_hpp */

