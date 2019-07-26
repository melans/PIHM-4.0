
#include "Model_Data.hpp"
#include "is_sm_et.hpp"

Model_Data::Model_Data(){
}

Model_Data::~Model_Data(){
    FreeData();
}
void Model_Data::TimeSpent(){
#ifdef _PIHMOMP
    double toc = omp_get_wtime();
    double dt = toc - tic;
    screeninfo("\n\tNumber of calls of f function:\t %ld \n", nFCall);
    printf("\n\tTime used by model:\t %.3f seconds.\n", dt);
    screeninfo("\n\nThe successful end. \n\n");
    
#else
    clock_t toc = (double)clock();
    double dt = (toc - tic) / CLOCKS_PER_SEC;
    screeninfo("\n\tNumber of calls of f function:\t %ld \n", nFCall);
//    screeninfo("\n\tNumber of calls of f function:\t %ld \n", nFCall2);
    screeninfo("\n\tTime used by model:\t %.3f seconds.\n", dt);
    screeninfo("\n\nThe successful end. \n\n");
#endif
    
}
void Model_Data::modelSummary(FileIn * fileIn, int end){
    char str[MAXLEN];
    screeninfo("\n========================================================\n");
    screeninfo("Summary:\n");
    screeninfo("\tProject name:\t %s\n ", fileIn->projectname);
    screeninfo("\tInput path:\t %s\n", fileIn->inpath);
    screeninfo("\tOutput path:\t %s\n", fileIn->outpath);
    screeninfo("\tCalibration file:\t %s\n", fileIn->file_calib);
    screeninfo("\tParameter file:\t %s\n", fileIn->file_para);
    screeninfo("\tModel starts at: %.2f day\n", CS.StartTime / 1440);
    screeninfo("\tModel ends at: %.2f day\n", CS.EndTime / 1440);
    screeninfo("\tModel time step(max): %.2f minutes\n", CS.MaxStep);
    screeninfo("\tModel total number of steps(minimum): %d \n", CS.NumSteps);
    sprintf(str,"\tSize of model: \tNcell = %d \tNriver = %d\t NSeg = %d", NumEle, NumRiv, NumSegmt);
    screeninfo(str);
#ifdef _PIHMOMP
    screeninfo("\n\n\tOpenMP enable. No of threads = %d\n", CS.num_threads);
    screeninfo("\n========================================================\n");
    if (end) {
        TimeSpent();
    } else {
        tic = omp_get_wtime();
        nFCall = 0;
        screeninfo("\nModel Starting ... \n\n");
    }
#else
    screeninfo("\n\tOpenMP disable");
    screeninfo("\n========================================================\n");
    if (end) {
        TimeSpent();
    } else {
        tic = (double)clock();
        nFCall = 0;
//        nFCall2 = 0;
        screeninfo("\nModel Starting ... \n\n");
    }
#endif
}

void Model_Data::allocateMemory()
{
    /* allocate memory storage to flux terms */
    QeleSurf    = new double *[NumEle];
    QeleSub     = new double *[NumEle];
    QeleSurfTot = new double[NumEle];
    QeleSubTot  = new double[NumEle];
    QoutSurf    = new double[NumEle]; // 5
    
    Qe2r_Surf = new double[NumEle]; //5.1
    Qe2r_Sub  = new double[NumEle]; // 5.2
    
    qEleET      = new double *[NumEle];
    qElePrep    = new double[NumEle];
    qEleTF      = new double[NumEle];
    qEleETP     = new double[NumEle];
    qEleETA     = new double[NumEle];
    qEleETloss  = new double[NumEle]; //10
    
    qEleNetPrep = new double[NumEle];
    qEleInfil   = new double[NumEle];
    qEleRecharge = new double[NumEle]; //13
    
    yEleIS      = new double[NumEle];
    yEleISmax   = new double[NumEle];
    yEleISsnowmax = new double[NumEle]; //16
    
    yEleSnow    = new double[NumEle];
    yEleSnowGrnd = new double[NumEle];
    yEleSnowCanopy = new double[NumEle];
    yEleGW      = new double[NumEle];    //20
    
    yEleSurf    = new double[NumEle];
    yEleUnsat   = new double[NumEle];
    yEleWetFront = new double[NumEle];  //23
    
    yRivStg     = new double[NumRiv];
    QrivSurf    = new double[NumRiv];
    QrivSub     = new double[NumRiv];   //26
    QrivDown    = new double[NumRiv];
    QrivUp      = new double[NumRiv]; //28
    QsegSurf    = new double[NumSegmt];
    QsegSub     = new double[NumSegmt];
    
    if (NumLake > 0){
        yLakeStg    = new double[NumLake];
        QLakeSurf   = new double[NumLake];  //30
        QLakeSub    = new double[NumLake];
        QLakeRiv    = new double[NumLake];
        qLakePrcp   = new double[NumLake];
        qLakeEvap   = new double[NumLake];  //34
    }
    
//    NumY1 = 3 * NumEle;
//    NumY2 = 1 * NumRiv + 1 * NumLake;
    NumY = 3 * NumEle + 1 * NumRiv + 1 * NumLake;
    
    uYsf = new double[NumEle];
    uYus = new double[NumEle];
    uYgw = new double[NumEle];
//    uYele = new double[NumY1];  // 35
    if(NumRiv > 0){
        uYriv = new double[NumRiv];  // 35.1
    }
    
    for (int i = 0; i < NumEle; i++) {
        QeleSurf[i] = new double[3];
        QeleSub[i] = new double[3];
        qEleET[i] = new double[3];
    }
    
    t_prcp  = new double[NumEle];  //
    t_temp  = new double[NumEle];  //
    t_rh    = new double[NumEle];  //
    t_wind  = new double[NumEle];  //
    t_rn    = new double[NumEle];  //
    t_vp    = new double[NumEle];  //
    t_lai   = new double[NumEle];  //
    t_mf    = new double[NumEle];  //
    t_rl    = new double[NumEle];  //
}

void Model_Data::copyCalib(){
    for (int i = 0; i < NumSoil; i++) {
        Soil[i].applyCalib(&(gc.csoil));
    }
    for (int i = 0; i < NumGeol; i++) {
        Geol[i].applyCalib(&(gc.cgeol));
    }
    for (int i = 0; i < NumLC; i++) {
        LandC[i].applyCalib(&(gc.clandc));
    }
    for (int i = 0; i < NumRivType; i++) {
        Riv_Type[i].applyCalib(&(gc.criv));
    }
}
void Model_Data::InitFloodAlert(const char *fn){
    flood->InitAlert(NumRiv, NumRivType);
    flood->InitPointer(yRivStg, QrivDown);
    flood->InitPara(Riv_Type);
    for(int i =  0; i < NumRiv; i++){
        flood->pushRiverType(i, Riv[i].type);
    }
    flood->InitFile(fn);
}
double Model_Data::updateArea(){
    WatershedArea = 0.;
    for(int i = 0; i < NumEle; i++){
        WatershedArea += Ele[i].area;
    }
//    for(int i = 0; i < NumLake; i++){
//        WatershedArea += Lake[i].area;
//    }
    return WatershedArea;
}
double Model_Data::getArea(){
    return WatershedArea;
}
void Model_Data::rmSinks(){
    double zmin, z;
    int inabr;
    for(int i = 0; i < NumEle; i++){
        z =  Ele[i].zmax;
        zmin = 1.0e200;
        for(int j = 0; j < 3; j++){
            inabr = Ele[i].nabr[j] - 1;
            if(inabr >= 0){ /* Nabr exists */
                zmin = min(zmin, Ele[inabr].zmax);
            }
        }
        if( zmin > z + Ele[i].AquiferDepth){
            fprintf(stderr, "Warning: remove sink on %d, from %.2f to %.2f. dz = %.2f\n", i+1, z, zmin, zmin - z);
            Ele[i].zmax = zmin;
            Ele[i].zmin = zmin - Ele[i].AquiferDepth;
        }
    }
    
    for (int i = 0; i < NumEle; i++) {
        Ele[i].InitElement();
    }
}

void Model_Data::debugData(const char *outdir){
    char fn[MAXLEN];
    char str[MAXLEN];
    sprintf(str, "%s/%s", outdir, "Debug_Table");
    if(NumEle > 0){
        sprintf(fn, "%s%s", str, "_Element.csv");
        ElementTable(fn);
    }
    if(NumRiv > 0){
        sprintf(fn, "%s%s", str, "_River.csv");
        RiverTable(fn);
    }
    if(NumLake > 0){
        sprintf(fn, "%s%s", str, "_Lake.csv");
        LakeTable(fn);
    }
}
void Model_Data::ElementTable(const char *fn){
    FILE *fp = fopen(fn, "w");
    Ele[0].printHeader(fp);
    for(int i = 0; i < NumEle; i++){
        Ele[i].printInfo(fp);
    }
    fclose(fp);
}
void Model_Data::RiverTable(const char *fn){
    FILE *fp = fopen(fn, "w");
    Riv[0].printHeader(fp);
    for(int i = 0; i < NumRiv; i++){
        Riv[i].printInfo(fp);
    }
    fclose(fp);
}
void Model_Data::LakeTable(const char *fn){
    FILE *fp = fopen(fn, "w");
    fclose(fp);
}

