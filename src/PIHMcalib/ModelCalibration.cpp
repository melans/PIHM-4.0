#include "ModelCalibration.hpp"
ModelCalibration::ModelCalibration(){
    
}
ModelCalibration::~ModelCalibration(){
    delete[] fns_calib;
    delete[] outdirs;
    delete[] gc;
//    for(int i = 0; i < N; i++){
//        delete[] xrange[i];
//    }
//    delete xrange;
}
void ModelCalibration::init(int nsample, const char *outdir, const char *cbfile, const char *cdir){
    strcpy(cmaes_dir, cdir);
    sprintf(fn_config, "%s/cmaes.cfg", cmaes_dir);
    sprintf(fn_range, "%s/range.txt", cmaes_dir);
    sprintf(fn_cinit, "%s/initials.par", cmaes_dir);
    sprintf(fn_csignal, "%s/signals.par", cmaes_dir);
    sprintf(fn_result, "%s/allcmaes.txt", cmaes_dir);
//    readConfig(); /* read the configuration. esp. lambda and sigma */
    readRange();
    lamda = nsample;
    if(lamda < 5){
        fprintf(stderr, "WARNING: lambda value is too small %d\n\n", lamda);
    }
    fns_calib = new string[lamda];
    outdirs = new string[lamda];
    gc = new globalCal[lamda];
    for(int i = 0; i < lamda; i++){
        gc[i] = gc0;
    }
    
    outpath = outdir;
    file_calibdata = outpath + "/" + "cmaes_DataSets.csv";
    file_popdata = outpath + "/" + "cmaes_PopData.csv";
    file_final = outpath + "/" + "Final.calib";
    file_ov = outpath + "/" + "cmaes_ObjVal.csv";
    nBest = 0;
    creatFile(file_calibdata.c_str());
    creatFile(file_ov.c_str());
    
    setFiles(cbfile, outpath.c_str());
    
    
    FILE *fp = fopen(file_popdata.c_str(), "w");
    CheckFile(fp, file_popdata.c_str());
    fprintf(fp, "%s\t%s\t", "Round", "Process");
    for(int j = 0; j < N; j++){
        fprintf(fp, "%s\t",  vn[j].c_str() ) ;
    }
    fprintf(fp, "\n");
    fclose(fp);
    
    fp = fopen(file_calibdata.c_str(), "w");
    CheckFile(fp, file_calibdata.c_str());
    fprintf(fp, "%s\t%s\t", "Round", "Process");
    for(int j = 0; j < N; j++){
        fprintf(fp, "%s\t",  vn[j].c_str() ) ;
        
    }
    fprintf(fp, "\n");
    fclose(fp);
}
void ModelCalibration::outBestData(int igc){
    gcbest.copy(&(gc[igc]) ); // update to gcbest
    gc0.copy(&(gc[igc])  ); // update to gc0, which will be use for next generation POP.
    nBest++;
    file_steps = outpath + "/cmaes_step" + to_string(nBest) + ".calib";
//    gc0.write(file_steps.c_str());
    gcbest.write(file_steps.c_str());
}

void ModelCalibration::readRange(){
    FILE *fp;
    char str[MAXLEN], optstr[MAXLEN];
    fp = fopen (fn_range, "r");
    CheckFile(fp, fn_range);
    /* start reading calib_file */
    double x[2];
    N = 0;
    while (fgets(str, MAXLEN, fp)) {
        if (str[0] == '#' || str[0] == '\n' || str[0] == '\0' || str[0] == ' '){
            continue;
        }
        sscanf (str, "%s %lf %lf %d", optstr, x, x+1, &(logscale[N]) );
        cMin.push(optstr, x[0]);
        cMax.push(optstr, x[1]);
        vn[N] = optstr;
        N++;
    } // end of while gets()
    fclose (fp);
}
double ModelCalibration::ConvertRange(const char *var, double xsig, int ilog){
    double xmin, xmax, xlim;
    double ret;
    xmin = cMin.getValue(var);
    xmax = cMax.getValue(var);
    double x0 = gc0.getValue(var);
    double sig = xsig;
    if (ilog){
        /* Parameters range in log scale. Such as, Ksat ranges from
         0.001 to 1000, i.e. 10^-3 ~ 10^+3 */
        if(sig > 0){
            xlim = fabs( log10(xmax) - log10(x0) );
        }else{
            xlim = fabs(log10(x0) - log10(xmin) );
        }
        //        xlim = (log10(xmax) - log10(xmin) ) * 0.5;
        ret = pow(10, log10(x0) + sig * xlim );
        ret = ret;
    }else{
        if(sig > 0){
            xlim = fabs(xmax - x0);
        }else{
            xlim = fabs(x0 - xmin);
        }
        //        xlim = (xmax - xmin ) * 0.5;
        ret = x0 + sig * xlim;
        ret = ret;
    }
    if(ret < xmin){
        ret = xmin;
    }
    if(ret > xmax){
        ret = xmax;
    }
    return ret;
}
void ModelCalibration::writefiles(){
    for(int i = 0; i < lamda; i++) {
//        printf("%s\n", fns_calib[i].c_str());
        gc[i].write(fns_calib[i].c_str());
    }
}
void ModelCalibration::setIO(FileIn *fin, FileOut *fout, int index){
    fin->setCalibFile(fns_calib[index].c_str());
    fin->setOutpath(outdirs[index].c_str());
    fout->setOutpath(outdirs[index].c_str());
}
void ModelCalibration::printObjValue(double ov, int index){
    FILE *fp = fopen(file_ov.c_str(), "a+");
    CheckFile(fp, file_ov.c_str());
    if(nRound <= 1){
        fprintf(fp, "Round\tProcess\tObjVal\n");
    }
    nRound++;
    fprintf(fp, "%d\t%d\t%f\n", nRound, index + 1, ov);
    fclose(fp);
}
void ModelCalibration::setCalibration(double *const *pop){
    double y;
    for(int j = 0; j < N ; j++){
        printf("%-10s: ", vn[j].c_str());
        for(int i = 0; i < lamda; i++) {
            y = ConvertRange(vn[j].c_str(), pop[i][j], logscale[j]);
            gc[i].push(vn[j].c_str(), y);
//            printf("%s: %.4f(%.4f)\t", vn[j].c_str(), y, pop[i][j]);
            printf("%7.2f(%5.2f)\t", y, pop[i][j] - 0.5);
        }
        printf("\n");
    }
    writefiles();
    backupCalibdata();
    FILE *fp = fopen(file_popdata.c_str(), "a+");
    for(int i = 0; i < lamda; i++){
        fprintf(fp, "%d\t%d\t", nRound, i + 1);
        for(int j = 0; j < N; j++){
            fprintf(fp, "%f\t",  pop[i][j] );
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}
void ModelCalibration::setFiles(const char *str, const char *dir){
    for(int i = 0; i < lamda; i++){
        fns_calib[i] = str + to_string(i);
        outdirs[i] = dir; // + to_string(i);
    }
}
void ModelCalibration::backupCalibdata(){
    FILE *fp = fopen(file_calibdata.c_str(), "a+");
    for(int i = 0; i < lamda; i++){
        fprintf(fp, "%d\t%d\t", nRound, i + 1);
        for(int j = 0; j < N; j++){
            fprintf(fp, "%f\t",  gc[i].getValue(vn[j].c_str()) );
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}
//void ModelCalibration::readConfig(){
//    FILE *fp;
//    char str[MAXLEN], optstr[MAXLEN];
//    double val;
//    fp = fopen (fn_config.c_str(), "r");
//    CheckFile(fp, fn_config.c_str());
//    /* start reading calib_file */
//    N = 0;
//    while (fgets(str, MAXLEN, fp)) {
//        if (str[0] == '#' || str[0] == '\n' || str[0] == '\0' || str[0] == ' '){
//            continue;
//        }
//        //        sscanf (str, "%s %lf %lf %lf", optstr, x, x+1, x+2 );
//        sscanf (str, "%s %lf", optstr, &val);
//        if (strcasecmp ("lambda", optstr) == 0)
//            lamda = (int) val;
//        else if (strcasecmp ("sigma", optstr) == 0)
//            sigma = val;
//        /* Unrecognized Parameter Flag */
//        else{
//            printf
//            ("\n  Parameter: %s cannot be recognized. Please see User's Manual for more details!\n",
//             optstr);
//            exit (1);
//        }//end ifelse
//    } // end of while gets()
//    fclose (fp);
//}
