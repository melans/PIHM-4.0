
#include "ObsSim.hpp"

ObsnSim::ObsnSim(){}
ObsnSim::~ObsnSim(){
    delete[] time;
    delete[] sim;
    delete[] obs;
    delete gof;
}
void ObsnSim::init(int n){
    time = new double[n];
    sim = new double[n];
    obs = new double[n];
    for(int i = 0; i < n; i++){
        time[i] = 0.;
        sim[i] = 0.;
        obs[i] = 0.;
    }
}
void ObsnSim::readobs(const char *fn){
    int i;
    FILE *fp = fopen(fn, "r");
    CheckFile(fp, fn);
    char str[MAXLEN];
    fgets(str, MAXLEN, fp); // Dimension of the table
    sscanf(str, "%d %d %d", &nTotal, &nSpin, &targetID);
    targetID = targetID - 1;
    if(targetID < 0 || targetID > 100000){
        fprintf(stderr, "ERROR: Index of simulation target is %d.\n",targetID);
        myexit(ERRFileIO);
    }
    init(nTotal);
    SimMean = 0.;
    fgets(str, MAXLEN, fp); // the header
    for (i=0; i < nTotal && fgets(str, MAXLEN, fp) ; i++){
        sscanf(str, "%lf %lf", &(time[i]), &(obs[i]) );
//        printf("%g\t%g\n", time[i], obs[i]);
    }
    if( i < nTotal){
        fprintf(stderr, "Wrong length of observation data. %d vs %d\n", i+1, nTotal);
        myexit(ERRDATAIN);
    }
    DayStart = time[0];
    DayEnd = time[nTotal - 1];
    nCalib = nTotal - nSpin;
    fclose(fp);
}
void ObsnSim::Pointer2Sim(double *x, double ct1, double ct2){
    ptr = &x[targetID];
    ctm = ct1;
    cta = ct2;
}
void ObsnSim::pushsim(double t){
//    printf("%ld\n", ptr);
    double tday = t / UNIT_C;
    static int nr = 0;
    nr ++;
    sim[iNow] += *ptr * ctm + cta;
    
    if( time[iNow] + 1 <= tday ){
//        printf("Ntimes = %d\n" , nr);
        sim[iNow] = sim[iNow] / nr;
        iNow ++;
        if(iNow < nTotal){
            sim[iNow] = 0.0; /* Initial value*/
        }
        nr = 0;
    }else{
    }
    if( iNow > nTotal ){
        fprintf(stdout, "End of calibration period\n");
    }
}
double ObsnSim::getNSE(){
    double ret = 0.;
    callGOF();
    ret = gof->gof_NSE();
    return ret;
}
void ObsnSim::callGOF(){
    gof = new GoodOfFit;
    gof->init( &(sim[nSpin]), &(obs[nSpin]), nTotal - nSpin);
}
void ObsnSim::printGOF(FILE *fp){
    gof->print_gofname(fp);
    gof->print_gof(fp);
}
void ObsnSim::printData(const char *fn){
    FILE *fp = fopen(fn, "w");
    CheckFile(fp, fn);
    
    fprintf(fp, "%s\t%s\t%s\n", "TIME", "Observation", "Simulation");
    for(int i = 0; i < nTotal; i++){
        fprintf(fp, "%g\t%g\t%g\n", time[i], obs[i], sim[i]);
    }
    fclose(fp);
}
