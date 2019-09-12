#include "f.hpp"
int f(double t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS){
    double       *Y, *DY;
    Model_Data      * MD;
    MD = (Model_Data *) DS;
#ifdef _PIHMOMP
    Y = NV_DATA_OMP(CV_Y);
    DY = NV_DATA_OMP(CV_Ydot);
    MD->f_update_omp(Y, DY, t);
    MD->f_loop_omp(Y, DY, t);
    MD->f_applyDY_omp(DY, t);
#else
    Y = NV_DATA_S(CV_Y);
    DY = NV_DATA_S(CV_Ydot);
    MD->f_update(Y, DY, t);
    MD->f_loop(Y, DY, t);
    MD->f_applyDY(DY, t);
#endif
    MD->nFCall++;
    return 0;
}

