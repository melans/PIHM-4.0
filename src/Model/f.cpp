#include "f.hpp"
extern int debug_mode;
extern double *dy;
int f(double t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS){
    double       *Y, *DY;
    Model_Data      * MD;
#ifdef _PIHMOMP
    Y = NV_DATA_OMP(CV_Y);
    DY = NV_DATA_OMP(CV_Ydot);
#else
    Y = NV_DATA_S(CV_Y);
    DY = NV_DATA_S(CV_Ydot);
#endif
    MD = (Model_Data *) DS;
    MD->f_update(Y, DY, t);
    MD->f_loop(Y, DY, t);
    MD->f_applyDY(DY, t);
    MD->nFCall++;
    return 0;
}
