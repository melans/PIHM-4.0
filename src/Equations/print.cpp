#include "print.hpp"
void PIHMlogo(void){
    printf ("\n");
    printf ("\t########  ####  ##     ##  ##     ##             #    \n");
    printf ("\t##     ##  ##   ##     ##  ###   ###             #    \n");
    printf ("\t##     ##  ##   ##     ##  #### ####     #    ####### \n");
    printf ("\t########   ##   #########  ## ### ##     #       #    \n");
    printf ("\t##         ##   ##     ##  ##     ##  #######    #    \n");
    printf ("\t##         ##   ##     ##  ##     ##     #            \n");
    printf ("\t##        ####  ##     ##  ##     ##     # Verion 4.0 \n");
    
    printf ("\n\t\tThe Penn State Integrated Hydrologic Model v4.0\n");
    
#ifdef _CALIBMODE
    printf("\nCalibration mode enable.\n");
#endif
#ifdef _DEBUG
    printf("\nModel-debug mode enable.\n");
#endif
    
#ifdef _PIHMOMP
    printf("\nopenMP enabled. Maximum Threads = %d\n", omp_get_max_threads());
#else
    printf("\nopenMP disabled.\n");
#endif
}


