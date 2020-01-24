#ifndef f_hpp
#define f_hpp
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ModelConfigure.hpp"
#include "Macros.hpp"
#include "functions.hpp"
#include "Model_Data.hpp"

int f(double t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS);
#endif /* f_h */
