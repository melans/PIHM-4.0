//  print.hpp
//  Created by Lele Shu (lele.shu@gmail.com) on 2018.
//  Copyright © 2018 Lele Shu. All rights reserved.
//
#ifndef print_hpp
#define print_hpp
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Macros.hpp"
#include "IO.hpp"
#include "functions.hpp"
#include "Model_Data.hpp"

void initialize_output(FileOut *fout, Model_Data * DS, Control_Data  *CS);
void fun_printASCII(Print_Ctrl PCtrl, double t, double dt);
void fun_printBINARY(Print_Ctrl PCtrl, double t, double dt);
void close_output(Control_Data *CS);

void PrintDataNew (Print_Ctrl PCtrl, double tmpt, double dt, Control_Data *CS);
void PrintInit (Model_Data * DS, char *fn);
void PIHMlogo(void);

#endif /* print_hpp */
