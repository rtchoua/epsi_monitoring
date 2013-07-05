#ifndef __MINMAX_H__
#define __MINMAX_H__

#include "common.h"

bool minmax_open(char *path, bool readmode); // mode: 1=read, 0=append 
void minmax_close();
void minmax_print(char *varname, int time, double min, double max);
void minmax_get(char *varname, double undef_value, double *gmin, double *gmax);

#endif
