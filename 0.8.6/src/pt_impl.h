#ifndef PTIMPL_H
#define PTIMPL_H

#include "mc_global.h"

//void SA_read_params(double *, double *, double *, int *);
void PT_set_temp(double);
void PT_create_params(double *, int);
void PT_write_params(double *, int, int);
#endif
