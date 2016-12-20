#ifndef SAIMPL_H
#define SAIMPL_H

#include "mc_global.h"

void SA_read_params(double *, double *, double *, int *, int *, double*, double*, int*,int);
void SA_set_temp(double);
void SA_set_mc_trials(double, double, double);
void SA_rescale_mc_trials(double *, double *, double *, double);
void SA_init_mc_trials(double *, double *, double *, double, double , double);
void SA_write_params(double, double, double, int, int, double, double, int, int);
double sa_pow(double , int);

#endif
