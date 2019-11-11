#ifndef SAIMPL_H
#define SAIMPL_H

#include "mc_global.h"
#include "mc_utils.h"
void SA_read_params(double *, double *, double *, int *, int *, double*, double*, int*,int);
void SA_set_temp(double);
void SA_set_mc_trials(double, double, double);
void SA_rescale_mc_trials(double *, double *, double *, double);
void SA_init_mc_trials(double *, double *, double *, double, double , double);
void SA_write_params(double, double, double, int, int, double, double, int, int);
double sa_pow(double , int);
void SA_params_to_arr(double, double, double, int, int, double, double, int, double *);
void SA_arr_to_params(double *, double *, double *, int *, int *, double*, double*, int*,double *);
#endif
