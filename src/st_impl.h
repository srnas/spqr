#ifndef STIMPL_H
#define STIMPL_H

#include "mc_global.h"
#include "mc_utils.h"

void ST_read_parameters(int, int, int *, double *, double *, int *, double **, double **, int *);
void ST_write_parameters(int, int, int , double , double, int, double *, double *, int );
void ST_update_temperature(double);


#endif
