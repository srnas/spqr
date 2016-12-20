#ifndef CHECKPOINTS_H
#define CHECKPOINTS_H

#define DEFAULT_CHK_FREQ 50
#include "mc_global.h"
#include "mc_utils.h"

FILE *mc_configs;

extern int mc_chk_freq;

void MC_initialize_save_configs(int, int, int);
void MC_close_configs();
void MC_save_configuration(int, double *, double*, double*);
void MC_save_xyz_configuration(int, double *, double *, double *, double, int);
void MC_save_current_configuration(int, double *, double *, double *, double, int, int);
void MC_min_energ_xyz_configuration(int, double *, double *, double *, double, int );
#endif
