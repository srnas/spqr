#ifndef CHECKPOINTS_H
#define CHECKPOINTS_H

#define DEFAULT_CHKP_STEPS 1000
#define DEFAULT_TRAJ_STEPS 100
#include "mc_global.h"
#include "mc_energies.h"
#include "mc_utils.h"
#include "mc_initialization.h"
#include <sys/stat.h>
FILE *mc_configs;
FILE *pdb_configs;

extern int mc_chkp_steps;
extern int mc_traj_steps;
extern int PDB_OUTPUT;

void MC_initialize_save_configs(int, int);
void MC_write_pdb(char*, int, double*, double*, double*, double, int);
void MC_close_configs();
void MC_save_configuration(int, double *, double*, double*, double,int);
void MC_save_xyz_configuration(int, double *, double *, double *, double, double, int);
void MC_save_current_configuration(int, double *, double *, double *, double, int, double, int);
void MC_min_energ_xyz_configuration(int, double *, double *, double *, double, double, int );
int MC_read_checkpoint(int *, double **, double **, double **, int *, int, char *, int, double *);
void MC_save_checkpoint(int, double *, double *, double *, int, double, int, double *);
void MC_append_pdb(int, double *, double *, double *, double,int);
#endif
