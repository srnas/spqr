#ifndef MCINIT_H
#define MCINIT_H
#include "mc_global.h"
#include "mc_cells.h"
#include "mc_energies.h"
#include "mc_verlet_lists.h"
#include "mc_bond_lists.h"
#include "mc_integrate.h"
#include "mc_utils.h"
#include "mc_checkpoints.h"
#include "mc_base_pairing.h"
//#include "../mpc_interface.h"

#define MAX_BUFFER 100

#define READ_MC_CONF 1
#define DONT_READ_MC_CONF 0
#define READ_MC_CHECK 2

void MC_initialize(int, double **, double **, double **, double, double, double, int, int, int, int);
void MC_initialize_global(int, double, double, double, int, int);
void MC_print_parameters(int);
void MC_read_nsolute(int *, int);
void MC_initialize_positions(int, double **, double **, double **, int, int);
void MC_initialize_energies(int, double *, double *, double *);
void MC_print_positions(int, double **, double **, double **);
void MC_print_params(double, double, double, int, int);
void MC_read_params(double *, double *, double *, int *, int *, int);
void MC_initialize_sugars(int, double **, double **, double **);
#ifdef FROZEN
void MC_init_frozen_file(int);
#endif
#ifdef SECONDSTC
void MC_init_sec_str_c_lists(int);
#endif

#endif
