#ifndef MCCELLS_H
#define MCCELLS_H
#include "mc_global.h"

extern int *mc_cells;
extern int mc_n_linked_cells;
extern int mc_nc[DIM];
extern double mc_linked_cell_l;
void add_particle_to_mc_cell(int, int, int);
void remove_particle_from_mc_cell(int);
void MC_update_linked_lists(int, double *, double *, double *);
void calculate_cell_index(int *, double, double, double);
#endif
