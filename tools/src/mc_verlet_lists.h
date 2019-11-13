#ifndef VLISTS_H
#define VLISTS_H
#include "mc_global.h"
#include "mc_energies.h"
/* this has to be consistent with the dimensionality of the system */
#define N_OFFSET 14


extern int *vl_n_pairs;
extern int vl_n_pairs_max;
//int vl_flag;
extern int vl_neighbors_per_atom;
extern double vl_max_displ;

extern double vl_skin;
extern double vl_rmax;
extern double vl_rmax_sq;
extern double vl_max_vel;

extern int **vl_neighbor_tables;
extern int *vl_count;
extern int vl_ncrit;

void MC_add_vl_count(int);
void MC_initialize_verlet_lists(int, double *, double *, double *);
void MC_build_verlet_lists(int, double *, double *, double *, int);
void MC_free_verlet_lists();
//void MC_boundary_condition(int *, int *, int *);
void MC_set_adjacent_cell(int*, int*, int);
#ifdef EXCLUSIONS
int MC_are_neighbors(int, int);
#endif

#endif
