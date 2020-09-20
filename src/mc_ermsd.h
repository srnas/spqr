#ifndef MCERMSD_H
#define MCERMSD_H
#define ERMSDX 5.0
#define ERMSDY 5.0
#define ERMSDZ 3.0

#define ERMSD_AU_x 2.6389 // u-a projected on a
#define ERMSD_AU_y 4.93495
#define ERMSD_AU_z -0.0446288
#define ERMSD_UA_x 2.70755
#define ERMSD_UA_y 4.89191
#define ERMSD_UA_z -0.240751
#define ERMSD_CG_x 2.76488
#define ERMSD_CG_y 4.82442
#define ERMSD_CG_z -0.251235
#define ERMSD_GC_x 2.35092
#define ERMSD_GC_y 5.04514
#define ERMSD_GC_z -0.0492806

#define LOOP_HP 0
#define LOOP_ST 1
#define LOOP_IL 2


#include "mc_global.h"
#include "mc_utils.h"

extern FILE *ermsd_obs;
extern double ERMSD_SQ;
extern double DELTA_ERMSD_SQ;
extern double ERMSD_ENERG;
extern double DELTA_ERMSD_ENERG;
extern int ERMSD_FLAG;
extern double *ERMSD_PREF;
extern double ERMSD_CUTOFF;
extern double *ermsdref_X;
extern double *ermsdref_Y;
extern double *ermsdref_Z;
extern double *ermsdref_B;
extern double *ermsdref_S;
extern double *ermsdref_P;
extern double ***G_curr;
extern double ***G_trial;
extern double ***G_ref;
extern int **G_groups;
extern int ERMSD_N_GROUPS;
extern int ERMSD_NNT;
extern double **ERMSD_SSTRUCT;

extern double *phpull_K;
extern double **phpull_p0;
extern int **loop1;
extern int **loop2;
extern int *nnt_loop1;
extern int *nnt_loop2;
extern int *in_link;
extern double Dlnksq;
extern int mc_N_links;
extern int mc_N_loops;
extern int **mc_links;
extern int **mc_loops;
extern int *mc_loop_size;
extern int *mc_loop_type;
extern double **mc_loop_CM;
extern double **mc_loop_clpair;
extern double LOOP_K_cmcm;
extern double LOOP_K_cmclp;

void MC_init_ermsd_restr(int);
void MC_init_ermsd_out(int);
void MC_copy_single_ermsd_g(int, int, int);
void MC_copy_ermsd_g(int, int);
void MC_update_ermsd_g(int, int);
double ermsd_norm(double *);
double get_ermsd();
double get_ermsd_energ();
double get_first_ermsd(double **, double **, double **, int, double*, double*);
double MC_get_pair_ermsd(double, double, double, double, double, double, double);
void build_ermsd_g(double **, double **, double **, double ***, int);
void get_ermsd_g_pair(int, int, double **, double **, double **, double ***);
void MC_assign_ermsd_g_pair(int, int, double *, double ***);
void MC_get_ermsd_pair_type(int, int, double *, int);
void MC_write_ermsd_obs(int, double);


//void MC_init_def_phpull(int);
void MC_set_linked_loops(int, double*, double *, double*);
void mc_update_loop(int, int, double*, double*, double*);
int nt_is_in_link(int, int);
int nt_is_in_loop(int, int);
int nts_in_same_link_but_different_loops(int, int);
void get_real_sugpos(int, int, double*, double*, double*, double*);
double calc_link_energy(int, double*, double*, double*);

double MC_wall_energy(double , double, double);
#endif
