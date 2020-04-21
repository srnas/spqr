#ifndef MCINTEGRATE_H
#define MCINTEGRATE_H
#include "mc_global.h"
#include "mc_energies.h"
#include "mc_verlet_lists.h"
#include "mc_utils.h"
#include "mc_base_pairing.h"
#include "mc_ermsd.h"


#define MAP_SUG_GL_X_A_R -1.56584
#define MAP_SUG_GL_Y_A_R -4.48467 
#define MAP_SUG_GL_Z_A_R 0.709821
#define MAP_SUG_GL_X_H_R -1.56584
#define MAP_SUG_GL_Y_H_R -4.48467 
#define MAP_SUG_GL_Z_H_R 0.709821
/* #define MAP_SUG_GL_X_H_R -1.45095  */
/* #define MAP_SUG_GL_Y_H_R -4.72806 */
/* #define MAP_SUG_GL_Z_H_R -0.2771 */
#define MAP_SUG_GL_X_S_R -0.0632599
#define MAP_SUG_GL_Y_S_R -4.34523
#define MAP_SUG_GL_Z_S_R 0.0595195
#define MAP_SUG_GL_X_A_Y 1.05136
#define MAP_SUG_GL_Y_A_Y -3.42329
#define MAP_SUG_GL_Z_A_Y 0.652454
#define MAP_SUG_GL_X_H_Y 1.05136
#define MAP_SUG_GL_Y_H_Y -3.42329
#define MAP_SUG_GL_Z_H_Y 0.652454
/* #define MAP_SUG_GL_X_H_Y 1.29557 */
/* #define MAP_SUG_GL_Y_H_Y -3.63314  */
/* #define MAP_SUG_GL_Z_H_Y -0.162178 */
#define MAP_SUG_GL_X_S_Y 2.47236
#define MAP_SUG_GL_Y_S_Y -2.90077 
#define MAP_SUG_GL_Z_S_Y 0.154339 


extern double MAP_SUG_GL_X[N_BASES][N_GLYC_STATES];
extern double MAP_SUG_GL_Y[N_BASES][N_GLYC_STATES];
extern double MAP_SUG_GL_Z[N_BASES][N_GLYC_STATES];

double MC_integrate(int, double **, double **, double **);
void MC_update_positions(double **, double **, double **, int);
void MC_update_glypuck(int);
double MC_eval_displacement(int, double **, double **, double **, int, double, double, double);
void MC_perform_nt_rotation(int,int);
void MC_perform_bb_rotation(int);
void MC_perform_nt_translation(int);
void MC_perform_p_translation(int);
int MC_perform_trial_movement(int);
int MC_calculate_local_energy(double *, double *, double *, int, double *, int, int);
void MC_select_rand_nt(int, double *, double *, double *, int*);
void MC_copy_nt(int, double *, double *, double *);
void MC_eval_total_energy(int, double *, double *, double *);
void MC_print_contact_list(int, double *, double *, double *);
void MC_write_contact_list(int, double *, double *, double *, int);
double SA_eval_total_energy(int, double *, double *, double *);
double MC_get_energy(int, double *, double *, double *, int);
void MC_get_sp_bonded_energy(int, int, double*, double* , double*);


void MC_map_sugar_temp(int);
void MC_map_sugar(int , double **, double **, double **);

int is_purine(int);
int is_pyrimidine(int);

int MC_init_flip(int, int, int, int, double *, double *);
int MC_get_temp_glp(int, int);
void MC_assign_temp_glp(int, int);
void MC_perform_base_flip(int, int);
/* double calc_min_dist(double, double, double, double, double, double); */
/* void calc_min_vec(double, double, double, double, double, double, double *, double *); */
int MC_eval_AH_inter(int, double*, int *, int, int, int, int, double, double, double, double, int *);

void MC_print_radius_of_gyration(int, double*, double*, double*);
void MC_print_secondary_structure(int, double*, double*, double *);
int is_canonical(int, int);

/* double boltzmann_dist(double, double); */
/* void fold_coordinate(double *, double, double); */
/* double vec_norm(double, double, double); */
/* float nrran2(); */
/* double max_2(double, double); */
/* double rand_d(double); */
/* int rand_i(int); */
/* double dot_prod(double *, double *); */
/* void svec_prod(double *, double *, double *, double *); */
/* void vec_prod(double *, double *, double *); */
/* void norex_vec_prod(double, double, double, double, double ,double, double *); */


#endif
