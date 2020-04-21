#ifndef MCFORCES_H
#define MCFORCES_H

#include "mc_global.h"
#include "mc_utils.h"
#include "mc_verlet_lists.h"
//#include "mc_bond_lists.h"
#include "mc_checkpoints.h"
#include "mc_integrate.h"
#include "mc_ermsd.h"

#define MC_CAP 10000.0
#define N_MAX_TYPES 4
#define MAX_BOND_TYPES 6460
#define MAX_ANG_TYPES 1
#define MAX_DIH_TYPES 1
#define DEFAULT_LJ_EPS  1
#define DEFAULT_LJ_SIG  1
#define DEFAULT_HARM_K  5000
#define DEFAULT_HARM_R  1
#define DEFAULT_ANG_K 0
#define DEFAULT_ANG_TH 0
#define DEFAULT_DIH_K 0
#define DEFAULT_DIH_PHI 0
#define DEFAULT_DIH_N 1

//STACKING FACES
#define STFACE3 0
#define STFACE5 1
#define STFACES33 0
#define STFACES35 1
#define STFACES53 2
#define STFACES55 3
//STACKING CRITICAL ANGLE : 23 DEG
#define STANG 0.4

//FLIP FLAGS
#define NT_UNFLIPPED 0
#define NT_FLIPPED 1

//GLYCOSIDIC STATES
#define GLYC_A 0
#define GLYC_H 1
#define GLYC_S 2
#define GLYCS_AA 0
#define GLYCS_AH 1
#define GLYCS_AS 2
#define GLYCS_HA 3
#define GLYCS_HH 4
#define GLYCS_HS 5
#define GLYCS_SA 6
#define GLYCS_SH 7
#define GLYCS_SS 8


//PUCKERS 
#define PUCK_3 0 
#define PUCK_2 1 
 #define PUCKS_33 0 
 #define PUCKS_32 1 
 #define PUCKS_23 2 
 #define PUCKS_22 3 

//EXCLUDED VOLUME
#define EV_GLOB_RCUT 7.0

#define EV_SUGSUG1 0.0625
#define EV_SUGSUG2 0.0625
#define EV_SUGSUG3 0.0625
#define EV_SUGSUG4 0
#define EV_SUGSUG5 0
#define EV_SUGSUG6 0


#define EV_SUGPHO1 0.0516529
#define EV_SUGPHO2 0.0516529
#define EV_SUGPHO3 0.0516529
#define EV_SUGPHO4 0
#define EV_SUGPHO5 0
#define EV_SUGPHO6 0



#define EV_PHOSUG1 0.0516529
#define EV_PHOSUG2 0.0516529
#define EV_PHOSUG3 0.0516529
#define EV_PHOSUG4 0
#define EV_PHOSUG5 0
#define EV_PHOSUG6 0

//#define EV_SUGPHOR 3.5
//#define EV_PHOSUGR 3.5
#define EV_SUGPHOB 0.0816327
#define EV_PHOSUGB 0.0816327

#define EV_PHOS1 0.0416493
#define EV_PHOS2 0.0416493
#define EV_PHOS3 0.0416493
#define EV_PHOS4 0
#define EV_PHOS5 0
#define EV_PHOS6 0

#define EV_PURSUG1 0.063225
#define EV_PURSUG2 0.063225
#define EV_PURSUG3 0.063225
#define EV_PURSUG4 0
#define EV_PURSUG5 0
#define EV_PURSUG6 0

#define EV_PYRSUG1 0.0692521
#define EV_PYRSUG2 0.0692521
#define EV_PYRSUG3 0.0692521
#define EV_PYRSUG4 0
#define EV_PYRSUG5 0
#define EV_PYRSUG6 0


#define EV_SUGPUR1 0.0513906283333 
#define EV_SUGPUR2 0.0379621127778 
#define EV_SUGPUR3 0.0801443361111 
#define EV_SUGPUR4 0.00206081777778 
#define EV_SUGPUR5 -0.00597076944444 
#define EV_SUGPUR6 0.00844257611111
//rescaled parameters 0.05103, 0.0449899, 0.071579, 0.0208221, -0.00299439, 0.000222708
#define EV_SUGPYR1 0.05103
#define EV_SUGPYR2 0.0449899
#define EV_SUGPYR3 0.071579
#define EV_SUGPYR4 0.0208221
#define EV_SUGPYR5 -0.00299439
#define EV_SUGPYR6 0.000222708
/* #define EV_SUGPYR1 0.047731901724 */
/* #define EV_SUGPYR2 0.040628712808 */
/* #define EV_SUGPYR3 0.071561573448 */
/* #define EV_SUGPYR4 0.0246147088 */
/* #define EV_SUGPYR5 -0.003234166 */
/* #define EV_SUGPYR6 0.0004984296 */


#define EV_BASE1 0.04
#define EV_BASE2 0.04
#define EV_BASE3 0.1111111

#define EV_PURPHO1 0.0612685
#define EV_PURPHO2 0.0612685
#define EV_PURPHO3 0.0612685
#define EV_PURPHO4 0
#define EV_PURPHO5 0
#define EV_PURPHO6 0

#define EV_PYRPHO1 0.0657462
#define EV_PYRPHO2 0.0657462
#define EV_PYRPHO3 0.0657462
#define EV_PYRPHO4 0
#define EV_PYRPHO5 0
#define EV_PYRPHO6 0


#define EV_PHOPUR1 0.045221901 
#define EV_PHOPUR2 0.040814517 
#define EV_PHOPUR3 0.0780271643333 
#define EV_PHOPUR4 0.0135575276667 
#define EV_PHOPUR5 0.00325463733333 
#define EV_PHOPUR6 -0.00704502666667

//bonded parameters
/* #define EV_PHOPYR1 0.0502646633333 */
/* #define EV_PHOPYR2 0.0387524473333 */
/* #define EV_PHOPYR3 0.0720512393333 */
/* #define EV_PHOPYR4 0.01204247 */
/* #define EV_PHOPYR5 0.009543366 */
/* #define EV_PHOPYR6 0.01187689 */

/* #define EV_PHOPYR1 0.0511947311703  */
/* #define EV_PHOPYR2 0.0403455973541 */
/* #define EV_PHOPYR3 0.0745895213784 */
/* #define EV_PHOPYR4 0.0182336901685 */
/* #define EV_PHOPYR5 -0.0125180762 */
/* #define EV_PHOPYR6 0.0232693618919 */
//rescaled parameters
#define EV_PHOPYR1 0.0544948
#define EV_PHOPYR2 0.0472098
#define EV_PHOPYR3 0.0765752
#define EV_PHOPYR4 0.0134742 
#define EV_PHOPYR5 -0.00995821 
#define EV_PHOPYR6 0.0195775


// VERLET LISTS
#define DEFAULT_MC_RCUT 21
#define DEFAULT_MC_WC_RCUT 11
#define DEFAULT_MC_BPH_RCUT 8
#define DEFAULT_MC_NB_RCUT 21

#define DEFAULT_MC_N_TYPES 1
#define DEFAULT_MC_N_BOND_TYPES 1
#define DEFAULT_VL_SKIN 2


//NUMERICS
#define TINY_SIN_VAL 1E-5



extern int mc_n_types;
extern int mc_n_bond_types;
extern int mc_n_ang_types;
extern int mc_n_dih_types;

extern double mc_r_cut, mc_r_cut_sq;
extern double mc_wc_rcut, mc_wc_rcut_sq;
extern double mc_bph_rcut, mc_bph_rcut_sq;
extern double mc_nb_rcut, mc_nb_rcut_sq;

extern double **mc_lj_sig;
extern double **mc_lj_eps;

extern double *mc_harm_k;
extern double *mc_harm_r;

extern int *mc_glyc;
extern int *mc_temp_glyc;

extern int *mc_puck;
extern int *mc_temp_puck;

extern int *glp_is_flippable;

//double *mc_mass;
extern double *mc_ang_k;
extern double *mc_ang_th;
extern double *mc_dih_k, *mc_dih_phi, *mc_dih_n;
extern FILE *mc_force_file;

/* TABULATED POTENTIAL WELL  - FOR EACH PAIR OF NUCLEOTIDES*/
extern double nb_st_well[N_BASES_SQ];
extern double nb_wc_well[N_BASES_SQ][WC_FACES_SQ];
extern double nb_wc_well_F[N_BASES_SQ][WC_FACES_SQ];
extern double b_st_well[N_BASES_SQ][N_STFACES];
extern double nb_bp_well[N_BASES][WC_FACES];
extern double nb_bp_spec_well[N_BASES][3];
extern double BB_PREF;
extern double BB_PREF_A;
extern double glp_well_R[N_GLYC_STATES][N_PUCK_STATES];
//extern double puck_well_R[N_PUCK_STATES];
extern double glp_well_Y[N_GLYC_STATES][N_PUCK_STATES];
//extern double puck_well_Y[N_PUCK_STATES];
extern double wc_secdih_min[N_BASES_SQ][WC_FACES_SQ];
extern double wc_secdih_max[N_BASES_SQ][WC_FACES_SQ];
extern double wc_secdih_min_F[N_BASES_SQ][WC_FACES_SQ];
extern double wc_secdih_max_F[N_BASES_SQ][WC_FACES_SQ];
/* TABULATED ENERGIES */
extern double *table_ssB1_33;
extern double *table_ssB1_32;
extern double *table_ssB1_23;
extern double *table_ssB1_22;
extern double table_ssB1_params_33[3];
extern double table_ssB1_params_32[3];
extern double table_ssB1_params_23[3];
extern double table_ssB1_params_22[3];
extern int table_ssB1_N_33;
extern int table_ssB1_N_32;
extern int table_ssB1_N_23;
extern int table_ssB1_N_22;


extern double **table_nnB_0_s33;
extern double table_nnB_params_0_s33[N_BASES_SQ][DIM][3];
extern int table_nnB_N_0_s33[N_BASES_SQ][DIM];
extern double **table_nnB_0_s35;
extern double table_nnB_params_0_s35[N_BASES_SQ][DIM][3];
extern int table_nnB_N_0_s35[N_BASES_SQ][DIM];
extern double **table_nnB_0_s53;
extern double table_nnB_params_0_s53[N_BASES_SQ][DIM][3];
extern int table_nnB_N_0_s53[N_BASES_SQ][DIM];
extern double **table_nnB_0_s55;
extern double table_nnB_params_0_s55[N_BASES_SQ][DIM][3];
extern int table_nnB_N_0_s55[N_BASES_SQ][DIM];


extern double **table_nnB_1_s33;
extern double table_nnB_params_1_s33[N_BASES_SQ][3];
extern int table_nnB_N_1_s33[N_BASES_SQ];
extern double **table_nnB_1_s35;
extern double table_nnB_params_1_s35[N_BASES_SQ][3];
extern int table_nnB_N_1_s35[N_BASES_SQ];
extern double **table_nnB_1_s53;
extern double table_nnB_params_1_s53[N_BASES_SQ][3];
extern int table_nnB_N_1_s53[N_BASES_SQ];
extern double **table_nnB_1_s55;
extern double table_nnB_params_1_s55[N_BASES_SQ][3];
extern int table_nnB_N_1_s55[N_BASES_SQ];


extern double **table_nnN_0s3;
extern double table_nnN_params_0s3[N_BASES_SQ][DIM][3];
extern int table_nnN_N_0s3[N_BASES_SQ][DIM];
extern double **table_nnN_0s5;
extern double table_nnN_params_0s5[N_BASES_SQ][DIM][3];
extern int table_nnN_N_0s5[N_BASES_SQ][DIM];

extern double **table_nnN_2;
extern double table_nnN_params_2[N_BASES_SQ][DIM][3];
extern int table_nnN_N_2[N_BASES_SQ][4];
extern double **table_nnN_2_inv;
extern double table_nnN_params_2_inv[N_BASES_SQ][DIM][3];
extern int table_nnN_N_2_inv[N_BASES_SQ][4];
extern double **table_nnN_3;
extern double table_nnN_params_3[N_BASES_SQ][3];
extern int table_nnN_N_3[N_BASES_SQ];

extern double **table_nnN_2_F;
extern double table_nnN_params_2_F[N_BASES_SQ][DIM][3];
extern int table_nnN_N_2_F[N_BASES_SQ][4];
extern double **table_nnN_2_inv_F;
extern double table_nnN_params_2_inv_F[N_BASES_SQ][DIM][3];
extern int table_nnN_N_2_inv_F[N_BASES_SQ][4];
extern double **table_nnN_3_F;
extern double table_nnN_params_3_F[N_BASES_SQ][3];
extern int table_nnN_N_3_F[N_BASES_SQ];

extern double **table_bpI_A3;
extern double **table_bpI_H3;
extern double **table_bpI_S3;
extern double table_bpI_params_A3[N_BASES][DIM][3];
extern double table_bpI_params_H3[N_BASES][DIM][3];
extern double table_bpI_params_S3[N_BASES][DIM][3];
extern int table_bpI_N_A3[N_BASES][DIM];
extern int table_bpI_N_H3[N_BASES][DIM];
extern int table_bpI_N_S3[N_BASES][DIM];
extern double **table_bpI_A2;
extern double **table_bpI_H2;
extern double **table_bpI_S2;
extern double table_bpI_params_A2[N_BASES][DIM][3];
extern double table_bpI_params_H2[N_BASES][DIM][3];
extern double table_bpI_params_S2[N_BASES][DIM][3];
extern int table_bpI_N_A2[N_BASES][DIM];
extern int table_bpI_N_H2[N_BASES][DIM];
extern int table_bpI_N_S2[N_BASES][DIM];

extern double **table_bpB_A3;
extern double **table_bpB_H3;
extern double **table_bpB_S3;
extern double table_bpB_params_A3[N_BASES_SQ][DIM][3];
extern double table_bpB_params_H3[N_BASES_SQ][DIM][3];
extern double table_bpB_params_S3[N_BASES_SQ][DIM][3];
extern int table_bpB_N_A3[N_BASES_SQ][DIM];
extern int table_bpB_N_H3[N_BASES_SQ][DIM];
extern int table_bpB_N_S3[N_BASES_SQ][DIM];
extern double **table_bpB_A2;
extern double **table_bpB_H2;
extern double **table_bpB_S2;
extern double table_bpB_params_A2[N_BASES_SQ][DIM][3];
extern double table_bpB_params_H2[N_BASES_SQ][DIM][3];
extern double table_bpB_params_S2[N_BASES_SQ][DIM][3];
extern int table_bpB_N_A2[N_BASES_SQ][DIM];
extern int table_bpB_N_H2[N_BASES_SQ][DIM];
extern int table_bpB_N_S2[N_BASES_SQ][DIM];


extern double **table_npN_0;
extern double table_npN_params_0[N_BASES][DIM][3];
extern int table_npN_N_0[N_BASES][DIM];

#ifdef WARMUP
extern double ph_pintra_ave[N_BASES][3][2][DIM];
extern double ph_pinter_ave[N_BASES][3][2][DIM];
extern double ss_ang_ave;
#endif

void calc_rel_pos(double *, double *);
void calc_rel_pos_inv(double *, int, double *, double *, double *, double *);
double calc_nnB_tab_stacking(double *, int, int);
double calc_nnB_tab_stacking_inv(double *, int, int);
double calc_st_psdihedr(double *, double *, double *, double*, double*, double*, int, int);
double calc_nnB_tab_st_psdihedr(double, int, int);

double calc_nnN_tab_stacking_s3(double *, int);
double calc_nnN_tab_stacking_s5(double *, int);

double calc_ssB1_sugars(double, int);

void MC_initialize_tabulated_energies(int);
double MC_calc_nnB_stacking(int, int, double*, double*, double, int *);
double MC_calc_nnN_stacking(int, int, double *, double *, int *);
//double MC_calc_nb_stacking_no_dih(int, int, double*, double*, int *);
//double MC_calc_nb_watscric_no_dih(int, double*, double*, int *);


/* ENERGIES */
void MC_initialize_energy_parameters(int);
void MC_read_energy_parameters();
void MC_read_bin_energy_tables();
void MC_read_write_energy_tables();
double MC_calc_bonded_energy(int, double *, double *, double *, int*, int, int,double , double,
int, int,double , double,
			     int *, int*, double *, double *);
double MC_calc_non_bonded_energy(int, double *, double *, double *, int, double *, double, int *);

//double energy_hardcore(double *, double, double, double);

double MC_calc_BP_intra(int, double *, int, int, int*);
double MC_calc_BP_inter(int, double *, int, int, int*);
double MC_calc_intra_energy(int, int*, int *, double *, double *);

double energy_excludedvol(double, double);
double energy_harmonic(double, double, double);
double energy_angle(double, double, double);
//double dihedral_prefactor(double, double, double, double);

void MC_free_energy_params();

/* GLYCOSIDICS AND PUCKERS */
int MC_get_glycs(int, int);
int MC_get_pucks(int, int);
void MC_init_glycs_and_pucks(int, int);


static inline double energy_hardcore(double *vec, double sxx, double syy, double szz, double sxy, double sxz, double syz){
  double RET=MC_CAP+1;
#ifdef WARMUP
  //RET = (0.5*MC_CAP*exp(-sqrt(vec[0]*vec[0]*sxx+vec[1]*vec[1]*syy+vec[2]*vec[2]*szz)));
  //RET=MC_CAP/(vec[0]*vec[0]*sxx + vec[1]*vec[1]*syy + vec[2]*vec[2]*szz + 2*vec[0]*vec[1]*sxy + 2*vec[0]*vec[2]*sxz + 2*vec[1]*vec[2]*syz);
  RET=-MC_CAP*((vec[0]*vec[0]*sxx + vec[1]*vec[1]*syy + vec[2]*vec[2]*szz + 2*vec[0]*vec[1]*sxy + 2*vec[0]*vec[2]*sxz + 2*vec[1]*vec[2]*syz)-1)+MC_CAP;
#endif
  if( vec[0]*vec[0]*sxx + vec[1]*vec[1]*syy + vec[2]*vec[2]*szz + 2*vec[0]*vec[1]*sxy + 2*vec[0]*vec[2]*sxz + 2*vec[1]*vec[2]*syz   > 1.0 )
    RET= 0;
  //else printf("%lf %lf\n", vec[0]*vec[0]*sxx + vec[1]*vec[1]*syy + vec[2]*vec[2]*szz + 2*vec[0]*vec[1]*sxy + 2*vec[0]*vec[2]*sxz + 2*vec[1]*vec[2]*syz ,
  //	      sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]) );
  return RET;
}
#endif
