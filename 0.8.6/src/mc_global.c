#include "mc_global.h"

/* global*/
double mc_target_temp;
double box_l[DIM];
int **mc_pbox;
int *mc_types;
double *mc_temp_x;
double *mc_temp_y;
double *mc_temp_z;
int **mc_temp_pbox;
double MC_NT_XYZ;
double MC_PH_XYZ;
double MC_NT_ANGLE;
double MC_BB_ANGLE;
double MC_NT_ANGLE_COS, MC_NT_ANGLE_SIN;
double MC_BB_ANGLE_COS, MC_BB_ANGLE_SIN;
#ifdef FROZEN
int *fr_is_mobile;
#endif
double bia_lambda;
double bia_cap;
int *mc_glyc;
int *mc_temp_glyc;
double MAP_SUG_GL_X[N_BASES][N_GLYC_STATES];
double MAP_SUG_GL_Y[N_BASES][N_GLYC_STATES];
double MAP_SUG_GL_Z[N_BASES][N_GLYC_STATES];

int *mc_puck;
int *mc_temp_puck;

/* mc_bond_lists */
int **mc_bondlist, **mc_anglelist, **mc_dihedrallist, **mc_nbonds, **mc_tab_bonds_list;
int **mc_anglecenter;
FILE *mc_bond_file;

/* mc_cells */
int *mc_cells;
int mc_n_linked_cells;
int mc_nc[DIM];
double mc_linked_cell_l;

/* mc_checkpoints */
int mc_chk_freq;

/* mc_energies */
int mc_n_types;
int mc_n_bond_types;
int mc_n_ang_types;
int mc_n_dih_types;

double mc_r_cut, mc_r_cut_sq;
double mc_wc_rcut, mc_wc_rcut_sq;
double mc_bph_rcut, mc_bph_rcut_sq;
double mc_nb_rcut, mc_nb_rcut_sq;

double **mc_lj_sig;
double **mc_lj_eps;

double *mc_harm_k;
double *mc_harm_r;


//double *mc_mass;
double *mc_ang_k;
double *mc_ang_th;
double *mc_dih_k, *mc_dih_phi, *mc_dih_n;
FILE *mc_force_file;

/* TABULATED POTENTIAL WELL  - FOR EACH PAIR OF NUCLEOTIDES*/
double nb_st_well[N_BASES_SQ];
double nb_wc_well[N_BASES_SQ][WC_FACES_SQ];
double nb_wc_well_F[N_BASES_SQ][WC_FACES_SQ];
double b_st_well[N_BASES_SQ][N_STFACES];
double BB_PREF;
double BB_PREF_A;
double nb_bp_well[N_BASES][WC_FACES];
double nb_bp_spec_well[N_BASES][3];
double glyc_well[N_GLYC_STATES];
double puck_well[N_PUCK_STATES];
double wc_secdih_min[N_BASES_SQ][WC_FACES_SQ];
double wc_secdih_max[N_BASES_SQ][WC_FACES_SQ];
double wc_secdih_min_F[N_BASES_SQ][WC_FACES_SQ];
double wc_secdih_max_F[N_BASES_SQ][WC_FACES_SQ];
/* TABULATED ENERGIES */
double *table_ssB1_33;
double *table_ssB1_32;
double *table_ssB1_23;
double *table_ssB1_22;
double table_ssB1_params_33[3];
double table_ssB1_params_32[3];
double table_ssB1_params_23[3];
double table_ssB1_params_22[3];
int table_ssB1_N_33;
int table_ssB1_N_32;
int table_ssB1_N_23;
int table_ssB1_N_22;


double **table_nnB_0_s33;
double table_nnB_params_0_s33[N_BASES_SQ][DIM][3];
int table_nnB_N_0_s33[N_BASES_SQ][DIM];
double **table_nnB_0_s35;
double table_nnB_params_0_s35[N_BASES_SQ][DIM][3];
int table_nnB_N_0_s35[N_BASES_SQ][DIM];
double **table_nnB_0_s53;
double table_nnB_params_0_s53[N_BASES_SQ][DIM][3];
int table_nnB_N_0_s53[N_BASES_SQ][DIM];
double **table_nnB_0_s55;
double table_nnB_params_0_s55[N_BASES_SQ][DIM][3];
int table_nnB_N_0_s55[N_BASES_SQ][DIM];


double **table_nnB_1_s33;
double table_nnB_params_1_s33[N_BASES_SQ][3];
int table_nnB_N_1_s33[N_BASES_SQ];
double **table_nnB_1_s35;
double table_nnB_params_1_s35[N_BASES_SQ][3];
int table_nnB_N_1_s35[N_BASES_SQ];
double **table_nnB_1_s53;
double table_nnB_params_1_s53[N_BASES_SQ][3];
int table_nnB_N_1_s53[N_BASES_SQ];
double **table_nnB_1_s55;
double table_nnB_params_1_s55[N_BASES_SQ][3];
int table_nnB_N_1_s55[N_BASES_SQ];


double **table_nnN_0s3;
double table_nnN_params_0s3[N_BASES_SQ][DIM][3];
int table_nnN_N_0s3[N_BASES_SQ][DIM];
double **table_nnN_0s5;
double table_nnN_params_0s5[N_BASES_SQ][DIM][3];
int table_nnN_N_0s5[N_BASES_SQ][DIM];

double **table_nnN_2;
double table_nnN_params_2[N_BASES_SQ][DIM][3];
int table_nnN_N_2[N_BASES_SQ][4];
double **table_nnN_2_inv;
double table_nnN_params_2_inv[N_BASES_SQ][DIM][3];
int table_nnN_N_2_inv[N_BASES_SQ][4];
double **table_nnN_3;
double table_nnN_params_3[N_BASES_SQ][3];
int table_nnN_N_3[N_BASES_SQ];

double **table_nnN_2_F;
double table_nnN_params_2_F[N_BASES_SQ][DIM][3];
int table_nnN_N_2_F[N_BASES_SQ][4];
double **table_nnN_2_inv_F;
double table_nnN_params_2_inv_F[N_BASES_SQ][DIM][3];
int table_nnN_N_2_inv_F[N_BASES_SQ][4];
double **table_nnN_3_F;
double table_nnN_params_3_F[N_BASES_SQ][3];
int table_nnN_N_3_F[N_BASES_SQ];

double **table_bpI_A3;
double **table_bpI_H3;
double **table_bpI_S3;
double table_bpI_params_A3[N_BASES][DIM][3];
double table_bpI_params_H3[N_BASES][DIM][3];
double table_bpI_params_S3[N_BASES][DIM][3];
int table_bpI_N_A3[N_BASES][DIM];
int table_bpI_N_H3[N_BASES][DIM];
int table_bpI_N_S3[N_BASES][DIM];
double **table_bpI_A2;
double **table_bpI_H2;
double **table_bpI_S2;
double table_bpI_params_A2[N_BASES][DIM][3];
double table_bpI_params_H2[N_BASES][DIM][3];
double table_bpI_params_S2[N_BASES][DIM][3];
int table_bpI_N_A2[N_BASES][DIM];
int table_bpI_N_H2[N_BASES][DIM];
int table_bpI_N_S2[N_BASES][DIM];

double **table_bpB_A3;
double **table_bpB_H3;
double **table_bpB_S3;
double table_bpB_params_A3[N_BASES_SQ][DIM][3];
double table_bpB_params_H3[N_BASES_SQ][DIM][3];
double table_bpB_params_S3[N_BASES_SQ][DIM][3];
int table_bpB_N_A3[N_BASES_SQ][DIM];
int table_bpB_N_H3[N_BASES_SQ][DIM];
int table_bpB_N_S3[N_BASES_SQ][DIM];
double **table_bpB_A2;
double **table_bpB_H2;
double **table_bpB_S2;
double table_bpB_params_A2[N_BASES_SQ][DIM][3];
double table_bpB_params_H2[N_BASES_SQ][DIM][3];
double table_bpB_params_S2[N_BASES_SQ][DIM][3];
int table_bpB_N_A2[N_BASES_SQ][DIM];
int table_bpB_N_H2[N_BASES_SQ][DIM];
int table_bpB_N_S2[N_BASES_SQ][DIM];

double **table_npN_0;
double table_npN_params_0[N_BASES][DIM][3];
int table_npN_N_0[N_BASES][DIM];
/* mc_utils */
double *idum;

#ifdef WARMUP
double ph_pintra_ave[N_BASES][3][2][DIM];
double ph_pinter_ave[N_BASES][3][2][DIM];
double ss_ang_ave;
#endif
#ifdef SECONDSTC
int *wc_sstruct_N;
int **wc_sstruct_neigh;
#endif
/* mc_verlet_lists */
int *vl_n_pairs;
int vl_n_pairs_max;
//int vl_flag;
int vl_neighbors_per_atom;
double vl_max_displ;

double vl_skin;
double vl_rmax;
double vl_rmax_sq;
double vl_max_vel;

int **vl_neighbor_tables;
int *vl_count;
int vl_ncrit;

/* mc_base_pairing */
double mc_sameface_max_aWW[N_BASES][N_BASES];
double mc_sameface_max_aSS[N_BASES][N_BASES];
double mc_sameface_max_pSS[N_BASES][N_BASES];
double mc_sameface_min_aWW[N_BASES][N_BASES];
double mc_sameface_min_aSS[N_BASES][N_BASES];
double mc_sameface_min_pSS[N_BASES][N_BASES];
int **wc_face_ind, **wc_face_ind_trial;
double **wc_face_ene, **wc_face_ene_trial;
int **wc_is_paired, **wc_is_paired_trial;
double wc_face_angle[N_BASES_SQ][4];
double wc_face_angle_F[N_BASES_SQ][4];
int *bp_face_ind, *bp_face_ind_trial;
double *bp_face_ene, *bp_face_ene_trial;
int *bp_is_paired, *bp_is_paired_trial;
double bp_face_angle[N_BASES][4];
double bp_spec_ang[N_BASES][2];
