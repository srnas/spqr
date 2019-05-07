#ifndef MCBASEPAIRING_H
#define MCBASEPAIRING_H
#include "mc_global.h"
#include "mc_utils.h"
#include "mc_verlet_lists.h"
#include "mc_energies.h"

#define E_HB 27.5
#define RHB 3.1
#define CHB -0.65
#define WC_FACE_SUGAR 0
#define WC_FACE_WATSCRICK 1
#define WC_FACE_HOOGSTEEN 2
#define WC_FACE_PHOSPHATE 3

//ANGLES ARE IN DEGREES
#define WCFACEANG_00_0u -59.9167
#define WCFACEANG_00_1u 19.525
#define WCFACEANG_00_2u 122.6
#define WCFACEANG_00_3u 159.807
#define WCFACEANG_00_0f -54.3104
#define WCFACEANG_00_1f 39.3
#define WCFACEANG_00_2f 114.3
#define WCFACEANG_00_3f 172.5

#define WCFACEANG_01_0u -53.1893
#define WCFACEANG_01_1u 25.35
#define WCFACEANG_01_2u 110.5
#define WCFACEANG_01_3u 204.261
#define WCFACEANG_01_0f -45.2618
#define WCFACEANG_01_1f 40.6
#define WCFACEANG_01_2f 120.95
#define WCFACEANG_01_3f 175.789

#define WCFACEANG_02_0u -54.2106
#define WCFACEANG_02_1u 28.075
#define WCFACEANG_02_2u 121.8
#define WCFACEANG_02_3u 164.748
#define WCFACEANG_02_0f -35.7055
#define WCFACEANG_02_1f 39.25
#define WCFACEANG_02_2f 119.075
#define WCFACEANG_02_3f 172.5

#define WCFACEANG_03_0u -65.524
#define WCFACEANG_03_1u 11.475
#define WCFACEANG_03_2u 116.45
#define WCFACEANG_03_3u 166.035
#define WCFACEANG_03_0f -61.9277
#define WCFACEANG_03_1f 59.4
#define WCFACEANG_03_2f 113.75
#define WCFACEANG_03_3f 162.518

#define WCFACEANG_10_0u -75.5349
#define WCFACEANG_10_1u 28.25
#define WCFACEANG_10_2u 93.075
#define WCFACEANG_10_3u 148.583
#define WCFACEANG_10_0f -30.7686
#define WCFACEANG_10_1f 20.425
#define WCFACEANG_10_2f 105.05
#define WCFACEANG_10_3f 151.292

#define WCFACEANG_11_0u -25.8449
#define WCFACEANG_11_1u 24.5
#define WCFACEANG_11_2u 120.475
#define WCFACEANG_11_3u 151.851
#define WCFACEANG_11_0f -33.1508
#define WCFACEANG_11_1f 12.35
#define WCFACEANG_11_2f 112.1
#define WCFACEANG_11_3f 187.793

#define WCFACEANG_12_0u -20.8945
#define WCFACEANG_12_1u 18.5
#define WCFACEANG_12_2u 100.975
#define WCFACEANG_12_3u 150.39
#define WCFACEANG_12_0f -24.0556
#define WCFACEANG_12_1f 8.675
#define WCFACEANG_12_2f 122.575
#define WCFACEANG_12_3f 184.395

#define WCFACEANG_13_0u -24.0982
#define WCFACEANG_13_1u 23.575
#define WCFACEANG_13_2u 93.4
#define WCFACEANG_13_3u 270
#define WCFACEANG_13_0f -24.32
#define WCFACEANG_13_1f 37.225
#define WCFACEANG_13_2f 93.925
#define WCFACEANG_13_3f 187.792

#define WCFACEANG_20_0u -75.5349
#define WCFACEANG_20_1u 31.125
#define WCFACEANG_20_2u 95.175
#define WCFACEANG_20_3u 148.583
#define WCFACEANG_20_0f -56.0152
#define WCFACEANG_20_1f 26.45
#define WCFACEANG_20_2f 85.05
#define WCFACEANG_20_3f 169.079

#define WCFACEANG_21_0u -63.9792
#define WCFACEANG_21_1u 23.2
#define WCFACEANG_21_2u 128.125
#define WCFACEANG_21_3u 177.593
#define WCFACEANG_21_0f -60.9851
#define WCFACEANG_21_1f 19.475
#define WCFACEANG_21_2f 135.175
#define WCFACEANG_21_3f 182.816

#define WCFACEANG_22_0u -58.1633
#define WCFACEANG_22_1u 13.525
#define WCFACEANG_22_2u 108
#define WCFACEANG_22_3u 230.927
#define WCFACEANG_22_0f -40.9448
#define WCFACEANG_22_1f 20.725
#define WCFACEANG_22_2f 109.575
#define WCFACEANG_22_3f 209.558

#define WCFACEANG_23_0u -62.7025
#define WCFACEANG_23_1u 13.475
#define WCFACEANG_23_2u 111.725
#define WCFACEANG_23_3u 153.033
#define WCFACEANG_23_0f -59.8717
#define WCFACEANG_23_1f 9.075
#define WCFACEANG_23_2f 93.875
#define WCFACEANG_23_3f 164.348

#define WCFACEANG_30_0u -34.254
#define WCFACEANG_30_1u 22.5
#define WCFACEANG_30_2u 120
#define WCFACEANG_30_3u 158.34
#define WCFACEANG_30_0f -28.4141
#define WCFACEANG_30_1f 48.05
#define WCFACEANG_30_2f 121.475
#define WCFACEANG_30_3f 180.71

#define WCFACEANG_31_0u -3.22669
#define WCFACEANG_31_1u 36.9
#define WCFACEANG_31_2u 115.15
#define WCFACEANG_31_3u 148.572
#define WCFACEANG_31_0f -22.1082
#define WCFACEANG_31_1f 11.35
#define WCFACEANG_31_2f 107.5
#define WCFACEANG_31_3f 174.667

#define WCFACEANG_32_0u -54.2106
#define WCFACEANG_32_1u 18.45
#define WCFACEANG_32_2u 102.15
#define WCFACEANG_32_3u 270
#define WCFACEANG_32_0f -15.6727
#define WCFACEANG_32_1f 11.925
#define WCFACEANG_32_2f 102.575
#define WCFACEANG_32_3f 166.054

#define WCFACEANG_33_0u -34.254
#define WCFACEANG_33_1u 18.6
#define WCFACEANG_33_2u 115.225
#define WCFACEANG_33_3u 169.34
#define WCFACEANG_33_0f -29.7735
#define WCFACEANG_33_1f 29.75
#define WCFACEANG_33_2f 115.95
#define WCFACEANG_33_3f 179.402


#define BPFACEANG_0_0 -65.5
#define BPFACEANG_0_1 30
#define BPFACEANG_0_2 117
#define BPFACEANG_0_3 204.3
#define BPFACEANG_1_0 -54.2
#define BPFACEANG_1_1 19
#define BPFACEANG_1_2 105
#define BPFACEANG_1_3 230
#define BPFACEANG_2_0 -75.5
#define BPFACEANG_2_1 3
#define BPFACEANG_2_2 105
#define BPFACEANG_2_3 151.9
#define BPFACEANG_3_0 -34.3
#define BPFACEANG_3_1 22.2
#define BPFACEANG_3_2 117
#define BPFACEANG_3_3 238.6

#define BPSPECANG_2_0 30
#define BPSPECANG_2_1 52
#define BPSPECANG_3_0 151
#define BPSPECANG_3_1 180

#define BP_OP1X -0.9022973
#define BP_OP1Y -0.9411747
#define BP_OP1Z -0.5219181
#define BP_OP2X -0.10462739
#define BP_OP2Y 1.340639
#define BP_OP2Z -0.4210464


#define BP_CONTACT_G_S_X  2.68918
#define BP_CONTACT_G_S_Y  0.130544
#define BP_CONTACT_G_S_Z  -0.00852924
#define BP_CONTACT_G_W1_X 2.69001
#define BP_CONTACT_G_W1_Y 0.124714
#define BP_CONTACT_G_W1_Z 0.00423577
#define BP_CONTACT_G_W2_X 0.622095
#define BP_CONTACT_G_W2_Y 1.17305
#define BP_CONTACT_G_W2_Z -0.00995796
#define BP_CONTACT_U_W_X  0.680016
#define BP_CONTACT_U_W_Y  1.15549
#define BP_CONTACT_U_W_Z  0.0
#define BP_CONTACT_A_W_X  -1.34088
#define BP_CONTACT_A_W_Y  2.40597
#define BP_CONTACT_A_W_Z  -0.0135445
#define BP_CONTACT_A_H_X  -1.30871
#define BP_CONTACT_A_H_Y  2.41847
#define BP_CONTACT_A_H_Z  0.00488927
#define BP_CONTACT_C_W_X  -1.23645
#define BP_CONTACT_C_W_Y   2.38727
#define BP_CONTACT_C_W_Z  0.000929817
#define BP_CONTACT_C_H_X  -1.24184
#define BP_CONTACT_C_H_Y   2.38439
#define BP_CONTACT_C_H_Z  0.00223364

#define BP_HYDRO_G_S_X  3.29193
#define BP_HYDRO_G_S_Y  -0.695619
#define BP_HYDRO_G_S_Z  0.0148206
#define BP_HYDRO_G_W1_X 3.1827
#define BP_HYDRO_G_W1_Y 1.02883
#define BP_HYDRO_G_W1_Z 0.0325176
#define BP_HYDRO_G_W2_X 1.09028
#define BP_HYDRO_G_W2_Y 2.09202
#define BP_HYDRO_G_W2_Z -0.0585693
#define BP_HYDRO_U_W_X  1.21479
#define BP_HYDRO_U_W_Y  2.04385
#define BP_HYDRO_U_W_Z  0.00232977
#define BP_HYDRO_A_W_X  -0.853091
#define BP_HYDRO_A_W_Y  3.30294
#define BP_HYDRO_A_W_Z  -0.00778549
#define BP_HYDRO_A_H_X  -2.32178
#define BP_HYDRO_A_H_Y  2.54785
#define BP_HYDRO_A_H_Z  -0.0240658
#define BP_HYDRO_C_W_X  -0.707043
#define BP_HYDRO_C_W_Y  3.25897
#define BP_HYDRO_C_W_Z  0.000907987
#define BP_HYDRO_C_H_X  -2.2631
#define BP_HYDRO_C_H_Y  2.47282
#define BP_HYDRO_C_H_Z  0.0474806

extern double bp_contacts[N_BASES][WC_FACES][DIM];
extern double bp_hydros[N_BASES][WC_FACES][DIM];
extern double bp_G_cont_W2[DIM];
extern double bp_G_hydr_W2[DIM];



extern double mc_sameface_max_aWW[N_BASES][N_BASES];
extern double mc_sameface_max_aSS[N_BASES][N_BASES];
extern double mc_sameface_max_pSS[N_BASES][N_BASES];
extern double mc_sameface_min_aWW[N_BASES][N_BASES];
extern double mc_sameface_min_aSS[N_BASES][N_BASES];
extern double mc_sameface_min_pSS[N_BASES][N_BASES];







extern double wc_face_angle[N_BASES_SQ][4];
extern double wc_face_angle_F[N_BASES_SQ][4];

extern int **wc_face_ind, **wc_face_ind_trial;
extern double **wc_face_ene, **wc_face_ene_trial;
extern int **wc_is_paired, **wc_is_paired_trial;

// BASE-PHOSPHATE
extern double bp_face_angle[N_BASES][4];
extern double bp_spec_ang[N_BASES][2];

extern int *bp_face_ind, *bp_face_ind_trial;
extern double *bp_face_ene, *bp_face_ene_trial;
extern int *bp_is_paired, *bp_is_paired_trial;

#ifdef SECONDSTC
extern int *wc_sstruct_N;
extern int **wc_sstruct_neigh;
#endif

int eval_sameface_excl(int, int, int, int, double*, double *, int);

int MC_get_wc_face(int, double *, int);
int MC_get_bp_face(int, double *);
int MC_get_BPh_subface(int, double *, int);

void MC_update_wc_lists(int);
void MC_check_total_wc_list(int);
//void MC_check_local_wc_list(int, double *, double *, double *);
//int MC_is_wc_corresp_trial(int, int);
void MC_is_wc_corresp(int, int, int, int*);
void MC_find_min_wc(int, int, double *, double *, double *, int);
void MC_find_min_wc_per_nt(int, int, int, double*, double*, double*, int, int, int);
double MC_calculate_local_wc_energy(int, int, double *, double *, double *);
double MC_calculate_total_wc_energy(int,int, double *, double *, double *);
void MC_init_wc_arrays(int, double *, double *, double *);
void MC_calc_nnN_watscric(int, int, double *, double *, double , double, int, int, int, int);
double MC_sugar_penalization(int, int, int, int, int, int, int);
double calc_wc_psdihedr(double *, double *, double *, double *, double *, double *, int , int);
double calc_nnN_tab_wc_psdihedr(double , int, int);
double calc_nnN_tab_watscric(double *, int, int, int);
double calc_nnN_tab_watscric_inv(double *, int, int, int);
int MC_update_min_wc(int, int, double *, double *, double *, int);

void MC_calc_npN_phosbase(int, double *, int, int, int, double *, double *, double *, int);
double calc_npN_tab_phosbase(double *, int);
double calc_wc_secdihed(double *, double *, double *, double *, double *, double * , int, int);
int eval_wc_secdihed(double, int, int, int);
void  MC_get_op1op2(double *, double *, double *, double *, double *, double *);
int MC_check_bph_HB(double *, double *, int , int);

#endif
