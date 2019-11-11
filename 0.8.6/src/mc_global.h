/* THIS IS FROM mc_0.1 */
#ifndef MCGLOBAL_H
#define MCGLOBAL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
//#define MPIMC
#define EXCLUSIONS
//#define MCCELLS
//#define STWCBONDED
//#define MCVLISTS
#define MCDILUTE

//#define MCCAPPING
//#define DIHEDRALS
#define WC_FACES 3
#define WC_FACES_SQ 9
#define N_PARTS_PER_NT 5
#define N_BASES 4
#define N_BASES_SQ 16
#define N_GLYC_STATES 3
#define N_PUCK_STATES 3
#define DIM 3
#define IPHO 4
#define ISUG 3
#define IBAS 0
#define IX 1
#define IY 2
#define TYP_ADENINE  0
#define TYP_URACIL   1
#define TYP_GUANINE  2
#define TYP_CYTOSINE 3
#define N_STFACES 4

#define ERR_INIT 2
#define ERR_INTEG 3
#define ERR_ANALYSIS 4
#define ERR_INPUT 5
#define ERR_MC 6
#define ERR_NEIGHBOR 7
#define ERR_FORCES 8

/* DEFAULTS */
#define MC_NT_XYZ_DEF 0.3
#define MC_PH_XYZ_DEF 0.6
#define MC_NT_ANGLE_DEF 0.1
#define MC_BB_ANGLE_DEF 0.1
#define MC_TEMP_DEF 1.0
#define MC_RAND_SEED_DEF 1.0
#define MC_ITER_DEF 1000

#define FR_MOB_FULL 3
#define FR_MOB_PHOS 1
#define FR_MOB_BASE 2
#define FR_MOB_FROZ 0

/* MACROS */
#if 1==DIM
#define index(ic,nc) ((ic)[0])
#elif 2==DIM
#define index(ic,nc) ((ic)[0] + (nc)[0]*(ic)[1])
#elif 3==DIM
#define index(ic,nc) ((ic)[0] + (nc)[0]*((ic)[1]+ (nc)[1]*(ic)[2]))
#endif


extern double mc_target_temp;
extern double box_l[DIM];
extern int **mc_pbox;
extern int *mc_types;
extern double *mc_temp_x;
extern double *mc_temp_y;
extern double *mc_temp_z;
extern int **mc_temp_pbox;
extern double MC_NT_XYZ;
extern double MC_PH_XYZ;
extern double MC_NT_ANGLE;
extern double MC_BB_ANGLE;
extern double MC_NT_ANGLE_COS, MC_NT_ANGLE_SIN;
extern double MC_BB_ANGLE_COS, MC_BB_ANGLE_SIN;



#ifdef FROZEN
extern int *fr_is_mobile;
#endif
extern double bia_lambda;
extern double bia_cap;

#endif
