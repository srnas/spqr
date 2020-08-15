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
#define N_PUCK_STATES 2
#define DIM 3
#define DIMSQ 9
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
#define ERR_WRITING 9

#define MAX_BUFFER 1000
//bonds
#define N_BONDED_INTERACTIONS 2

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

#define GLP_FIXED  0
#define GLP_GLYC   1
#define GLP_PUCK   2
#define GLP_BOTH   3

#define FLIP_GLYC 0
#define FLIP_PUCK 1

#define ANN_NPARAMS 8

/* RANDOM NUMBERS */
/***From Numerical Recipes ***/
#define mc_ecapping -1
#define IM1 2147483563 
#define IM2 2147483399 
#define AM (1.0/IM1) 
#define IMM1 (IM1-1) 
#define IA1 40014 
#define IA2 40692 
#define IQ1 53668 
#define IQ2 52774 
#define IR1 12211 
#define IR2 3791 
#define NTAB 32 
#define NDIV (1+IMM1/NTAB) 
#define EPS 1.2e-7 
#define RNMX (1.0-EPS) 
#define SQ(x) ((x)*(x))
#define MAXSTR 10000


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

extern int **mc_bondlist;
extern int **mc_anglelist;
extern int **mc_dihedrallist;
extern int **mc_nbonds;
extern int **mc_tab_bonds_list;
extern int **mc_anglecenter;
extern FILE *mc_bond_file;

extern char ENERG_PATH[MAX_BUFFER];

#ifdef FROZEN
extern int *fr_is_mobile;
#endif
extern double bia_lambda;
extern double bia_cap;

#endif


//gyration radius
extern double KRG;
extern double RG_target;

//wall
extern double *wall_epsilon;
extern double *wall_sigma;
extern double *wall_A;
extern double *wall_B;
extern double *wall_C;
extern double *wall_D;
extern double *wall_MODSQ;
extern double WALL_ENERG;
extern double DELTA_WALL_ENERG;
extern int *WALL_TYPE;
extern int N_WALLS;

extern double UCM[DIM];
extern double UCMK0,UCMK1,UCMK2;
extern int UMBRELLA_TYPE, N_UMBRELLAS;
