#ifndef MCBONDEDF
#define MCBONDEDF

#include "mc_global.h"

#ifdef DIHEDRALS
#define N_BONDED_INTERACTIONS 3
#elif defined STWCBONDED
#define N_BONDED_INTERACTIONS 4
#else 
#define N_BONDED_INTERACTIONS 2
#endif



extern int **mc_bondlist;
extern int **mc_anglelist;
extern int **mc_dihedrallist;
extern int **mc_nbonds;
extern int **mc_tab_bonds_list;
extern int **mc_anglecenter;
extern FILE *mc_bond_file;

int MC_open_bondlist_file();
void MC_read_bondlist_from_file(int, int);
void MC_default_bonds(int);
void MC_close_bondlist_file();
#endif
