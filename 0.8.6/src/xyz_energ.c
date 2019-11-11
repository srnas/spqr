#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mc.h"

#define LX 30
#define LY 30
#define LZ 30
#define NSOLUTE 868
#define MC_ITER 20000
#define EVERYMILLION 500000

FILE *md_configs;
//char confname[21]="configs/confs.p00.mc";
char confname[256];

char checkpoint[256];
/* MC variables */
int mc_n;
double *rx,*ry,*rz;
int mc_iter, mc_ini;
/****************/




int main(int argc, char **argv) {
  int i,j;
  int base, face, stf;
  double TOT, BB1, BB2, ST_and_BB2, WCE_and_BPH, WCE, NST, BPH;
  int rand_a;
  double lx, ly, lz;
  int typ_e=atoi(argv[1]);
  /* initialization */
  mc_n=(int)NSOLUTE;
  int nt_n;
  if (argc==2){
    MC_read_nsolute(&mc_n,20);
    mc_ini=0;
    MC_read_params(&lx, &ly, &lz, &mc_iter, &rand_a,20);
    //MC_initialize(mc_n, &rx, &ry, &rz, lx, ly, lz, READ_MC_CONF, rand_a, mc_ini);
    MC_initialize_global(mc_n, lx, ly, lz, rand_a,20);
    //if(mc_read_conf_flag==READ_MC_CONF){
    MC_initialize_positions(mc_n,&rx,&ry,&rz,READ_MC_CONF,20);
    /* INITIALIZE CONFIGURATION FILE */
    //MC_initialize_save_configs(mc_n, ini);
    MC_initialize_energies(mc_n, rx, ry, rz);

    
    
  }  else {
    printf("We require the type of energy to calculate:\n-1  : Annotations\n0  : Total energy (detailed)\n1  : Total energy\n2  : Backbone energy\n3  : Bonded energy (Includes S-S and Stacking)\n4  : Non-bonded energy\n");
    exit(1);
  }
  /* before running */
  nt_n=mc_n/N_PARTS_PER_NT;
#ifdef NEW_BIA
  mc_target_temp=1;
#endif
  if(typ_e<0)
    MC_print_contact_list(nt_n, rx, ry, rz);
  else if(typ_e==0){
    
    TOT=MC_get_energy(nt_n, rx, ry, rz, 3);
    BB1=MC_get_energy(nt_n, rx, ry, rz, 0);
      
    
    WCE_and_BPH=MC_get_energy(nt_n, rx, ry, rz, 4);
    for(base=0;base<N_BASES;base++){
      for(face=0;face<WC_FACES;face++){
	nb_bp_well[base][face]=0;
      }
      nb_bp_spec_well[base][0]=0;
      nb_bp_spec_well[base][1]=0;
    }
    WCE=MC_get_energy(nt_n, rx, ry, rz, 4);
    BPH=WCE_and_BPH-WCE;
    NST=MC_get_energy(nt_n, rx, ry, rz, 2);
    ST_and_BB2=MC_get_energy(nt_n, rx, ry, rz, 1);
    
    for(base=0;base<N_BASES_SQ;base++){
      for(stf=0;stf<N_STFACES;stf++)
	b_st_well[base][stf]=0;
      nb_st_well[base]=0;
    }
    BB2=MC_get_energy(nt_n, rx, ry, rz, 1);
    
    printf("TOTAL ENERGY \t\t: %lf\nBackbone E    \t\t: %lf\nStacking(B) E     \t: %lf\nStacking(NB) E \t\t: %lf\nBase-pairing(WC) E\t: %lf\nBase-Phosphate(BPH) E  \t: %lf\n",
	   TOT,
	   BB1+BB2,
	   ST_and_BB2-BB2,
	   NST,
	   WCE,
	   BPH);
  }
  else if(typ_e==1)
    printf("Total E       : %lf\n    ",MC_get_energy(nt_n, rx, ry, rz, 3) );
  else if(typ_e==2 || typ_e==3){
    BB1=MC_get_energy(nt_n, rx, ry, rz, 0);
    ST_and_BB2=MC_get_energy(nt_n, rx, ry, rz, 1);
    for(base=0;base<N_BASES_SQ;base++){
      for(stf=0;stf<N_STFACES;stf++)
	b_st_well[base][stf]=0;
      nb_st_well[base]=0;
    }
    BB2=MC_get_energy(nt_n, rx, ry, rz, 1);
    if(typ_e==2)
      printf("S-P E         : %lf\n    ",BB1+BB2);
    else if(typ_e==3)
      printf("B stacking E  : %lf\n    ",ST_and_BB2 - BB2);
  }
  else if(typ_e==4)
    printf("NB stacking E : %lf\n    ",MC_get_energy(nt_n, rx, ry, rz, 2));
  else if(typ_e==5)
    printf("WC E          : %lf\n    ",MC_get_energy(nt_n, rx, ry, rz, 3));
  
  
  else{printf("We require the type of energy to calculate:\n-1  : Annotations\n0  : Total energy\n1  : Backbone energy\n2  : Bonded energy (Includes S-S and Stacking)\n3  : Non-bonded energy\n"); exit(1);}
  return 0;
}
