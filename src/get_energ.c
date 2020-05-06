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
#define BAD -10

FILE *md_configs;
char confname[256];

char checkpoint[256];
/* MC variables */
int mc_n;
double *rx,*ry,*rz;
int mc_iter;
/****************/



void print_help();

int main(int argc, char **argv) {
  int i,j,tempinit;
  int base, face, stf;
  double TOT, BB1, BB2, ST_and_BB2, WCE_and_BPH, WCE, NST, BPH;
  char targ[10];
  int rand_a;
  //  double lx, ly, lz;
  //int typ_e=-1;
  
  /* initialization */
  mc_n=(int)NSOLUTE;
  int nt_n;
  if (argc<3){
    print_help();
    exit(ERR_INPUT);
    //strcpy(argv[2],"-t");
    //check if file is binary(.mc) or ascii(.pdb)
    
  }
  
  int len=strlen(argv[1]);
  const char *last_three=&argv[1][len-3];
  if(!strcmp(last_three, ".mc")){
    MC_read_params(&mc_iter, &rand_a,-1);
    tempinit=MC_read_checkpoint(&mc_n, &rx, &ry, &rz, &rand_a, -1, argv[1], 0, NULL);
    MC_initialize_energies(mc_n, rx, ry, rz);
  }
  else if(!strcmp(last_three, "pdb")){
    MC_read_nsolute(&mc_n,-1, argv[1]);
    MC_read_params(&mc_iter, &rand_a,-1);
    //MC_initialize(mc_n, &rx, &ry, &rz, lx, ly, lz, READ_MC_CONF, rand_a, mc_ini);
    MC_initialize_global(mc_n, rand_a,-1);
    MC_initialize_arrays(mc_n,&rx,&ry,&rz);
    MC_read_pdb(mc_n, &rx, &ry, &rz, -1, argv[1]);
    MC_initialize_energies(mc_n, rx, ry, rz);
  }
  else{
    printf("Unrecognized file extension. It must be pdb or mc.\n");
    exit(ERR_INPUT);
  }
  /* before running */
  nt_n=mc_n/N_PARTS_PER_NT;
  if(!strcmp(argv[2], "-t")){
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
  else if(!strcmp(argv[2], "-s"))
    MC_print_secondary_structure(nt_n, rx, ry, rz);
  else if(!strcmp(argv[2], "-r"))
    MC_print_radius_of_gyration(nt_n, rx, ry, rz);
  else if(!strcmp(argv[2], "-a"))
    MC_print_contact_list(nt_n, rx, ry, rz);
  else if(!strcmp(argv[2], "-f"))
    printf("Total E       : %lf\n    ",MC_get_energy(nt_n, rx, ry, rz, 3) );
  else if(!strcmp(argv[2], "-b") ||  !strcmp(argv[2], "-B")){
    BB1=MC_get_energy(nt_n, rx, ry, rz, 0);
    ST_and_BB2=MC_get_energy(nt_n, rx, ry, rz, 1);
    for(base=0;base<N_BASES_SQ;base++){
      for(stf=0;stf<N_STFACES;stf++)
	b_st_well[base][stf]=0;
      nb_st_well[base]=0;
    }
    BB2=MC_get_energy(nt_n, rx, ry, rz, 1);
    if( !strcmp(argv[2], "-b"))
      printf("S-P E         : %lf\n    ",BB1+BB2);
    else if(!strcmp(argv[2], "-B"))
      printf("B stacking E  : %lf\n    ",ST_and_BB2 - BB2);
  }
  else if(!strcmp(argv[2], "-n"))
    printf("NB stacking E : %lf\n    ",MC_get_energy(nt_n, rx, ry, rz, 2));
  else if(!strcmp(argv[2], "-w"))
    printf("WC E          : %lf\n    ",MC_get_energy(nt_n, rx, ry, rz, 3));
  
  else if(!strcmp(argv[2], "") || !strcmp(argv[2], "-h")){
    print_help();
    exit(ERR_INPUT);
  }
  
  else{
    //printf("We require the type of energy to calculate:\n-1  : Annotations\n0  : Total energy\n1  : Backbone energy\n2  : Bonded energy (Includes S-S and Stacking)\n3  : Non-bonded energy\n"); exit(1);
    print_help();
    exit(ERR_INPUT);
}
  if(strcmp(argv[2],"-s")) printf("No clashes found.\n");
  return 0;
}

void print_help(){
  printf("SPQR annotation/energy calculation.\n");
  printf("This requires the file params.spqr to be in the current directory.\n");
  printf("Usage: SPQR_ENERG <input_config> <arg>\n");
  printf("<input_config> can be a .pdb file in the spqr format or a .mc binary configuration.\n");
  printf("Values of <arg> can be specified with \"-c\", where c is a character:\n");
  printf("(s)  : Secondary structure\n");
  printf("(a)  : Full annotations\n");
  printf("(t)  : Total energy (detailed)\n");
  printf("(f)  : Total energy\n");
  printf("(b)  : Backbone energy\n");
  printf("(B)  : Bonded energy (Includes S-S and Stacking)\n");
  printf("(n)  : Non-bonded energy\n");
  printf("(w)  : Base-pairing energy\n");
  printf("(r)  : Radius of gyration\n");

}
