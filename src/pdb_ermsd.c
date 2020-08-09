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
int mc_iter;
/****************/




int main(int argc, char **argv) {
  int i,j;
  int base, face, stf,tempinit=0;
  int rand_a;
  double lx, ly, lz;
  /* initialization */
  mc_n=(int)NSOLUTE;
  int nt_n;
  int len=strlen(argv[1]);
  const char *last_three=&argv[1][len-3];
  if(!strcmp(last_three, ".mc")){
    MC_read_params(&mc_iter, &rand_a,-1);
    tempinit=MC_read_checkpoint(&mc_n, &rx, &ry, &rz, &rand_a, -1, argv[1], 0, NULL);
  }
  else if(!strcmp(last_three, "pdb")){
    MC_read_nsolute(&mc_n,-1, argv[1]);
    MC_read_params(&mc_iter, &rand_a,-1);
    MC_initialize_global(mc_n, rand_a,-1);
    MC_initialize_arrays(mc_n,&rx,&ry,&rz);
    MC_read_pdb(mc_n, &rx, &ry, &rz, -1, argv[1]);
  }
  else{
    printf("Unrecognized file extension. It must be pdb or mc.\n");
    exit(ERR_INPUT);
  }


  nt_n=mc_n/N_PARTS_PER_NT;
  MC_init_ermsd_restr(nt_n);
  double ERMSD_SQ, temp;
  get_first_ermsd(&rx, &ry, &rz, nt_n, &ERMSD_SQ, &temp);
  DELTA_ERMSD_SQ=0;
  printf("ERMSD: %lf\n", sqrt(ERMSD_SQ));
  //#endif
  /* INITIALIZE CONFIGURATION FILE */
  
  //MC_initialize_energies(mc_n, rx, ry, rz);
  return 0;
}
