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

int cmapflag;

void MC_open_analysis_files(int);
void MC_read_ith_configuration();
void MC_close_analysis_files();

int main(int argc, char **argv) {
  int i,j;
  int rand_a;
  double lx, ly, lz;
  int init=atoi(argv[1]);
  /* initialization */
  mc_n=(int)NSOLUTE;
  int nt_n, nt_c;
  if (argc<6){
    MC_read_nsolute(&mc_n,0);
    mc_ini=0;
    MC_read_params(&lx, &ly, &lz, &mc_iter, &rand_a,0);
    //MC_initialize(mc_n, &rx, &ry, &rz, lx, ly, lz, READ_MC_CONF, rand_a, mc_ini);
    MC_initialize_global(mc_n, lx, ly, lz, rand_a,0);
    //if(mc_read_conf_flag==READ_MC_CONF){
    MC_initialize_positions(mc_n,&rx,&ry,&rz,READ_MC_CONF,0);
    /* INITIALIZE CONFIGURATION FILE */
    //MC_initialize_save_configs(mc_n, ini);
    MC_initialize_energies(mc_n, rx, ry, rz);
    MC_open_analysis_files(init);
    
    
  }  else {
    printf("Default case not implemented.\n");
    exit(1);
  }
  double tl, tc;
  printf("nargs= %d \n", argc);
    
  /* before running */
  printf("Analyze %d MC steps every %d trials.\nStarting from iteration %d\n", mc_iter, mc_chk_freq, mc_ini);
  /******************/
  int base;
    
  nt_n=mc_n/N_PARTS_PER_NT;
  for(i=mc_ini;i<mc_iter/cmapflag;++i) {
    MC_read_ith_configuration();
    for(nt_c=0;nt_c<nt_n;nt_c++)
      printf("%d: %d   %lf %lf %lf\n", i, nt_c, rx[N_PARTS_PER_NT*nt_c], rx[N_PARTS_PER_NT*nt_c], rx[N_PARTS_PER_NT*nt_c]);
  }
  
  MC_close_analysis_files();
  return 0;
}

void MC_open_analysis_files(int it){
  int i;
  int temp_n;
  sprintf(confname , "configs/confs.p%02d.mc", it);
  if((md_configs=fopen(confname, "rb"))==NULL){
    fprintf(stderr,"File %s not found!\n", confname);
    exit(1);
  }
  
  fread(&temp_n, sizeof(int),1,md_configs);
  int *temp_type;
  temp_type=(int*)malloc(sizeof(int)*mc_n);
  if(temp_n != mc_n){
    printf("Number of particles does not match the previous value!\n");
    exit(1);
  }
  for(i=0;i<temp_n;i++){
    fread(&temp_type[i], sizeof(int),1,md_configs);
  }
}

void MC_read_ith_configuration(){
  /* printf("Reading comfiguration...\n"); */
  /* printf("box_l = %lf   %lf   %lf\n", box_l[0], box_l[1], box_l[2]); */
  /* fflush(stdout); */
  int i,d;
  double tempp[DIM];
  double dtempp[DIM];
  for(i=0;i<mc_n;i++){
    fread(tempp,sizeof(double),DIM,md_configs);
    //printf("%f %f %f\n", tempp[0], tempp[1], tempp[2]);fflush(stdout);
    for(d=0;d<DIM;d++){
      //rx[i]=(double)tempp[d];
      dtempp[d]=(double)tempp[d];
      mc_pbox[i][d]=0;
#ifdef PBC
      mc_pbox[i][d]=(int)floor(dtempp[d]/box_l[d]);
      fold_coordinate(&(dtempp[d]),box_l[d],0);
#endif
    }
    rx[i]=dtempp[0];ry[i]=dtempp[1];rz[i]=dtempp[2];
  }
}

void MC_close_analysis_files(){
  fclose(md_configs);
}
