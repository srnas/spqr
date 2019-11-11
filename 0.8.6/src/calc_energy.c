#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mc.h"
//#include "bia_impl.h"
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
  int nt_n;
  if (argc<4){
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
  //printf("nargs= %d \n", argc);
  if(argc==3){
    cmapflag=atoi(argv[2]);
  }
  /* else if (argc==4){  */
  /*   cmapflag=1; */
  /*   tl=atof(argv[2]); */
  /*   tc=atof(argv[3]); */
  /*   printf("lambda = %lf, cap = %lf\n", tl, tc); */
  /*   BIA_set_params(tl, tc); */
  /* } */
  /* else{ */
  /*   cmapflag=1; */
  /*   tl=1; */
  /*   tc=-1000000; */
  /*   BIA_set_params(tl, tc); */
  /* } */
  
  /* before running */
  printf("Analyze %d MC steps every %d trials.\nStarting from iteration %d\n", mc_iter, mc_chk_freq, mc_ini);
  /******************/
  double energ_bb1, energ_bb2, energ_nb, energ_to, energ_wc, energ_bb2_st;
  int base, stf;
  double BACK_BST[N_BASES_SQ][N_STFACES], BACK_NST[N_BASES_SQ];
  for(base=0;base<N_BASES_SQ;base++){
    for(stf=0;stf<N_STFACES;stf++)
      BACK_BST[base][stf]=b_st_well[base][stf];
    BACK_NST[base]=nb_st_well[base];
  }
  
  
  nt_n=mc_n/N_PARTS_PER_NT;
  for(i=mc_ini;i<mc_iter/cmapflag;++i) {
    if(argc==3 && i==0){
      MC_read_ith_configuration();
      printf("\n");
      MC_print_contact_list(nt_n, rx, ry, rz);
    }
    for(j=0;j<cmapflag;j++){
      MC_read_ith_configuration();
      //if(i==0 && j==0 && argc==2){
      for(base=0;base<N_BASES_SQ;base++){
	for(stf=0;stf<N_STFACES;stf++)
	  b_st_well[base][stf]=BACK_BST[base][stf];
	nb_st_well[base]=BACK_NST[base];
      }
      
      printf("%d\t", i*cmapflag+j+1);
      energ_bb1=MC_get_energy(nt_n, rx, ry, rz, 0);
      energ_bb2_st=MC_get_energy(nt_n, rx, ry, rz, 1);
      energ_nb=MC_get_energy(nt_n, rx, ry, rz, 2);
      energ_to=MC_get_energy(nt_n, rx, ry, rz, 3);
      energ_wc=MC_get_energy(nt_n, rx, ry, rz, 4);
      
      for(base=0;base<N_BASES_SQ;base++){
	for(stf=0;stf<N_STFACES;stf++)
	  b_st_well[base][stf]=0;
	nb_st_well[base]=0;
      }
      energ_bb2=MC_get_energy(nt_n, rx, ry, rz, 1);
      
      printf("%lf %lf %lf %lf   %lf\n", energ_bb1+energ_bb2, energ_bb2_st-energ_bb2, energ_nb, energ_wc, energ_to);
      //
      for(base=0;base<N_BASES_SQ;base++){
	for(stf=0;stf<N_STFACES;stf++)
	  b_st_well[base][stf]=BACK_BST[base][stf];
	nb_st_well[base]=BACK_NST[base];
      }
      
    }
    if(argc==3) printf("# ");
        if(argc==3){
      printf("\n");
      MC_print_contact_list(nt_n, rx, ry, rz);
    }
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
