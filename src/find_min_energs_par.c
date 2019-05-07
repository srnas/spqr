#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mc.h"

#ifdef MPIMC
#include <mpi.h>
#endif

#define LX 30
#define LY 30
#define LZ 30
#define NSOLUTE 868
#define MC_ITER 20000
#define EVERYMILLION 500000

FILE *md_configs;
char confname[256];

char checkpoint[256];
/* MC variables */
int mc_n;
double *rx,*ry,*rz;
int mc_iter, mc_ini;
/****************/

int cmapflag=1;
//double temperature;
//double energy_t;
void MC_open_analysis_files(int);
void MC_read_ith_configuration(double *, double *);
void MC_close_analysis_files();

int main(int argc, char **argv) {
  int i,j;
  int rand_a;
  double lx, ly, lz;
  //int init=atoi(argv[1]);
  double minenerg=0, energ;
  int mc_block;
  double energy_t;
  double temperature;
  /* initialization */
    /*MPI STUFF*/
  int mpi_id=0;
#ifdef MPIMC
  int mpi_count;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_count);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
  double buff[2], buff1[2],buff2[2];
  int init=mpi_id;
#endif
  /***********/
  
  mc_n=(int)NSOLUTE;
  int nt_n;
  if (argc==3){
    MC_read_nsolute(&mc_n,0, NULL);
    mc_ini=atoi(argv[1]);
    mc_block=atoi(argv[2]);
    MC_read_params(&lx, &ly, &lz, &mc_iter, &rand_a,0);
    //MC_initialize(mc_n, &rx, &ry, &rz, lx, ly, lz, READ_MC_CONF, rand_a, mc_ini);
    MC_initialize_global(mc_n, lx, ly, lz, rand_a,0);
    //if(mc_read_conf_flag==READ_MC_CONF){
    MC_initialize_arrays(mc_n,&rx,&ry,&rz,READ_MC_CONF,0);
    /* INITIALIZE CONFIGURATION FILE */
    //MC_initialize_save_configs(mc_n, ini);
    MC_initialize_energies(mc_n, rx, ry, rz);
    MC_open_analysis_files(init);
  } else {
    printf("Arguments expected:\n[initial conf] [number of confs]\n");
    exit(1);
  }
  /* before running */
  printf("Analyze %d MC steps.\nStarting from iteration %d\n", mc_block, mc_ini);
  /******************/
  
  nt_n=mc_n/N_PARTS_PER_NT;
  for(i=0;i<mc_ini;++i)
    MC_read_ith_configuration(&energy_t, &temperature);
  for(i=mc_ini;i<mc_ini+mc_block;++i) {
    
    MC_read_ith_configuration(&energy_t, &temperature);
   
    //energ=MC_get_energy(nt_n, rx, ry, rz, 3);
    energ=energy_t;
    if(energ<minenerg){
      minenerg=energ;
      MC_min_energ_xyz_configuration(mc_n, rx, ry, rz, minenerg, temperature, init);
      MC_write_contact_list(nt_n, rx, ry, rz, init);
    }
  }
  MC_close_analysis_files();
#ifdef MPIMC
  MPI_Finalize();
#endif
  
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

void MC_read_ith_configuration(double *energy_t, double *temperature){
  int i,d, nt;
  int glp;
  double tempp[DIM];
  double dtempp[DIM];
  int nt_n=mc_n/N_PARTS_PER_NT;
  
  fread(energy_t,sizeof(double),1,md_configs);
  fread(temperature,sizeof(double),1,md_configs);
  //printf("energy = %lf  , temperature = %lf \n", *energy_t, *temperature);
  for(i=0;i<mc_n;i++){
    fread(dtempp,sizeof(double),DIM,md_configs);
    rx[i]=dtempp[0];ry[i]=dtempp[1];rz[i]=dtempp[2];
    //printf("%d %lf %lf %lf\n", i, rx[i], ry[i], rz[i]);
  }
  
  for(nt=0;nt<nt_n;nt++){
    //GLYC
    fread(&glp,sizeof(int),1,md_configs);
    mc_glyc[nt]=glp;
    //PUCK    
    fread(&glp,sizeof(int),1,md_configs);
    mc_puck[nt]=glp;
  }
}

void MC_close_analysis_files(){
  fclose(md_configs);
}
