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


/****************/

int main(int argc, char **argv) {
  int i,j;
  int rand_a;
  double lx, ly, lz;
  char checkpoint[256];
  /* MC variables */
  int mc_n;
  double *rx,*ry,*rz;
  int mc_iter, mc_ini;
  /* initialization */
  
  /*MPI STUFF*/
  int mpi_id=0;
#ifdef MPIMC
  int mpi_count;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_count);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
#endif
  /***********/
  
  
  
  mc_n=(int)NSOLUTE;
#ifdef MPIMC
  printf("Launching from %d \n" , mpi_id);
#else
  printf("Launching from single processor \n");
#endif
  
  int nt_n;
  if (argc==2){
    MC_read_nsolute(&mc_n, mpi_id);
    mc_ini=0;
    MC_read_params(&lx, &ly, &lz, &mc_iter, &rand_a, mpi_id);
    
    /* MPI STUFF */
    //rand_a+=mpi_id;
    MC_initialize(mc_n, &rx, &ry, &rz, lx, ly, lz, READ_MC_CONF, rand_a, mc_ini, mpi_id);
  }  else {
    printf("Default case not implemented.\n");
    exit(1);
  }
  
  /* before running */
  if(mpi_id==0)
    printf("Run %d MC steps. Saving every %d trials.\n", mc_iter, mc_chk_freq);
  /******************/
  
  nt_n=mc_n/N_PARTS_PER_NT;
  for(i=mc_ini;i<mc_iter;++i) {
    for(j=0;j<mc_chk_freq;++j) {
      //printf("%d, %d\n", i,j);
      /* do your MC here */
      MC_integrate(mc_n, &rx, &ry, &rz);
    }
    
    MC_save_configuration(mc_n, rx, ry, rz);
    if(((i+1)%(EVERYMILLION))==0){
      MC_save_current_configuration(mc_n, rx, ry, rz, i*j, i+1, mpi_id);
    }
    //MC_eval_total_energy(nt_n, rx, ry, rz);
  }
  printf("Simulation %d ended properly.\n", mpi_id);
  MC_end(mc_n, rx, ry, rz, i, mpi_id);
#ifdef MPIMC
  MPI_Finalize();
#endif
  return 0;
}
