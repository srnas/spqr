#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mc.h"

#ifdef MPIMC
#include <mpi.h>
#endif

/* #define LX 30 */
/* #define LY 30 */
/* #define LZ 30 */
#define NSOLUTE 868
#define MC_ITER 20000
#define EVERYMILLION 500000


/****************/

int main(int argc, char **argv) {
  int i,j;
  int tempinit;
  int rand_a=0;
  //double lx, ly, lz;
  char checkpoint[256];
  /* MC variables */
  int mc_n, nt_n;
  double *rx,*ry,*rz;
  double d_energ, energy_t;
  int mc_iter=MC_ITER, mc_i0=0;
  int openflag=-1;
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
  openflag=MC_detect_initial_condition(mpi_id);
#ifdef MPIMC
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  /**********initialization************/
  tempinit=MC_initialize(&mc_n, &rx, &ry, &rz, &mc_iter, &rand_a, mpi_id, openflag, NULL);
  /************************************/
  
#ifdef MPIMC
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  /* before running */
  
  /******************/
  int n, face;
  nt_n=mc_n/N_PARTS_PER_NT;
  //IF -R IS GIVEN, Mc_i0=0
  if(argc>2 && !strcmp(argv[2],"-r")) mc_i0=tempinit;
  else mc_i0=0;
  
  if(mpi_id==0)
    printf("Run %d MC steps. Saving trajectory and checkpointing every %d and %d swaps.\n", mc_iter, mc_traj_steps, mc_chkp_steps);
  double jcenergy;
  energy_t=MC_get_energy(nt_n, rx, ry, rz, 3);
  printf("INITIAL ENERGY IS %lf at step %d\n", energy_t, mc_i0);
#ifdef ERMSDR
  energy_t+=0.5*ERMSD_PREF*ERMSD_SQ;
#endif
  for(i=mc_i0;i<mc_iter;++i) {
    /* do your MC here */
    d_energ=MC_integrate(mc_n, &rx, &ry, &rz);
    energy_t+=d_energ;
    
    /* jcenergy=MC_get_energy(nt_n, rx, ry, rz, 3)+0.5*ERMSD_PREF*ERMSD_SQ; */
    /* printf("Step %d  %d  ,   acc = %lf  get_energ = %lf\n", i,j,energy_t, jcenergy); */
    /* if(fabs(jcenergy-energy_t)>0.01){ */
    /*   printf("STOP ! %d  %d  ,   acc = %lf  get_energ = %lf\n", i,j,energy_t, jcenergy); */
    /*   exit(1); */
    /* } */
    
#ifdef ERMSDR
    //energy_t=MC_get_energy(nt_n, rx, ry, rz, 3);
    //ERMSD_SQ=get_current_ermsd(&rx, &ry, &rz, nt_n);
    //energy_t+=0.5*ERMSD_PREF*ERMSD_SQ;
    if((i+1)%mc_traj_steps==0)
      MC_write_ermsd_obs(i+1,energy_t);
#endif
    if((i+1)%mc_traj_steps==0)
      MC_save_configuration(mc_n, rx, ry, rz, energy_t); //this writes in binary file
    if(((i+1)%(mc_chkp_steps))==0){
      //MC_save_current_configuration(mc_n, rx, ry, rz, i*j, i+1, energy_t, mpi_id); //this writes xyz and glp files which allow restarting
      MC_save_checkpoint(mc_n, rx, ry, rz, (i+1), energy_t, mpi_id,NULL); //this writes a binary checkpoint
    }
  }
  
  printf("Simulation %d ended properly.\n", mpi_id);
  MC_write_pdb("final", mc_n, rx, ry, rz, energy_t, mpi_id);
  MC_save_checkpoint(mc_n, rx, ry, rz, -1, energy_t, mpi_id,NULL); //this writes a binary checkpoint
  MC_end(mc_n, rx, ry, rz, i, energy_t, mpi_id);
  //printf("Checkpoints written.\n");
#ifdef MPIMC
  MPI_Finalize();
#endif
  return 0;
}
