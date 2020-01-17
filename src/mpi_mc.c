#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
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
  int mpi_id=0,argmax;
#ifdef MPIMC
  int mpi_count;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_count);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
  #else
  int opt;
  while((opt = getopt(argc, argv, "hi:")) != -1)  
    {  
      switch(opt)  
        {  
	case 'i':
	  mpi_id=atoi((char *)optarg);
	  //printf("Job_id = %d\n",mpi_id);
	  break;
	case 'h':
	  printf("Usage: ./SPQR_MC [-i job_id]\nRemember that params.pms and pdb_inits, with a proper initial condition must be present in the simulation directory.\n");
	  exit(ERR_INPUT);
	  break;
	case '?':
	  if (optopt == 'i')
	    printf ("Option -%c requires an argument.\n", optopt);
	  else if (isprint (optopt))
	    printf ("Unknown option `-%c'.\n", optopt);
	  else
	    printf ("Unknown option character `\\x%x'.\n",   optopt);
	  exit(ERR_INPUT);
	  break;
	}
    }  
#endif
  /***********/
  
  mc_n=(int)NSOLUTE;
  
#ifdef MPIMC
  printf("Launching from %d \n" , mpi_id);
#else
  printf("Launching from single processor. Job_id = %d \n", mpi_id);
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
  int n;
  nt_n=mc_n/N_PARTS_PER_NT;
  //IF -R IS GIVEN, Mc_i0=0
  if(argc>2 && !strcmp(argv[2],"-r")) mc_i0=tempinit;
  else mc_i0=0;
  
  //if(mpi_id==0)
  printf("Simulation %d . Run %d MC steps. Saving trajectory and checkpointing every %d and %d swaps.\n", mpi_id, mc_iter, mc_traj_steps, mc_chkp_steps);
  energy_t=MC_get_energy(nt_n, rx, ry, rz, 3);
  printf("INITIAL ENERGY IS %lf at step %d\n", energy_t, mc_i0);
#ifdef ERMSDR
  //energy_t+=0.5*ERMSD_PREF*ERMSD_SQ;
  energy_t+=ERMSD_ENERG;
#endif
  for(i=mc_i0;i<mc_iter;++i) {
    /* do your MC here */
    d_energ=MC_integrate(mc_n, &rx, &ry, &rz);
    energy_t+=d_energ;
    
#ifdef ERMSDR
    if((i+1)%mc_traj_steps==0)
      MC_write_ermsd_obs(i+1,energy_t);
#endif
    if((i+1)%mc_traj_steps==0)
      MC_save_configuration(mc_n, rx, ry, rz, energy_t,i+1); //this writes in binary file and possibly the pdb output
    if(((i+1)%(mc_chkp_steps))==0){
      MC_save_checkpoint(mc_n, rx, ry, rz, (i+1), energy_t, mpi_id,NULL); //this writes a binary checkpoint
    }
  }
  
  printf("Simulation %d ended properly.\n", mpi_id);
  MC_write_pdb("final", mc_n, rx, ry, rz, energy_t, mpi_id);
  MC_save_checkpoint(mc_n, rx, ry, rz, -1, energy_t, mpi_id,NULL); //this writes a binary checkpoint
  MC_end(mc_n, rx, ry, rz, i, energy_t, mpi_id);
#ifdef ERMSDR
  fclose(ermsd_obs);
#endif
  
#ifdef MPIMC
  MPI_Finalize();
#else
  return 0;
#endif
}
