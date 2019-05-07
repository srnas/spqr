#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mc.h"
//#include "st_impl.h"

#define LX 30
#define LY 30
#define LZ 30
#define NSOLUTE 868
#define MC_ITER 20000
#define EVERYMILLION 1000000
//#define EVERYMILLION 100000
#define EVERYTHOUSAND 10000


/****************/

int main(int argc, char **argv) {
  int i,j;
  int rand_a;
  double lx, ly, lz;
  char checkpoint[256];
  /* MC variables */
  int mc_n;
  double *rx,*ry,*rz;
  double d_energ, energy_t;
;
  int mc_iter, mc_ini;
  /* initialization */
  
  /*MPI STUFF*/
  int mpi_count, mpi_id;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_count);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
  
  /***********/
  
  mc_n=(int)NSOLUTE;
  printf("Launching from %d \n" , mpi_id);
  int nt_n;
   /* SIMULATED TEMPERING */
  int ST_ind, ST_ind_trial, ST_N;
  double ST_DELTA, *ST_g, *ST_Ts, ST_C;
  double ST_expfac;
  /************************/
  
  if (argc==2){
    MC_read_nsolute(&mc_n, mpi_id,NULL);
    mc_ini=0;
    MC_read_params(&lx, &ly, &lz, &mc_iter, &rand_a, mpi_id);
    /* GENERAL INIT */
    MC_initialize(mc_n, &rx, &ry, &rz, lx, ly, lz, READ_MC_CONF, rand_a, mc_ini, mpi_id);//this initializes wrong the name of the confs.mc for reinitialized simulations!
    /* SIMULATED TEMPERING */
    /* each processor reads its ST parameters */
    ST_read_parameters(mpi_id, mpi_count, &ST_ind, &ST_DELTA, &ST_C, &mc_ini, &ST_g, &ST_Ts , &ST_N);
    MPI_Barrier(MPI_COMM_WORLD);
    fflush(stdout);
    /************************/

    
  }  else {
    printf("Default case not implemented.\n");
    exit(1);
  }
  
  
  
  //exit(0);  
  /* before running */
  
  /******************/
  
  nt_n=mc_n/N_PARTS_PER_NT;
  if(mpi_id==0)
    printf("Run %d MC steps. Saving trajectory and checkpointing every %d and %d swaps.\n", mc_iter, mc_traj_steps, mc_chkp_steps);
  
  energy_t=MC_get_energy(nt_n, rx, ry, rz, 3);
  double st_energy=energy_t;
  MPI_Barrier(MPI_COMM_WORLD);
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  
  /***************************************/
  FILE *TEMPFILE;
  char tfname[256];
  sprintf(tfname, "temps.p%02d.dat", mpi_id);
  TEMPFILE=fopen(tfname, "w");
  /***************************************/
  
  printf("Processor %d, INITIAL TEMPERATURE:  %lf   INITIAL ENERGY IS %lf\n", mpi_id, ST_Ts[ST_ind], energy_t);
  fflush(stdout);
  for(i=0;i<mc_iter;++i) {
    //first, we update the temperature
    if((i%mc_traj_steps)==0)
      ST_update_temperature(ST_Ts[ST_ind]);
    //for(j=0;j<mc_chk_freq;++j) {
    /* do your MC here */
      d_energ=MC_integrate(mc_n, &rx, &ry, &rz);
      energy_t+=d_energ;
#ifdef ERMSDR
      if((i+1)%mc_traj_steps==0)
	MC_write_ermsd_obs(i+1,energy_t);
#endif
      st_energy=energy_t-0.5*ERMSD_PREF*ERMSD_SQ;
      //attempt to change the temperature
      if(rand_d(1.0)<0.5) ST_ind_trial=ST_ind+1;
      else ST_ind_trial=ST_ind-1;
      
      if(ST_ind_trial<ST_N && ST_ind_trial>=0){
	ST_expfac=(1.0/ST_Ts[ST_ind_trial]-1.0/ST_Ts[ST_ind])*st_energy-(ST_g[ST_ind_trial]-ST_g[ST_ind]);
	if(rand_d(1.0) < exp(-ST_expfac)){
	  ST_ind=ST_ind_trial;
	}
      }
      else 
	ST_ind_trial=ST_ind; // for safety
      
      
      //ST_g[ST_ind]-=ST_DELTA*(1.0/(1.0+ST_C*sqrt((double)i)));
      ST_g[ST_ind]-=ST_DELTA*(1.0/(1.0+ST_C*((double)i)));
      //temperature is updated at the beginning of the integration loop
      
      /***************************************/
      if(i%mc_traj_steps==0) {
	//fprintf(TEMPFILE, "%d %lf %lf   %lf \n", i, mc_target_temp, energy_t,MC_get_energy(nt_n, rx, ry, rz, 3));
	fprintf(TEMPFILE, "%d %lf %lf\n", i, mc_target_temp, st_energy);
      }
      /***************************************/
      
      if(((i+1)%mc_traj_steps)==0)
	MC_save_configuration(mc_n, rx, ry, rz, energy_t); //this writes in binary file
      if(((i+1)%(mc_chkp_steps))==0){
	//MC_save_current_configuration(mc_n, rx, ry, rz, i+1, energy_t, mpi_id); //this writes xyz and glp files which allow restarting
	MC_save_checkpoint(mc_n, rx, ry, rz, (i+1), energy_t, mpi_id,NULL); //this writes a binary checkpoint
      }
  }
  ST_write_parameters(mpi_id, mpi_count, ST_ind, ST_DELTA, ST_C, i, ST_g, ST_Ts , ST_N);
  printf("Simulation %d ended properly.\n", mpi_id);
  MC_end(mc_n, rx, ry, rz, i, energy_t, mpi_id);
  
  /***************************************/
  fclose(TEMPFILE);
  /***************************************/
  MPI_Finalize();
  return 0;
}
