#ifdef MPIMC
#include <mpi.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "mc.h"
#include "sa_impl.h"

#define LX 30
#define LY 30
#define LZ 30
#define NSOLUTE 868
#define MC_ITER 20000
#define EVERYMILLION 500000


/****************/

int main(int argc, char **argv) {
  int i,j;
  int rand_a=-1;
  double lx, ly, lz;
  char checkpoint[256];
  /* MC variables */
  int mc_n,nt_n;
  double *rx,*ry,*rz;
  double SA_DATA[ANN_NPARAMS];
  int sa_read_flag=0;
  double d_energ;
  int mc_iter, mc_i0, tempinit=0;
  int openflag=-1;
  double energy_t;
  int eff_iter;
  /* initialization */
  
  /*MPI STUFF*/
  int mpi_count, mpi_id=0, job_id, argmax;
#ifdef MPIMC
  argmax=2;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_count);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
#else
  int opt;
  while((opt = getopt(argc, argv, "hi:r")) != -1)  
    {  
      switch(opt)  
        {  
	case 'i':
	  mpi_id=atoi((char *)optarg);
	  //printf("Job_id = %d\n",mpi_id);
	  break;
	case 'h':
	  printf("Usage: SPQR_SA [-i job_id] [-r]\nRemember that params.pms and pdb_inits, with a proper initial condition must be present in the simulation directory.\n");
	  exit(ERR_INPUT);
	  break;
	case 'r':
	  sa_read_flag=1;
	  break;
	case '?':
	  if (optopt == 'i')
	    printf ("Option -%c requires an argument.\n", optopt);
	  else if (optopt == 'r')
	    printf ("Option -%c must be entered if the checkpoint contains previous parameters.\n", optopt);
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
  printf("Launching from %d \n" , mpi_id);
  openflag=MC_detect_initial_condition(mpi_id);
  
#ifdef MPIMC
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  tempinit=MC_initialize(&mc_n, &rx, &ry, &rz, &mc_iter, &rand_a, mpi_id, openflag, sa_read_flag, SA_DATA);
#ifdef MPIMC
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  nt_n=mc_n/N_PARTS_PER_NT;
  mc_i0=0;
  if(sa_read_flag){
    mc_i0=(int)(tempinit-floor(tempinit/mc_iter)*mc_iter);
  }
  energy_t=MC_get_energy(nt_n, rx, ry, rz, 3);
  
  /* ANNEALING */
  double sa_tmax=1.0, sa_tfac=1.0, sa_tmin=1, sa_sfac=1.0;
  int sa_NT=1;
  int ann_step=0;
  int sa_ini=0;
  double sa_temp;
  double sa_prev_energ, sa_this_energ, sa_temp_energ;
  double smc_nt_xyz, smc_ph_xyz, smc_nt_ang;
  int sa_resc_times=0;
  int sa_flag=1;
  int sa_mark=-1;
  //we read from params.spqr or update the data found in the checkpoint file
  SA_read_params(&sa_tmax, &sa_tmin, &sa_tfac, &sa_ini, &sa_NT, &sa_prev_energ, &sa_sfac, &sa_resc_times, mpi_id);
  if(sa_read_flag){
    SA_arr_to_params(&sa_tmax, &sa_tmin, &sa_tfac, &sa_ini, &sa_NT, &sa_prev_energ, &sa_sfac, &sa_resc_times, SA_DATA);
    printf("Reading SA parameters from checkpoint.\n");
  }
#ifdef MPIMC
  if(mpi_id==0)
#endif
    {
      printf("/******* SIMULATED ANNEALING **********/\n");
      printf("Starting from anneal step %d\n", sa_ini);
      printf("Maximum temperature:                %lf \n", sa_tmax);
      printf("Minimum temperature:                %lf \n", sa_tmin);
      printf("Temperature factor:                 %lf \n", sa_tfac);
      printf("Step factor:                        %lf \n", sa_sfac);
      printf("Number of iterations:                %d \n", sa_NT);
      printf("Previous energy:                     %lf \n", sa_prev_energ);
      printf("Each run is %d MC steps. Saving every %d trials.\n", mc_iter, mc_traj_steps);
    }
  
  printf("Processor (job) %d starting from time %d\n", mpi_id, mc_i0);
  sa_temp=sa_tmax;
  SA_init_mc_trials(&smc_nt_xyz, &smc_ph_xyz, &smc_nt_ang, MC_NT_XYZ*sa_pow(sa_sfac,sa_resc_times),MC_PH_XYZ*sa_pow(sa_sfac,sa_resc_times),MC_NT_ANGLE*sa_pow(sa_sfac,sa_resc_times));

#ifdef MPIMC
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  //INITIALIZE SIMULATION ENERGY : PREVIOUS ENERGY IS READ, OR INITIALIZED IN CASE IT IS THE FIRST STEP
  energy_t=MC_get_energy(nt_n, rx, ry, rz, 3);
  if(sa_ini==0)
    sa_prev_energ=energy_t;
  printf("JOB %d , INITIAL ENERGY IS %lf\n", mpi_id, energy_t);
#ifdef ERMSDR
  //energy_t+=0.5*ERMSD_PREF*ERMSD_SQ;
  energy_t+=ERMSD_ENERG;
#endif
  for(ann_step=sa_ini;ann_step<sa_NT;ann_step++){
    printf("Step %d on processor %d, temperature set to %lf\n", ann_step, mpi_id, sa_temp);
    SA_set_temp(sa_temp);
    SA_set_mc_trials(smc_nt_xyz, smc_ph_xyz, smc_nt_ang);
    sa_temp_energ=energy_t;
#ifdef ERMSDR
    sa_temp_energ-=ERMSD_ENERG;
#endif
    sa_this_energ=sa_temp_energ;
    eff_iter=mc_iter;
    if(sa_temp<=6.0) {
      eff_iter=mc_iter/10;
      if(sa_mark==-1) sa_mark=ann_step;
    }
    for(i=mc_i0;i<eff_iter;++i) {
      /* do your MC here */
      d_energ=MC_integrate(mc_n, &rx, &ry, &rz);
      energy_t+=d_energ;
      //WE DON'T CONSIDER THE PULLING IN THE ANNEALING ENERGY!
      sa_temp_energ=energy_t;
#ifdef ERMSDR
      sa_temp_energ-=ERMSD_ENERG;
      //0.5*ERMSD_PREF*ERMSD_SQ;
      if((i+1)%mc_traj_steps==0)
	MC_write_ermsd_obs(i+1,energy_t);
#endif
      //if(i==mc_i0+1) sa_this_energ=sa_temp_energ;
      if(sa_this_energ > sa_temp_energ) //find the minimum energy
	sa_this_energ = sa_temp_energ;
      if((i+1)%mc_traj_steps==0){
	MC_save_configuration(mc_n, rx, ry, rz, energy_t,(i+1)+(ann_step-sa_mark)*eff_iter + sa_mark*mc_iter);
      }
      
      if(((i+1)%(mc_chkp_steps))==0 && i<eff_iter){
	SA_params_to_arr(sa_temp, sa_tmin, sa_tfac, ann_step, sa_NT, sa_prev_energ, sa_sfac, sa_resc_times, SA_DATA);
	MC_save_checkpoint(mc_n, rx, ry, rz, (i+1)+(ann_step-sa_mark)*eff_iter + sa_mark*mc_iter, energy_t, mpi_id, SA_DATA); //this writes a binary checkpoint
      }
    }
    printf("Step %d on %d, energy is %lf . It was %lf\n", ann_step, mpi_id, sa_this_energ, sa_prev_energ);
    if(sa_temp <sa_tmin){
      sa_temp=0;
    }
    else{      
      if(sa_this_energ>=sa_prev_energ){ // energy is not dropping - reset temperature
	sa_temp=sa_temp*sa_tfac;
	if(sa_temp<1.0){
	  sa_resc_times++;
	  SA_rescale_mc_trials(&smc_nt_xyz, &smc_ph_xyz, &smc_nt_ang, sa_sfac);
	}
      }
    }
    sa_prev_energ=sa_this_energ;
    SA_params_to_arr(sa_temp, sa_tmin, sa_tfac, ann_step+1, sa_NT, sa_prev_energ, sa_sfac, sa_resc_times, SA_DATA);
    if((eff_iter+(ann_step-sa_mark)*eff_iter + sa_mark*mc_iter) %(mc_chkp_steps)==0)
      MC_save_checkpoint(mc_n, rx, ry, rz, eff_iter+(ann_step-sa_mark)*eff_iter + sa_mark*mc_iter, energy_t, mpi_id, SA_DATA); //this writes a binary checkpoint
    //printf("%d %d  %d %d\n", mc_iter, eff_iter, sa_mark,eff_iter+(ann_step-sa_mark)*eff_iter + sa_mark*mc_iter );
    //SA_save_params(sa_temp, sa_tmin, sa_tfac, ann_step+1, sa_NT, sa_prev_energ, sa_sfac, sa_resc_times, mpi_id);
  }
  printf("Simulation %d ended properly. Checkpoint saved with temperature %lf and minimum energy %lf .\n", mpi_id, sa_temp, sa_this_energ);
  //MC_write_pdb("final", mc_n, rx, ry, rz, sa_this_energ, mpi_id);
  MC_write_pdb("final", mc_n, rx, ry, rz, energy_t, mpi_id);
  MC_save_checkpoint(mc_n, rx, ry, rz, -1, energy_t, mpi_id, SA_DATA); //this writes a binary checkpoint
  MC_end(mc_n, rx, ry, rz, i, energy_t, mpi_id);

#ifdef MPIMC
  MPI_Finalize();
#endif
  return 0;
}
