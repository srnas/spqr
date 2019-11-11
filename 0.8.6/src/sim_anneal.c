#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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
  int rand_a;
  double lx, ly, lz;
  char checkpoint[256];
  /* MC variables */
  int mc_n;
  double *rx,*ry,*rz;
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
  
  MPI_Barrier(MPI_COMM_WORLD);
  /* before running */
  
  /******************/
  
  nt_n=mc_n/N_PARTS_PER_NT;
  
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
  double bia_pref=1;
  double temp_e;
  //SA_init_mc_trials(&smc_nt_xyz, &smc_ph_xyz, &smc_nt_ang);
  SA_read_params(&sa_tmax, &sa_tmin, &sa_tfac, &sa_ini, &sa_NT, &sa_prev_energ, &sa_sfac, &sa_resc_times, mpi_id);
  
  if(mpi_id==0){
    printf("/******* SIMULATED ANNEALING **********/\n");
    printf("Starting from anneal step %d\n", sa_ini);
    printf("Maximum temperature:                %lf \n", sa_tmax);
    printf("Minimum temperature:                %lf \n", sa_tmin);
    printf("Temperature factor:                 %lf \n", sa_tfac);
    printf("Step factor:                        %lf \n", sa_sfac);
    printf("Number of iterations:                %d \n", sa_NT);
    printf("Initial energy:                     %lf \n", sa_prev_energ);
  }
  if(mpi_id==0)
    printf("Each run is %d MC steps. Saving every %d trials.\n", mc_iter*mc_chk_freq, mc_chk_freq);
  sa_temp=sa_tmax;
  SA_init_mc_trials(&smc_nt_xyz, &smc_ph_xyz, &smc_nt_ang, MC_NT_XYZ*sa_pow(sa_sfac,sa_resc_times),MC_PH_XYZ*sa_pow(sa_sfac,sa_resc_times),MC_NT_ANGLE*sa_pow(sa_sfac,sa_resc_times));
  MPI_Barrier(MPI_COMM_WORLD);
  if(mpi_id==0){
#ifdef NEW_BIA
    bia_pref=mc_target_temp;
    mc_target_temp=1.0;
#endif
    temp_e=MC_get_energy(nt_n, rx, ry, rz, 3);
#ifdef NEW_BIA
    mc_target_temp=bia_pref;
    /* bia_pref=1; */
    /* if(mc_target_temp>1) bia_pref=mc_target_temp; */
    /* temp_e/=bia_pref; */
#endif
    printf("INITIAL ENERGY IS %lf\n", temp_e);
  }
  for(ann_step=sa_ini;ann_step<sa_NT;ann_step++){
    printf("Step %d on %d, temperature set to %lf\n", ann_step, mpi_id, sa_temp);
    SA_set_temp(sa_temp);
    SA_set_mc_trials(smc_nt_xyz, smc_ph_xyz, smc_nt_ang);
    sa_temp_energ=0;
    if(ann_step==0){
#ifdef NEW_BIA
      bia_pref=mc_target_temp;
      mc_target_temp=1.0;
#endif
      sa_prev_energ=MC_get_energy(nt_n, rx, ry, rz, 3);
#ifdef NEW_BIA
      mc_target_temp=bia_pref;
#endif
    }
    for(i=mc_ini;i<mc_iter;++i) {
      for(j=0;j<mc_chk_freq;++j) {
	/* do your MC here */
	MC_integrate(mc_n, &rx, &ry, &rz);
      }
#ifdef NEW_BIA
      bia_pref=mc_target_temp;
      mc_target_temp=1.0;
#endif   
      sa_temp_energ=MC_get_energy(nt_n, rx, ry, rz, 3);
#ifdef NEW_BIA
      mc_target_temp=bia_pref;
#endif
      if(i==0) sa_this_energ=sa_temp_energ;
      if(sa_this_energ > sa_temp_energ) //find the minimum energy
	sa_this_energ = sa_temp_energ;
      
      MC_save_configuration(mc_n, rx, ry, rz);
    }
    MC_save_current_configuration(mc_n, rx, ry, rz, i*j, ann_step*mc_iter, mpi_id);
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
    
    
    
    //SA_save_params(ann_step,sa_temp, sa_tmin, sa_tfac, ann_step+1, sa_NT, sa_prev_energ, sa_sfac, sa_resc_times, mpi_id);
    SA_save_params(sa_temp, sa_tmin, sa_tfac, ann_step+1, sa_NT, sa_prev_energ, sa_sfac, sa_resc_times, mpi_id);
  }
  
  printf("Simulation %d ended properly.\n", mpi_id);
  MC_end(mc_n, rx, ry, rz, i, mpi_id);
  MPI_Finalize();
  return 0;
}
