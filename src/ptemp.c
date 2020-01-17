#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "mc.h"
#include "pt_impl.h"

#define NSOLUTE 868
#define MC_ITER 20000
#define EVERYMILLION 1000000


/****************/

int main(int argc, char **argv) {
  int i,j;
  int tempinit;
  int rand_a=0;
  double lx, ly, lz;
  char checkpoint[256];
  /* MC variables */
  int mc_n, nt_n;
  double *rx,*ry,*rz;
  double d_energ, energy_t;
  int mpi_root=0;
  int openflag;
  int mc_iter, mc_ini;
  /* initialization */
  
  /*MPI STUFF*/
  int mpi_count, mpi_id;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_count);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
  
  /***********/
  mc_n=(int)NSOLUTE;
  
  
   /* PARALLEL TEMPERING */
  int PT_ind, PT_freq, PT_N;
  
  int pt_sel;
  double energ_sel, energ_next, chosen_temp, next_temp, beta_sel, beta_next, PT_expfac, PT_fac, dtemp;
  FILE *PT_FILE;
  PT_ind=mpi_id;
  PT_N=mpi_count;
  
  /************************/


  openflag=MC_detect_initial_condition(mpi_id);
  MPI_Barrier(MPI_COMM_WORLD);
  tempinit=MC_initialize(&mc_n, &rx, &ry, &rz, &mc_iter, &rand_a, mpi_id, openflag, NULL);
  
  /* PARALLEL TEMPERING */
  /* each processor reads its PT parameters */
  PT_read_params(mpi_id, mpi_count,&PT_FILE);
  MPI_Barrier(MPI_COMM_WORLD);
  PT_freq=10*mc_traj_steps;
  fflush(stdout);
  /************************/

  /* before running */
  
  /******************/
  
  nt_n=mc_n/N_PARTS_PER_NT;
  if(mpi_id==mpi_root)
    printf("PARALLEL TEMPERING:\nRun %d MC steps. Saving trajectory and checkpointing every %d and %d swaps.\n", mc_iter, mc_traj_steps, mc_chkp_steps);
  
  energy_t=MC_get_energy(nt_n, rx, ry, rz, 3);
  double pt_energy=energy_t;
  MPI_Barrier(MPI_COMM_WORLD);
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  
  /***************************************/
  printf("Processor %d, INITIAL TEMPERATURE:  %lf   INITIAL ENERGY IS %lf\n", mpi_id, mc_target_temp, energy_t);
  fflush(stdout);
  

  /**** INTEGRATION LOOP *****/
  for(i=0;i<mc_iter;++i) {
    /* do your MC here */
    d_energ=MC_integrate(mc_n, &rx, &ry, &rz);
    energy_t+=d_energ;
    
    pt_energy=energy_t;

    /******SAVE ENERGY AND TEMPERATURE*******/
    
    /**** WRITE CONFS AND CHECKPOINTS ****/
    if(((i+1)%mc_traj_steps)==0){
      MC_save_configuration(mc_n, rx, ry, rz, energy_t,i+1); //this writes in binary file
      PT_write(i+1, mc_target_temp, energy_t,&PT_FILE);
    }
    if(((i+1)%(mc_chkp_steps))==0){
      MC_save_checkpoint(mc_n, rx, ry, rz, (i+1), energy_t, mpi_id,NULL); //this writes a binary checkpoint
    }
    
    /** HERE SWAP THE TEMPERATURES **/
    //select randomly the processor
    if(((i+1)%PT_freq)==0){
      printf("performing swap!\n");
      pt_sel=0;
      if(mpi_id==mpi_root){
	//only one processor sets this
	pt_sel=rand_i(PT_N-1);
	printf("SELECTED %d\n", pt_sel);
      }
      MPI_Bcast(&pt_sel, 1, MPI_INT, mpi_root, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      //we copy the temperatures and energies of processors pt_sel and pt_sel+1
      if(mpi_id==mpi_root){
	//we have to treat separately the case where mpi_root == pt_sel
	if(mpi_root != pt_sel){
	  MPI_Recv(&chosen_temp,1,MPI_DOUBLE,pt_sel  ,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  MPI_Recv(&energ_sel  ,1,MPI_DOUBLE,pt_sel  ,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
	else{
	  chosen_temp=mc_target_temp;
	  energ_sel=energy_t;
	}
	MPI_Recv(&next_temp  ,1,MPI_DOUBLE,pt_sel+1,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Recv(&energ_next ,1,MPI_DOUBLE,pt_sel+1,5,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	beta_sel=1.0/chosen_temp;
	beta_next=1.0/next_temp;
	PT_expfac=(beta_sel-beta_next)*(energ_next-energ_sel);
	PT_fac=exp(PT_expfac);
	if(rand_d(1.0)<PT_fac){
	  //we perform the swap of temperatures
	  dtemp=chosen_temp;
	  chosen_temp=next_temp;
	  next_temp=dtemp;
	}
	if(mpi_root != pt_sel)
	  MPI_Send(&chosen_temp,1,MPI_DOUBLE,pt_sel  ,2,MPI_COMM_WORLD);
	else
	  mc_target_temp=chosen_temp;
	MPI_Send(&next_temp  ,1,MPI_DOUBLE,pt_sel+1,3,MPI_COMM_WORLD);
	//temperature is updated here, immediately
      }
      else if (mpi_id == pt_sel) {
	if(mpi_root != pt_sel){
	  MPI_Send(&mc_target_temp,1,MPI_DOUBLE,mpi_root,0,MPI_COMM_WORLD);
	  MPI_Send(&energy_t      ,1,MPI_DOUBLE,mpi_root,4,MPI_COMM_WORLD);
	  MPI_Recv(&mc_target_temp,1,MPI_DOUBLE,mpi_root,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
      }
      else if (mpi_id == pt_sel+1) {
	MPI_Send(&mc_target_temp,1,MPI_DOUBLE,mpi_root,1,MPI_COMM_WORLD);
	MPI_Send(&energy_t      ,1,MPI_DOUBLE,mpi_root,5,MPI_COMM_WORLD);
	MPI_Recv(&mc_target_temp,1,MPI_DOUBLE,mpi_root,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
  
  printf("Simulation %d ended properly.\n", mpi_id);
  MC_write_pdb("final", mc_n, rx, ry, rz, energy_t, mpi_id);
  MC_save_checkpoint(mc_n, rx, ry, rz, -1, energy_t, mpi_id,NULL); //this writes a binary checkpoint
  MC_end(mc_n, rx, ry, rz, i, energy_t, mpi_id);
  
  /***************************************/
  PT_end(&PT_FILE);
  /***************************************/
  MPI_Finalize();
}
