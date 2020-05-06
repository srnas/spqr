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
  int mc_read_flag=0;
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
  
  int pt_sel,ptid,pt_next,id_sel, id_next, id_temp;
  double energ_sel, energ_next, chosen_temp, next_temp, beta_sel, beta_next, PT_expfac, PT_fac, dtemp;
  FILE *PT_FILE;
  PT_ind=mpi_id;
  PT_N=mpi_count;
  int *PT_list=malloc(sizeof(int)*PT_N);
  for(i=0;i<PT_N;i++)
    PT_list[i]=i;
  /************************/


  openflag=MC_detect_initial_condition(mpi_id);
  MPI_Barrier(MPI_COMM_WORLD);
  tempinit=MC_initialize(&mc_n, &rx, &ry, &rz, &mc_iter, &rand_a, mpi_id, openflag, mc_read_flag, NULL);
  
  /* PARALLEL TEMPERING */
  /* each processor reads its PT parameters */
  PT_read_params(mpi_id, mpi_count,&PT_FILE);
  MPI_Barrier(MPI_COMM_WORLD);
  PT_freq=2;
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
      pt_sel=0;
      if(mpi_id==mpi_root){
	//only one processor sets this
	pt_sel=rand_i(PT_N-1);
	pt_next=pt_sel+1;
	for(ptid=0;ptid<PT_N;ptid++){
	  if(PT_list[ptid]==pt_sel) id_sel=ptid;
	  if(PT_list[ptid]==pt_next) id_next=ptid;
	}
      }
      MPI_Bcast(&id_sel, 1, MPI_INT, mpi_root, MPI_COMM_WORLD);
      MPI_Bcast(&id_next, 1, MPI_INT, mpi_root, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      //we copy the temperatures and energies of processors pt_sel and pt_sel+1
      if(mpi_id==mpi_root){
	//we have to treat separately the case where mpi_root == id_sel or id_next
	if(mpi_root != id_sel){
	  MPI_Recv(&chosen_temp,1,MPI_DOUBLE,id_sel  ,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  MPI_Recv(&energ_sel  ,1,MPI_DOUBLE,id_sel  ,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
	else{
	  chosen_temp=mc_target_temp;
	  energ_sel=energy_t;
	}
	if(mpi_root != id_next){
	  MPI_Recv(&next_temp  ,1,MPI_DOUBLE,id_next,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  MPI_Recv(&energ_next ,1,MPI_DOUBLE,id_next,5,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
	else{
	  next_temp=mc_target_temp;
	  energ_next=energy_t;
	}
	
	beta_sel=1.0/chosen_temp;
	beta_next=1.0/next_temp;
	PT_expfac=(beta_sel-beta_next)*(energ_next-energ_sel);
	PT_fac=exp(-PT_expfac);
	if(rand_d(1.0)<PT_fac){
	  //we perform the swap of temperatures
	  dtemp=chosen_temp;
	  chosen_temp=next_temp;
	  next_temp=dtemp;
	  
	  //and perform the swap of indexes
	  //if(PT_list[ptid]==pt_sel) id_sel=ptid;
	  //if(PT_list[ptid]==pt_next) id_next=ptid;
	  PT_list[id_sel]=pt_next;
	  PT_list[id_next]=pt_sel;
	  
	}
	if(mpi_root != id_sel)
	  MPI_Send(&chosen_temp,1,MPI_DOUBLE,id_sel  ,2,MPI_COMM_WORLD);
	else
	  mc_target_temp=chosen_temp;
	if(mpi_root != id_next)
	  MPI_Send(&next_temp  ,1,MPI_DOUBLE,id_next,3,MPI_COMM_WORLD);
	else
	  mc_target_temp=next_temp;
	//temperature is updated here, immediately
	
      }
      else if (mpi_id == id_sel) {
	if(mpi_root != id_sel){
	  MPI_Send(&mc_target_temp,1,MPI_DOUBLE,mpi_root,0,MPI_COMM_WORLD);
	  MPI_Send(&energy_t      ,1,MPI_DOUBLE,mpi_root,4,MPI_COMM_WORLD);
	  MPI_Recv(&mc_target_temp,1,MPI_DOUBLE,mpi_root,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
      }
      else if (mpi_id == id_next) {
	if(mpi_root != id_next){
	  MPI_Send(&mc_target_temp,1,MPI_DOUBLE,mpi_root,1,MPI_COMM_WORLD);
	  MPI_Send(&energy_t      ,1,MPI_DOUBLE,mpi_root,5,MPI_COMM_WORLD);
	  MPI_Recv(&mc_target_temp,1,MPI_DOUBLE,mpi_root,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
      }
      //finally, for safety, we export the array of indexes
      MPI_Bcast(PT_list,PT_N,MPI_INT,mpi_root,MPI_COMM_WORLD);
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
