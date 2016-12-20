#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mc.h"



#define LX 30
#define LY 30
#define LZ 30
#define NSOLUTE 868
#define MC_ITER 20000
#define EVERYMILLION 500000


/****************/

int main(int argc, char **argv) {
  int i,j,ex;
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
  int mpi_count;
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
  if(mpi_id==0)
    printf("Run %d MC steps. Saving every %d trials.\n", mc_iter, mc_chk_freq);
  /******************/
  
  nt_n=mc_n/N_PARTS_PER_NT;
  
  
  /* PARALLEL TEMPERING */
  double *TEMP;
  //TEMP=(double *)malloc(sizeof(double)*mpi_count);
  PT_create_params(TEMP, mpi_count); //HERE WE ALLOCATE AND READ FILE, OR INIT THEM WITH DEFAULT VALUES
  //PT_set_temperatures(TEMP);
  PT_set_temp(TEMP[mpi_id]);
  double ENERG=MC_get_energy(nt_n, rx, ry, rz, 3);
  if(mpi_id==0){
    printf("/******* PARALLEL TEMPERING **********/\n");
    //printf("Maximum temperature:                %lf \n", sa_tmax);
    //printf("Minimum temperature:                %lf \n", sa_tmin);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  for(i=0;i<mpi_count;i++)
    printf("**  Proc %d\thas T=%lf  \t E=%lf\t**\n", i, TEMP[i], ENERG);
  MPI_Barrier(MPI_COMM_WORLD);
  
  //double pa_param=PT_init_params(mpi_id, mpi_count);
  double buff[2], buff1[2],buff2[2];
  int max_ex=100;
  double rand_ch, energ, ttemp;
  int rand_proc;

  
  
  for(ex=0;ex<max_ex;ex++){
    for(i=mc_ini;i<mc_iter;++i) {
      for(j=0;j<mc_chk_freq;++j) {
	/* do your MC here */
	//MC_integrate(mc_n, &rx, &ry, &rz, pa_param);
	MC_integrate(mc_n, &rx, &ry, &rz);
      }
      
      MC_save_configuration(mc_n, rx, ry, rz);
      if(((i+1)%(EVERYMILLION))==0){
	MC_save_current_configuration(mc_n, rx, ry, rz, i*j, i+1, mpi_id);
      }
    }
    energ=MC_get_energy(nt_n, rx, ry, rz, 3);
    buff[0]=energ;
    buff[1]=TEMP[mpi_id];
    
    MPI_Barrier(MPI_COMM_WORLD);
    if(mpi_id==0){
      //choose random replicas
      rand_proc=rand_i(mpi_count-1);
      //communicate rand_proc to everyone
      MPI_Bcast(&rand_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(mpi_id==rand_proc){
      MPI_Send(buff, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    if(mpi_id==rand_proc){
      MPI_Send(&buff, 2, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    if(mpi_id==0){
      MPI_Recv(&buff1, 2, MPI_DOUBLE, rand_proc  , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&buff2, 2, MPI_DOUBLE, rand_proc+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      rand_ch=rand_d(1.0);
      if(rand_ch < exp(-(buff1[0]-buff2[0])*(1.0/buff1[1]-1.0/buff2[1]))){
	//the exchange happens
	/* MPI_Send(buff1, 2, MPI_DOUBLE, rand_proc+1, 2, MPI_COMM_WORLD); */
	/* MPI_Send(buff2, 2, MPI_DOUBLE, rand_proc  , 3, MPI_COMM_WORLD); */
	ttemp=TEMP[rand_proc];
	TEMP[rand_proc]=TEMP[rand_proc+1];
	TEMP[rand_proc+1]=ttemp;
      }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(TEMP,mpi_count, MPI_DOUBLE, 0, MPI_COMM_WORLD); //THE ARRAY TEMP WAS CHANGED IF THE MOVE WAS ACCEPTED
    MPI_Barrier(MPI_COMM_WORLD);
    PT_set_temp(TEMP[mpi_id]);
    PT_write_params(TEMP, mpi_count, ex);
    }
  }
  printf("Simulation %d ended properly.\n", mpi_id);
  MC_end(mc_n, rx, ry, rz, i, mpi_id);
  MPI_Finalize();
  
  return 0;
}
  
