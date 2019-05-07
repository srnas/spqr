#include "st_impl.h"



void ST_read_parameters(int mpi_id, int mpi_count, int *ST_ind, double *ST_DELTA, double *ST_C, int *mc_ini,  double **ST_g, double **ST_Ts, int *ST_N ){
  char filename[256];
  sprintf(filename, "Sim_temp.pms");
  FILE *ST_PARAMS;
  double tmin, tmax, dt, DELTA, temp, C;
  int NTS, n, itemp, ini, ran;
  if((ST_PARAMS=fopen(filename, "r"))==NULL){
    printf("SIMULATED TEMPERING: Sim_temp.pms not found. Go create one.\n");
    exit(ERR_INPUT);
  }
  else{
    if(mpi_id==0)
      printf("SIMULATED TEMPERING: Sim_temp.pms file found. Taking the parameters from there for all processors.\n");
    fscanf(ST_PARAMS, "%lf %lf %lf %lf %lf %d", &tmin, &tmax, &dt, &DELTA, &C, &ini);
    NTS=(int)((tmax-tmin)/dt)+1;
    *ST_N=NTS;
    *ST_DELTA=DELTA;
    *ST_C=C;
    *mc_ini=ini;
    *ST_g=(double *)malloc(sizeof(double)*NTS);
    *ST_Ts=(double *)malloc(sizeof(double)*NTS);
    
    //for(n=0;n<mpi_count;n++){
    //fscanf(ST_PARAMS, "%d", &itemp);
    //if(mpi_id==n)
    //	*ST_ind=itemp;
    //}
    
    fscanf(ST_PARAMS, "%d", &itemp);
    //if(itemp<0) *ST_ind=mpi_id;
    if(itemp<0) {ran=rand_i(NTS); *ST_ind=ran;}//printf("mpiid = %d  ind = %d\n", mpi_id, ran); }
    else *ST_ind=itemp;
    
    for(n=0;n<NTS;n++){
      (*ST_Ts)[n]=tmin+n*dt;
      fscanf(ST_PARAMS, "%lf", &temp);
      (*ST_g)[n]=temp;
    }
    
  }
  fclose(ST_PARAMS);

  //printf("proc %d has index %d and temperature %lf. Tmin=%lf, Tmax=%lf, dT=%lf\n, NPROCS=%d, NTEMPS=%d and the corresponding g is %lf\t With ini %d , D=%lf and C=%lf\n", mpi_id, *ST_ind,(*ST_Ts)[*ST_ind], tmin, tmax, dt, mpi_count, NTS, (*ST_g)[*ST_ind], ini, DELTA, C);

}

void ST_write_parameters(int mpi_id, int mpi_count, int ST_ind, double ST_DELTA, double ST_C, int mc_ini, double *ST_g, double *ST_Ts, int ST_N ){
  char filename[256];
  sprintf(filename, "new_Sim_temp.p%02d.pms", mpi_id);
  FILE *ST_PARAMS;
  double tmin, tmax, dt, temp;
  int n, itemp;
  
  if((ST_PARAMS=fopen(filename, "w"))==NULL){
    printf("SIMULATED TEMPERING: Can't write parameters!.\n");
    exit(ERR_INPUT);
  }
  else{
    if(mpi_id==0)
      printf("SIMULATED TEMPERING: writing final parameters.\n");
    tmin=ST_Ts[0];
    tmax=ST_Ts[ST_N-1];
    dt=(tmax-tmin)/((double)(ST_N-1));
    fprintf(ST_PARAMS, "%lf %lf %lf %lf %lf %d\n%d\n", tmin, tmax, dt , ST_DELTA, ST_C, mc_ini, ST_ind);
    
    
    /* for(n=0;n<mpi_count;n++){ */
    /*   fprintf(ST_PARAMS, "%d", &itemp); */
    /*   if(mpi_id==n) */
    /* 	*ST_ind=itemp; */
    /* } */
    
    
    for(n=0;n<ST_N;n++){
      fprintf(ST_PARAMS, "%lf ", ST_g[n]);
      
    }
    fprintf(ST_PARAMS, "\n");
  }
  fclose(ST_PARAMS);
    
  //printf("proc %d has index %d and temperature %lf. Tmin=%lf, Tmax=%lf, dT=%lf\n, NPROCS=%d, NTEMPS=%d and the corresponding g is %lf\n", mpi_id, *ST_ind,(*ST_Ts)[*ST_ind], tmin, tmax, dt, mpi_count, NTS, (*ST_g)[*ST_ind]);

}

void ST_update_temperature(double temp){
  mc_target_temp=temp;
}
