#include "sa_impl.h"

void SA_read_params(double *tmax, double *tmin, double *tfac, int *sa_step,  int *NT, double *prev_energ, double *sfac, int *resc_times, int mpi_id){
  char filename[256];
  sprintf(filename, "sim_ann.p%02d.pms", mpi_id);
  
  FILE *SA_PARAMS;
  if((SA_PARAMS=fopen(filename, "r"))==NULL){
    printf("SIMULATED ANNEALING: sim_ann.p%02d.pms file not found!\n", mpi_id);
    printf("SIMULATED ANNEALING: Looking for generic sim_ann.pms file.\n");
    sprintf(filename, "sim_ann.pms");
    if((SA_PARAMS=fopen(filename, "r"))==NULL){
      printf("SIMULATED ANNEALING: sim_ann.pms not found either. Go create one.\n");
      exit(ERR_INPUT);
    }
    else printf("SIMULATED ANNEALING: sim_ann.pms file found. Taking the parameters from there for all processors.\n");
  }
  else printf("SIMULATED ANNEALING: Process %d reading parameters from file sim_ann.%02d.pms\n", mpi_id, mpi_id);
  double temp1, temp2, temp3, temp5, temp6;
  int temp4, temp7, temp8;
  fscanf(SA_PARAMS,"%lf %lf %lf %d %d %lf %lf %d", &temp1, &temp2, &temp3, &temp8, &temp4, &temp5, &temp6, &temp7);
  *tmax=temp1;
  *tmin=temp2;
  *tfac=temp3;
  *sa_step=temp8;
  *NT=temp4;
  *prev_energ=temp5;
  *sfac=temp6;
  *resc_times=temp7;
  fclose(SA_PARAMS);
}

void SA_set_temp(double temp){
  mc_target_temp=temp;
}

void SA_set_mc_trials(double smc_nt_xyz, double smc_ph_xyz, double smc_nt_ang){
  MC_NT_XYZ=smc_nt_xyz;
  MC_PH_XYZ=smc_ph_xyz;
  MC_NT_ANGLE=smc_nt_ang;
}

void SA_init_mc_trials(double *smc_nt_xyz, double *smc_ph_xyz, double *smc_nt_ang, double a, double b, double c){
  *smc_nt_xyz=a;
  *smc_ph_xyz=b;
  *smc_nt_ang=c;
}

void SA_rescale_mc_trials(double *smc_nt_xyz, double *smc_ph_xyz, double *smc_nt_ang, double sfac){
  (*smc_nt_xyz)*=sfac;
  (*smc_ph_xyz)*=sfac;
  (*smc_nt_ang)*=sfac;
}

void SA_save_params(double temp, double tmin, double tfac, int sa_step, int nt, double prevenerg, double sfac, int resc_times, int mpi_id){
  FILE *NEW_SA_PARAMS;
  char filename[256];
  sprintf(filename, "new_sim_ann.p%02d.pms", mpi_id);
  if((NEW_SA_PARAMS=fopen(filename, "w"))==NULL){
    printf("SIMULATED ANNEALING: new_sim_ann.p%02d.pms file could not be created!\n", mpi_id);
  }
  fprintf(NEW_SA_PARAMS,"%lf %lf %lf %d %d %lf %lf %d\n", temp, tmin, tfac, sa_step, nt, prevenerg, sfac, resc_times);
  fclose(NEW_SA_PARAMS);
}
double sa_pow(double base, int expo){
  int i;
  double res=1.0;
  for(i=0;i<expo;i++)
    res=res*base;
  return res;
}
