#include "pt_impl.h"

/* void SA_read_params(double *tmax, double *tmin, double *tfac, int *NT){ */
/*   FILE *SA_PARAMS; */
/*   if((SA_PARAMS=fopen("sim_ann.pms", "r"))==NULL){ */
/*     printf("SIMULATED ANNEALING: sim_ann.pms file not found!\n"); */
/*   } */
/*   double temp1, temp2, temp3; */
/*   int temp4; */
/*   fscanf(SA_PARAMS,"%lf %lf %lf %d", &temp1, &temp2, &temp3, &temp4); */
/*   *tmax=temp1; */
/*   *tmin=temp2; */
/*   *tfac=temp3; */
/*   *NT=temp4; */
/*   fclose(SA_PARAMS); */
/* } */

void PT_set_temp(double temp){
  mc_target_temp=temp;
}

/* double PT_init_params(int id, int tot){ */
/*   double Tmax=1.0; */
/*   double ans=Tmax*((double)(id))/((double)(tot)); */
/*   return ans; */
/* } */

void PT_create_params(double *temp, int nproc){
  int i;
  double ttemp;
  temp=(double *)malloc(nproc*sizeof(double));
  FILE *PT_PARAMS;
  if((PT_PARAMS=fopen("parall_temp.pms", "r"))==NULL){
    printf("No parallel tempering (parall_temp.pms) file found!\nGo create one.\n");
    exit(1);}
  for(i=0;i<nproc;i++){
    fscanf(PT_PARAMS,"%lf",&ttemp);
    temp[i]=ttemp;
  }
  
  fclose(PT_PARAMS);
}

void PT_write_params(double *temp, int nproc, int step){
  int i;
  FILE *NPT_PARAMS;
  char str[256];
  sprintf(str, "parall_temp.03%d.pms", step);
  NPT_PARAMS=fopen(str,"w");
  for(i=0;i<nproc;i++){
    fprintf(NPT_PARAMS, "%lf ",temp[i]);
  }
  fprintf(NPT_PARAMS, "\n");
  fclose(NPT_PARAMS);
}
