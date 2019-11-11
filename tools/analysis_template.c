#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DIM 3
#define N_PARTS_PER_NT 5

int MC_N, NT_N;
FILE *md_configs;
double *rx, *ry, *rz;
int *mc_glyc, *mc_puck, *mc_type;
char confname[9]="confs.mc";

void MC_open_analysis_files(char *);
void MC_read_ith_configuration();
void MC_close_analysis_files();

/* declare here analysis files and functions */


/*********************************************/

int main(int nargs, char **args){
  int i,j,d;
  double temperature, energy_t;
  if(nargs != 4){
    fprintf(stderr, "Analysis requires the number of configurations saved in file %s.\n", confname);
    exit(1);
  }
  
  int nconfs=atoi(args[2]);
  int confini=atoi(args[3]);
  MC_open_analysis_files(args[1]);
  
  for(i=0;i<confini;i++)
    MC_read_ith_configuration(&energy_t, &temperature);
  for(i=confini;i<nconfs;i++){
    MC_read_ith_configuration(&energy_t, &temperature);
    /* do analysis here */
      
      
      
  }
  MC_close_analysis_files();
  return 0;
}

void MC_open_analysis_files(char * name){
  int i;
  double tempf;
  if((md_configs=fopen(name, "rb"))==NULL){
    fprintf(stderr,"File %s not found!\n", name);
    exit(1);
  }
  
  fread(&MC_N, sizeof(int),1,md_configs);
  NT_N=MC_N/N_PARTS_PER_NT;
  
  rx=(double*)malloc(sizeof(double *)*MC_N);
  ry=(double*)malloc(sizeof(double *)*MC_N);
  rz=(double*)malloc(sizeof(double *)*MC_N);
  mc_type=(int*)malloc(sizeof(int)*MC_N);
  mc_glyc=(int *)malloc(sizeof(int)*NT_N);
  mc_puck=(int *)malloc(sizeof(int)*NT_N);
  
  
}

void MC_read_ith_configuration(double *energy_t, double *temperature){
  int i,d, nt, glp;
  double tempp[DIM], tempv[DIM];
  fread(energy_t,sizeof(double),1,md_configs);
  fread(temperature,sizeof(double),1,md_configs);
  
  for(i=0;i<MC_N;i++){
    fread(&(mc_type[i]), sizeof(int),1,md_configs);
  }
  
  for(i=0;i<MC_N;i++){
    fread(tempp,sizeof(double),DIM,md_configs);
    rx[i]=tempp[0];ry[i]=tempp[1];rz[i]=tempp[2];
  }
  for(nt=0;nt<NT_N;nt++){
    fread(&glp,sizeof(int),1,md_configs);
    mc_glyc[nt]=glp;
    fread(&glp,sizeof(int),1,md_configs);
    mc_puck[nt]=glp;
  }
}

void MC_close_analysis_files(){
  free(rx);
  free(ry);
  free(rz);
  free(mc_glyc);
  free(mc_puck);
  free(mc_type);
  fclose(md_configs);
}
