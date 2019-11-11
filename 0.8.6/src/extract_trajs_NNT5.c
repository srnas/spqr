#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DIM 3

FILE *md_configs;
double **pos;
int *type;
int MD_N;
int current_conf;
char confname[9]="confs.mc";

void MD_open_analysis_files(char *);
void MD_read_ith_configuration();
void MD_close_analysis_files();

/* declare here analysis files and functions */


/*********************************************/

int main(int nargs, char **args){
  int i,j,d;
  if(nargs != 4){
    fprintf(stderr, "Analysis requires the number of configurations saved in file %s.\n", confname);
    exit(1);
  }
  
  int nconfs=atoi(args[2]);
  int confini=atoi(args[3]);
  char ttype;
  char name[256];
  MD_open_analysis_files(args[1]);
  sprintf(name, "confs.xyz");
  FILE *tempconf;
  tempconf=fopen(name, "w");

  for(i=0;i<nconfs;i++){
    MD_read_ith_configuration();
    if(i>=confini){
      /* do analysis here */
      fprintf(tempconf,"%d\nconf %d\n",(int)MD_N, i);
      for(j=0;j<MD_N;j++){
	if(type[j]==-4)
	  ttype='P';
	if(type[j]==-3)
	  ttype='S';
	if(type[j]==-2)
	  ttype='Y';
	if(type[j]==-1)
	  ttype='X';
	if(type[j]==0)
	ttype='A';
	if(type[j]==1)
	  ttype='U';
	if(type[j]==2)
	  ttype='G';
	if(type[j]==3)
	  ttype='C';

	
	fprintf(tempconf, "%c ",ttype);
	
	for(d=0;d<DIM;d++)
	  fprintf(tempconf, "%f ", pos[j][d]);
	fprintf(tempconf, "\n");
      }
    }
  }
  MD_close_analysis_files();
  fclose(tempconf);

  return 0;
}

void MD_open_analysis_files(char * name){
  int i;
  double temp;
  current_conf=0;
  if((md_configs=fopen(name, "rb"))==NULL){
    fprintf(stderr,"File %s not found!\n", name);
    exit(1);
  }
  
  fread(&MD_N, sizeof(int),1,md_configs);
  
  pos=(double**)malloc(sizeof(double *)*MD_N);
  type=(int*)malloc(sizeof(int)*MD_N);
  for(i=0;i<MD_N;i++)
    pos[i]=(double *)malloc(sizeof(double)*DIM);
  
  for(i=0;i<MD_N;i++){
    fread(&type[i], sizeof(int),1,md_configs);
  }
  
  
}

void MD_read_ith_configuration(){
  int i,d;
  double tempp[DIM], tempv[DIM];
  for(i=0;i<MD_N;i++){
    fread(tempp,sizeof(double),DIM,md_configs);
    for(d=0;d<DIM;d++){
      pos[i][d]=(double)tempp[d];
    }
  }
}

void MD_close_analysis_files(){
  fclose(md_configs);
}
