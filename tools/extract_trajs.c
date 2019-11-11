#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DIM 3
#define N_PARTS_PER_NT 5
#define IPHO 4
#define ISUG 3
#define IBAS 0
#define IX 1
#define IY 2
#define TYP_ADENINE  0
#define TYP_URACIL   1
#define TYP_GUANINE  2
#define TYP_CYTOSINE 3
#define GLYC_A 0
#define GLYC_H 1
#define GLYC_S 2
#define PUCK_3 0 
#define PUCK_2 1


FILE *md_configs;
double *posx, *posy, *posz;
int *glyc, *puck;
int *type;
int MC_N, NT_N;
int current_conf;
char confname[9]="confs.mc";

void MC_open_analysis_files(char *);
void MC_read_ith_configuration(double *, double *);
void MC_close_analysis_files();

/* declare here analysis files and functions */
void MC_write_pdb(FILE * , int, double *, double *, double *, double , double , int *, int *, int *);

/*********************************************/

int main(int nargs, char **args){
  int i,j,d;
  if(nargs != 4){
    fprintf(stderr, "Analysis requires the file name, the number of configurations saved there and the initial config number.\n");
    exit(1);
  }
  
  int nconfs=atoi(args[2]);
  int confini=atoi(args[3]);
  char ttype;
  char name[256];
  double temperature, energy_t;
  
  MC_open_analysis_files(args[1]);
  sprintf(name, "confs.pdb");
  FILE *tempconf;
  tempconf=fopen(name, "w");

  for(i=0;i<nconfs;i++){
    MC_read_ith_configuration(&energy_t, &temperature);
    if(i>=confini){
      /* do analysis here */
      fprintf(tempconf, "MODEL %d\n", i);
      MC_write_pdb(tempconf,  MC_N,  posx, posy, posz,  energy_t, temperature, glyc, puck, type);
      

      fprintf(tempconf, "ENDMDL\n");

    }
  }
  MC_close_analysis_files();
  fclose(tempconf);

  return 0;
}


void MC_write_pdb(FILE * pdbfile, int mc_n, double *rx, double *ry, double *rz, double energy_t, double temperature, int *mc_glyc, int *mc_puck, int *mc_types){
  int at, nt_n=mc_n/N_PARTS_PER_NT, nt, tchain=0, nintdig;
  //FILE *pdbfile;
  char at_num[10], at_name[10], res_name[10], res_num[10], chain_id[10], xp[256], yp[256], zp[256], gly[10], pck[10], glp[10], frz[10];
  char pdbname[256];
  fprintf(pdbfile, "REMARK ENERGY %lf  TEMPERATURE %lf\n", energy_t, temperature);
  for(at=0;at<mc_n;at++){
    nt=(int)(at/N_PARTS_PER_NT);
    sprintf(at_num, "%5d", at);
    if(mc_types[at]==-1) sprintf(at_name,"XVEC");
    else if(mc_types[at]==-2) sprintf(at_name,"YVEC");
    else if(mc_types[at]==-3) sprintf(at_name,"SUGR");
    else if(mc_types[at]==-4) sprintf(at_name,"PHOS");
    else sprintf(at_name,"BASE");
    
    if(mc_types[nt*N_PARTS_PER_NT]==TYP_ADENINE)  sprintf(res_name, "  A");
    else if(mc_types[nt*N_PARTS_PER_NT]==TYP_URACIL)   sprintf(res_name, "  U");
    else if(mc_types[nt*N_PARTS_PER_NT]==TYP_GUANINE)  sprintf(res_name, "  G");
    else if(mc_types[nt*N_PARTS_PER_NT]==TYP_CYTOSINE) sprintf(res_name, "  C");
    else {printf("Residue type not recognized writing pdb\n"); exit(1);}
    //if(nt>0 && at%N_PARTS_PER_NT==0) {if(MC_are_neighbors(nt, nt-1)==0) tchain++;}
    tchain=0;
    sprintf(chain_id, "%d", tchain);
    sprintf(res_num, "%4d", nt);
      
      sprintf(xp, "%d", (int)rx[at]);
      nintdig=strlen(xp); if(nintdig>8){printf("ERROR: particle %d out of range for x coordinate in pdb format\n", at); exit(1);}if(rx[at]<0 && (int)rx[at]==0) nintdig++;
      sprintf(xp, "%.*f ", rx[at], 6-nintdig);
      sprintf(yp, "%d", (int)ry[at]);
      nintdig=strlen(yp); if(nintdig>8){printf("ERROR: particle %d out of range for y coordinate in pdb format\n", at); exit(1);}if(ry[at]<0 && (int)ry[at]==0) nintdig++;
      sprintf(yp, "%.*f ", ry[at], 6-nintdig);
      sprintf(zp, "%d", (int)rz[at]);
      nintdig=strlen(zp); if(nintdig>8){printf("ERROR: particle %d out of range for z coordinate in pdb format\n", at); exit(1);}if(rz[at]<0 && (int)rz[at]==0) nintdig++;
      sprintf(zp, "%.*f ", rz[at], 6-nintdig);
      /* sprintf(xp, "%8.3f", rx[at]); */
      /* sprintf(yp, "%8.3f", ry[at]); */
      /* sprintf(zp, "%8.3f", rz[at]); */
      
      fprintf(pdbfile, "ATOM  %s %s %s %s%s    %s%s%s", at_num, at_name, res_name, chain_id, res_num, xp, yp, zp);
      //printf( "ATOM  %s     %d  %d   %d", res_name,  at, nt, mc_types[nt*N_PARTS_PER_NT]);
      if(at%N_PARTS_PER_NT==0){
      	if(mc_glyc[nt]==GLYC_A) sprintf(gly,"A");
      	if(mc_glyc[nt]==GLYC_H) sprintf(gly,"H");
      	if(mc_glyc[nt]==GLYC_S) sprintf(gly,"S");
      	if(mc_puck[nt]==PUCK_3) sprintf(pck,"3");
      	if(mc_puck[nt]==PUCK_2) sprintf(pck,"2");
      	//if(glp_is_flippable[nt]==GLP_FIXED) sprintf(glp,"N");
      	//if(glp_is_flippable[nt]==GLP_GLYC) sprintf(glp,"G");
      	//if(glp_is_flippable[nt]==GLP_PUCK) sprintf(glp,"P");
      	//if(glp_is_flippable[nt]==GLP_BOTH) sprintf(glp,"A");
      	//if(fr_is_mobile[nt]==FR_MOB_FROZ) sprintf(frz,"N");
      	//if(fr_is_mobile[nt]==FR_MOB_BASE) sprintf(frz,"B");
      	//if(fr_is_mobile[nt]==FR_MOB_PHOS) sprintf(frz,"P");
      	//if(fr_is_mobile[nt]==FR_MOB_FULL) sprintf(frz,"A");
      	//fprintf(pdbfile, "  %s%s%s%s", gly, pck, glp, frz);
	fprintf(pdbfile, "  %s%sNN", gly, pck);
      }
      fprintf(pdbfile, "\n");
    }
}



void MC_open_analysis_files(char * name){
  int i;
  double temp;
  //double postemp[DIM];
  current_conf=0;
  if((md_configs=fopen(name, "rb"))==NULL){
    fprintf(stderr,"File %s not found!\n", name);
    exit(1);
  }
  
  fread(&MC_N, sizeof(int),1,md_configs);
  NT_N=MC_N/N_PARTS_PER_NT;
  posx=(double*)malloc(sizeof(double)*MC_N);
  posy=(double*)malloc(sizeof(double)*MC_N);
  posz=(double*)malloc(sizeof(double)*MC_N);
  type=(int*)malloc(sizeof(int)*MC_N);
  puck=(int*)malloc(sizeof(int)*NT_N);
  glyc=(int*)malloc(sizeof(int)*NT_N);
  //for(i=0;i<MC_N;i++){
  //postemp=(double *)malloc(sizeof(double)*DIM);
  //}
    
  for(i=0;i<MC_N;i++){
    fread(&type[i], sizeof(int),1,md_configs);
  }
}

void MC_read_ith_configuration(double *energy_t, double *temperature){
  int i,d, glp;
  double tempp[DIM], tempv[DIM];
  fread(energy_t,sizeof(double),1,md_configs);
  fread(temperature,sizeof(double),1,md_configs);
  for(i=0;i<MC_N;i++){
    fread(tempp,sizeof(double),DIM,md_configs);
    //pos[i][d]=(double)tempp[d];
    posx[i]=tempp[0];
    posy[i]=tempp[1];
    posz[i]=tempp[2];
  }
  for(i=0;i<NT_N;i++){
    fread(&glp,sizeof(int),1,md_configs);
    glyc[i]=glp;
    fread(&glp,sizeof(int),1,md_configs);
    puck[i]=glp;
  }
}

void MC_close_analysis_files(){
  fclose(md_configs);
}
