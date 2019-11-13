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
#define ANN_NPARAMS 8
#define FR_MOB_FULL 3
#define FR_MOB_PHOS 1
#define FR_MOB_BASE 2
#define FR_MOB_FROZ 0
#define GLP_FIXED  0
#define GLP_GLYC   1
#define GLP_PUCK   2
#define GLP_BOTH   3

#define FLIP_GLYC 0
#define FLIP_PUCK 1

FILE *FSECSTR;

FILE *mc_configs;
double *posx, *posy, *posz;
int *glyc, *puck;
int *type;
double add_data[ANN_NPARAMS];
int *is_flipp;
int *is_mob;
int *chains;
int MC_N, NT_N;
int current_conf;
char word[3];
double pidum;

void open_analysis_files(char *);
void MC_read_ith_configuration(double *, double *);
void close_analysis_files();
void print_help();
void write_checkpoint(int , char* ,  double *, double *, double *, int, double, double, double *);
/* declare here analysis files and functions */
//void MC_write_pdb(FILE * , int, double *, double *, double *, double , double , int *, int *, int *);

/*********************************************/

int main(int nargs, char **args){
  int i,j,d;
  if(nargs != 3){
    print_help();
    exit(1);
  }
  //int nconfs=atoi(args[2]);
  //int confini=atoi(args[3]);
  char ttype;
  char outname[256];
  double temperature, energy_t;
  
  char *lline=NULL, *pdbrectyp, *tmp, *tmp2, *stmp, *cpline;
  double dtemp1, dtemp2;
  int gr, l, ll, at=0, itemp, ERMSD_N_GROUPS;
  char s1[256], s2[256], s3[256];
  static size_t st_l=0;
  int group, *nnt_group, ntg, nt, NATS_SS;//, grflag;
  int **ntind_group;
    
  open_analysis_files(args[1]);
  MC_read_ith_configuration(&energy_t, &temperature);
  sprintf(outname, "out.mc");
  ///////////////DO HERE WHAT WE WANT TO DO WITH THE CONFIGURATION////////////////
  int *secstr;secstr=(int *)malloc(sizeof(int)*NT_N);for(i=0;i<NT_N;i++){secstr[i]=0;}
  
  if(!strcmp(args[2], "-f")){
    if((FSECSTR=fopen("ermsd_frags.lst", "r"))==NULL){
      printf("No secondary structure in file ermsd_frags.lst found.\n");
    }
    else{
      printf("Secondary structure taken from file ermsd_frags.lst\n");
      l=getline(&lline, &st_l, FSECSTR);
      sscanf(lline,"%s %s %s %d %lf %lf", s1, s2, s3,  &itemp, &dtemp1, &dtemp2);
      ERMSD_N_GROUPS=itemp;
      ntind_group=(int **) malloc(sizeof(int *)*ERMSD_N_GROUPS);
      nnt_group=(int *)malloc(sizeof(int)*ERMSD_N_GROUPS);
      NATS_SS=0;
      //read groups
      for(group=0;group<ERMSD_N_GROUPS;group++){
	l=getline(&lline, &st_l, FSECSTR);
	cpline=strndup(lline, 256);
	sscanf(cpline, "%s %s %s", s1, s2, s3);
	nnt_group[group]=0;
	tmp=strtok(cpline, " ");while(tmp!=NULL){if(atoi(tmp)!=0 ||(atoi(tmp)==0 && tmp[0]=='0')) {nnt_group[group]++;} tmp=strtok(NULL, " ");}free(cpline);
	ntind_group[group]=(int *)malloc(sizeof(int)*nnt_group[group]);
	NATS_SS+=nnt_group[group];
	ntg=0;
	tmp=strtok(lline, " ");while(tmp!=NULL){if(atoi(tmp)!=0 ||(atoi(tmp)==0 && tmp[0]=='0')) {ntind_group[group][ntg]=atoi(tmp);ntg++;} tmp=strtok(NULL, " ");};
      }
      for(gr=0;gr<ERMSD_N_GROUPS;gr++)
	for(i=0;i<nnt_group[gr];i++)
	  secstr[ntind_group[gr][i]]=1;
    }
    
    for(i=0;i<NT_N;i++){
      is_flipp[i]=GLP_FIXED;
      is_mob[i]=FR_MOB_FULL;
      if(secstr[i]==1)
	is_mob[i]=FR_MOB_FROZ;
    }
    
  }
  ////////////////////////////////////////////////////////////////////////////////
  current_conf=0;
  write_checkpoint(MC_N, outname, posx, posy, posz, current_conf, energy_t, temperature, add_data); //this writes a binary checkpoint
  close_analysis_files();
  
  return 0;
}

void write_checkpoint(int mc_n, char* confname,  double *rx, double *ry, double *rz, int iter, double energy_t, double temperature, double *add_data){
  int i, tempi, nt;
  double tempf;
  int nt_n=mc_n/N_PARTS_PER_NT;
  FILE *chkpnt;
  char word[3]="END";
  int tchain=0;
  
  if((chkpnt=fopen(confname,"wb"))!=NULL){
    fwrite(&mc_n, sizeof(int), 1, chkpnt);
    for(i=0;i<mc_n;i++)
      fwrite(&(type[i]), sizeof(int), 1, chkpnt);
    tempf=(double)(energy_t);
    fwrite(&tempf,sizeof(double),1,chkpnt);
    fwrite(&temperature,sizeof(double),1,chkpnt);
    //and then the configuration
    for(i=0;i<mc_n;i++){
      tempf=rx[i];
      fwrite(&tempf,sizeof(double),1,chkpnt);
      tempf=ry[i];
      fwrite(&tempf,sizeof(double),1,chkpnt);
      tempf=rz[i];
      fwrite(&tempf,sizeof(double),1,chkpnt);
    }
    for(nt=0;nt<nt_n;nt++){
      if(nt>0) {if(chains[nt]!=chains[nt-1]) tchain++;}
      tempi=(int)(glyc[nt]);
      fwrite(&tempi,sizeof(int),1,chkpnt);
      tempi=(int)(puck[nt]);
      fwrite(&tempi,sizeof(int),1,chkpnt);
      tempi=(int)(is_flipp[nt]);
      fwrite(&tempi,sizeof(int),1,chkpnt);
      tempi=(int)(is_mob[nt]);
      fwrite(&tempi,sizeof(int),1,chkpnt);
      fwrite(&tchain,sizeof(int),1,chkpnt);
    }
    //finally, the random seed
    fwrite(&pidum,sizeof(double),1,chkpnt);
    //and the time step
    fwrite(&iter,sizeof(int),1,chkpnt);
    
/* #ifdef TANN */
/*     //in this case, add_data is the annealing parameters */
/*     char word[3]="ANN"; */
/*     int tint=(int)ANN_NPARAMS; */
/*     fwrite(word, sizeof(char),3,chkpnt); */
/*     fwrite(&tint, sizeof(int),1,chkpnt); */
/*     fwrite(add_data, sizeof(double),ANN_NPARAMS,chkpnt); */
/* #endif */
    fwrite(word, sizeof(char),3,chkpnt);
    fclose(chkpnt);
  }
  else{
    printf("Unable to write checkpoint!\n");
    exit(2);
  }
}

void open_analysis_files(char * name){
  int i;
  double temp;
  //double postemp[DIM];
  current_conf=0;
  if((mc_configs=fopen(name, "rb"))==NULL){
    fprintf(stderr,"File %s not found!\n", name);
    exit(1);
  }
  
  fread(&MC_N, sizeof(int),1,mc_configs);
  NT_N=MC_N/N_PARTS_PER_NT;
  posx=(double*)malloc(sizeof(double)*MC_N);
  posy=(double*)malloc(sizeof(double)*MC_N);
  posz=(double*)malloc(sizeof(double)*MC_N);
  type=(int*)malloc(sizeof(int)*MC_N);
  puck=(int*)malloc(sizeof(int)*NT_N);
  glyc=(int*)malloc(sizeof(int)*NT_N);
  chains=(int*)malloc(sizeof(int)*NT_N);
  is_flipp=(int*)malloc(sizeof(int)*NT_N);
  is_mob=(int*)malloc(sizeof(int)*NT_N);
  for(i=0;i<MC_N;i++){
    fread(&type[i], sizeof(int),1,mc_configs);
  }
}

void MC_read_ith_configuration(double *energy_t, double *temperature){
  int i,d, att;
  double tempp[DIM], tempv[DIM],  datt;
  fread(energy_t,sizeof(double),1,mc_configs);
  fread(temperature,sizeof(double),1,mc_configs);
  for(i=0;i<MC_N;i++){
    fread(tempp,sizeof(double),DIM,mc_configs);
    //pos[i][d]=(double)tempp[d];
    posx[i]=tempp[0];
    posy[i]=tempp[1];
    posz[i]=tempp[2];
  }
  for(i=0;i<NT_N;i++){
    fread(&att,sizeof(int),1,mc_configs);
    glyc[i]=att;
    fread(&att,sizeof(int),1,mc_configs);
    puck[i]=att;
    fread(&att,sizeof(int),1,mc_configs);
    is_flipp[i]=att;
    fread(&att,sizeof(int),1,mc_configs);
    is_mob[i]=att;
    fread(&att,sizeof(int),1,mc_configs);
    chains[i]=att;
  }
  fread(&datt,sizeof(double),1,mc_configs); //*idum=att
  pidum=datt;
  fread(&att,sizeof(int),1,mc_configs); //tempiter
  current_conf=att;
  //in this case, add_data is the annealing parameters
  if(fread(word, sizeof(char),3,mc_configs)==3){
    fread(&att, sizeof(int),1,mc_configs);
    //add_data must already be allocated!
    fread(add_data, sizeof(double),ANN_NPARAMS,mc_configs);  
  }
}

void close_analysis_files(){
  fclose(mc_configs);
}

void print_help(){
  printf("SPQR: Modification of mc files\n");
  printf("Usage: MCMOD <input_config> <arg>\n");
  printf("<input_config> must be a single binary .mc file in the spqr format.\n");
  printf("Values of <arg> can be specified with \"-c\", where c is a character:\n");
  printf("(f)  : Freeze all the nucleotides involved in the secondary structure; releasing the positions of the rest.\n");
  printf("(F)  : Freeze all the nucleotides involved in the secondary structure; releasing the positions, puckers and glycosidic bond angles of the rest.\n");
  printf("(p)  : Dump to spqr .pdb file.\n");
  printf("REMEMBER THAT IT DELETES THE INFORMATION OF PREVIOUS ANNEALINGS AND RESETS THE TIME STEP!!\n");
}
    
