#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

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
#define MAXCHAR 256
#define FLIP_GLYC 0
#define FLIP_PUCK 1
#define ERR_INPUT 1
#define ERR_OUTPUT 2

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
void write_checkpoint(int , char* ,  double *, double *, double *, int, double, double, int *, int *, double *);
/* declare here analysis files and functions */
//void MC_write_pdb(FILE * , int, double *, double *, double *, double , double , int *, int *, int *);

/*********************************************/
static void show_usage(){
    fprintf(stderr, "Usage: MCEDIT -i <input_filename> -o <output_filename> [ -s <fasta> ] [ -L <list> ] [ -F <list> ] [ -p <offset> ]\n");
    fprintf(stderr, "-s allows change of glycosidic state and sugar pucker according to secondary structure.\n");
    //fprintf(stderr, "-G allows change of glycosidic state and sugar pucker according to list, while fixing the state of the rest.\n");
    fprintf(stderr, "-L crops the structure to the listed nucleotides.\n");
    fprintf(stderr, "-p is the index of the first nucleotide, for consistency with other formats.\n");
    fprintf(stderr, "-F fixes the whole molecule except for the nucleotides in the list.\n");
    fprintf(stderr, "List of nucleotides must be written as X-Y,Z, without spaces in between.\n");
    exit(ERR_INPUT);
}

int main(int argc, char **argv){
  int i,j,d;
  printf("MCEDIT : edit SPQR mc files.\n");
  //int nconfs=atoi(args[2]);
  //int confini=atoi(args[3]);
  char ttype;
  double temperature, energy_t;
  int ntfirst=-1,ntlast=-1;
  char *lline=NULL, *pdbrectyp, *tmp, *tmp2, *stmp, *cpline;
  double dtemp1, dtemp2;
  int gr, l, ll, at=0, itemp, ERMSD_N_GROUPS;
  char s1[MAXCHAR], s2[MAXCHAR], s3[MAXCHAR];
  double kermsd;
  static size_t st_l=0;
  int group, *nnt_group, ntg, nt, NATS_SS;//, grflag;
  int **ntind_group;
  int opt;
  char mcfilename[MAXCHAR], outfilename[MAXCHAR], ssfilename[MAXCHAR], ntrawlist[MAXCHAR], ntfixlist[MAXCHAR], ntglplist[MAXCHAR],ntfrzlist[MAXCHAR];
  sprintf(mcfilename," ");sprintf(outfilename,"out.mc");sprintf(ssfilename," ");sprintf(ntrawlist," ");sprintf(ntglplist, " ");sprintf(ntfrzlist," ");
  int verboseflag=0;
  int clistflag=0, ssflag=0, frzflag=0,glpflag=0;
  int *ntlist, *chainlist,*unfrzlist;
  int chainshift;
  int offset=0;
  while ((opt = getopt(argc, argv, "i:o:s:L:p:F:v")) != -1) {
    switch (opt)
      {
      case 'i':
        sprintf(mcfilename, "%s", optarg);
        break;
      case 'o':
        sprintf(outfilename, "%s", optarg);
        break;
      case 's':
        sprintf(ssfilename, "%s", optarg);
	ssflag=1;
	break;
      case 'L':
	sprintf(ntrawlist, "%s", optarg);
	clistflag=1;
	break;
      case 'F':
        sprintf(ntfrzlist, "%s", optarg);
	frzflag=1;
	break;
      case 'p':
	offset=atoi(optarg);
	break;
      case 'v':
        verboseflag=1;
        break;
      default:
	show_usage();
      }
  }
  open_analysis_files(mcfilename);
  MC_read_ith_configuration(&energy_t, &temperature);
  ///////////////DO HERE WHAT WE WANT TO DO WITH THE CONFIGURATION////////////////
  int *secstr;secstr=(int *)malloc(sizeof(int)*NT_N);for(i=0;i<NT_N;i++){secstr[i]=0;}
  int chini, chend;
  if(ntfirst<0) ntfirst=0; if(ntlast<0) ntlast=NT_N;
  ntlist=(int *)malloc(NT_N*sizeof(int));
  chainlist=(int *)malloc(NT_N*sizeof(int));
  unfrzlist=(int *)malloc(NT_N*sizeof(int));
  
  chainshift=0;
  for (i=0;i<NT_N;i++) chainlist[i]=0;
  for (i=0;i<NT_N;i++) ntlist[i]=1;
  for (i=0;i<NT_N;i++) unfrzlist[i]=1;
  
  if(frzflag==1){
    for (i=0;i<NT_N;i++) unfrzlist[i]=0;
    char *tok, *toktok, *temptok;
    int nt1, nt2, twoflag=0;
    tok=strtok(ntfrzlist,",");
    while (tok != NULL){
      //here we divide by hand the string
      chend=0;chini=0;
      twoflag=0;
      for(chend=0;chend<strlen(tok);chend++) {
	if(tok[chend]=='-'){
	  twoflag=1;
	  temptok=realloc(temptok, (chend-chini)*sizeof(char));
	  temptok[chend-chini]='\0';
	  strncpy(temptok,tok+chini,chend-chini);
	  nt1=atoi(temptok);
	  nt2=nt1;
	  chini=chend+1;
	}
      }
      //printf("realloc with %d   %d, %d\n", chend-chini, chend, chini);
      temptok=realloc(temptok, (chend-chini)*sizeof(char));
      temptok[chend-chini]='\0';
      strncpy(temptok,tok+chini,chend-chini);
      nt2=atoi(temptok);
      if(twoflag==0) nt1=nt2;
      //here apply the offset
      nt1-=offset;
      nt2-=offset;
      //printf("%d %d\n", nt1, nt2);
      if(nt1<0 || nt1 > NT_N || nt2<0 || nt2 > NT_N ){printf("ERROR: invalid nucleotide index %d or %d.\n", nt1, nt2); exit(ERR_INPUT);}
      tok=strtok(NULL,",");
      for(i=nt1;i<nt2+1;i++){
	unfrzlist[i]=1;
	chainlist[i]=chainshift;
      }
      chainshift++;
    }
    //here we freeze them
    for(i=0;i<NT_N;i++){
      if(unfrzlist[i]==0){
	is_mob[i]=FR_MOB_FROZ;
	//printf("%d is frozen\n", i);
      }
      else{
	is_mob[i]=FR_MOB_FULL;
	//printf("%d is mobile\n", i);
	is_flipp[i]=GLP_BOTH;
      }
    }
  }
  

  if(clistflag==1){
    for (i=0;i<NT_N;i++) ntlist[i]=0;
    char *tok, *toktok, *temptok;
    int nt1, nt2, twoflag=0;
    tok=strtok(ntrawlist,",");
    while (tok != NULL){
      //here we divide by hand the string
      chend=0;chini=0;
      twoflag=0;
      //printf("%s  %d\n", tok, strlen(tok));
      
      for(chend=0;chend<strlen(tok);chend++) {
	
	//printf("%c\n", tok[chend]);
	//printf("%d   %c\n", chend, tok[chend]);
	if(tok[chend]=='-'){
	  twoflag=1;
	  temptok=realloc(temptok, (chend-chini)*sizeof(char));
	  temptok[chend-chini]='\0';
	  strncpy(temptok,tok+chini,chend-chini);
	  //temptok[chend-chini]=NULL;
	  //printf("copied %s\n", temptok);
	  nt1=atoi(temptok);
	  nt2=nt1;
	  chini=chend+1;
	  //printf("part %d  %d\n", nt1,nt2);
	}
      }
      temptok=realloc(temptok, (chend-chini)*sizeof(char));
      temptok[chend-chini]='\0';
      strncpy(temptok,tok+chini,chend-chini);
      //temptok[chend-chini]=NULL;
      
      nt2=atoi(temptok);
      if(twoflag==0) nt1=nt2;
      //printf("%d  %d\n", nt1, nt2);
      //here apply the offset
      nt1-=offset;
      nt2-=offset;
      if(nt1<0 || nt1 > NT_N || nt2<0 || nt2 > NT_N ){printf("ERROR: invalid nucleotide index.\n"); exit(ERR_INPUT);}
      tok=strtok(NULL,",");
      for(i=nt1;i<nt2+1;i++){
	ntlist[i]=1;
	chainlist[i]=chainshift;
      }
      chainshift++;
    }
  }

  
  //for (i=0;i<NT_N;i++) printf("%d ", chainlist[i]);
  //for(i=0;i<NT_N;i++) printf("%d ", ntlist[i]); printf("\n");
  if((FSECSTR=fopen(ssfilename, "r"))==NULL){
    if(verboseflag)
      printf("No secondary structure in file fasta file found.\n");
  }
  else {
    printf("Secondary structure taken from file %s\n", ssfilename);
    l=getline(&lline, &st_l, FSECSTR);
    while(lline[0]=='#')    l=getline(&lline, &st_l, FSECSTR);
    l=getline(&lline, &st_l, FSECSTR);
    l=getline(&lline, &st_l, FSECSTR);
    int cnt=0;
    
    while(lline[cnt]!='\n'){
      printf("%c", lline[cnt]);
      if(lline[cnt]=='.'){
	if(verboseflag)
	  printf("%d is flippable\n", cnt);
	is_flipp[cnt]=GLP_BOTH;

      }
      else 
	is_flipp[cnt]=GLP_FIXED;
      cnt++;
    }
    /*secstr[ntind_group[gr][i]]=1;
      for(i=0;i<NT_N;i++){
      //is_flipp[i]=GLP_FIXED;
      // is_mob[i]=FR_MOB_FULL;
      if(secstr[i]!=1){
      printf("%d is flippable\n", i);
      is_flipp[i]=GLP_BOTH;
      }
      else 
      is_flipp[i]=GLP_FIXED;
      }
    */
  }
  
  ////////////////////////////////////////////////////////////////////////////////
  current_conf=0;
  write_checkpoint(MC_N, outfilename, posx, posy, posz, current_conf, energy_t, temperature, ntlist, chainlist, add_data); //this writes a binary checkpoint
  close_analysis_files();
  
  return 0;
}

void write_checkpoint(int mc_n, char* confname,  double *rx, double *ry, double *rz, int iter, double energy_t, double temperature, int *nlist, int *chlist, double *add_data){
  int i, tempi, nt;
  double tempf;
  int nt_n=mc_n/N_PARTS_PER_NT;
  FILE *chkpnt;
  char word[3]="END";
  int tchain=0;
  //int atini=nti*N_PARTS_PER_NT;
  //int atlas=ntl*N_PARTS_PER_NT;
  int effmc_n=0;
  for(i=0;i<nt_n;i++) if(nlist[i]==1) effmc_n+=N_PARTS_PER_NT;
  if((chkpnt=fopen(confname,"wb"))!=NULL){
    fwrite(&effmc_n, sizeof(int), 1, chkpnt);
    for(i=0;i<mc_n;i++)
      //for(i=atini;i<atlas;i++)
      if(nlist[(int)(i/N_PARTS_PER_NT)]==1)
	fwrite(&(type[i]), sizeof(int), 1, chkpnt);
    tempf=(double)(energy_t);
    fwrite(&tempf,sizeof(double),1,chkpnt);
    fwrite(&temperature,sizeof(double),1,chkpnt);
    //and then the configuration
    for(i=0;i<mc_n;i++){
      if(nlist[(int)(i/N_PARTS_PER_NT)]==1){
	tempf=rx[i];
	fwrite(&tempf,sizeof(double),1,chkpnt);
	tempf=ry[i];
	fwrite(&tempf,sizeof(double),1,chkpnt);
	tempf=rz[i];
	fwrite(&tempf,sizeof(double),1,chkpnt);
      }
    }
    int rchain;
    for(nt=0;nt<nt_n;nt++){
      if(nlist[nt]==1){
	rchain=0;
	//chains[nt]+chlist[nt];
	//if(nt>0) {if(chains[nt]!=chains[nt-1]) tchain++;}
	tempi=(int)(glyc[nt]);
	fwrite(&tempi,sizeof(int),1,chkpnt);
	tempi=(int)(puck[nt]);
	fwrite(&tempi,sizeof(int),1,chkpnt);
	tempi=(int)(is_flipp[nt]);
	fwrite(&tempi,sizeof(int),1,chkpnt);
	tempi=(int)(is_mob[nt]);
	fwrite(&tempi,sizeof(int),1,chkpnt);
	fwrite(&(rchain),sizeof(int),1,chkpnt);
      }
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
    exit(ERR_OUTPUT);
  }
}

void open_analysis_files(char * name){
  int i;
  double temp;
  //double postemp[DIM];
  current_conf=0;
  if((mc_configs=fopen(name, "rb"))==NULL){
    fprintf(stderr,"Invalid input file!\n");
    show_usage();
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
