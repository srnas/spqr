#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mc.h"
#include <unistd.h>

#define LX 30
#define LY 30
#define LZ 30
#define NSOLUTE 868
#define MC_ITER 20000
#define EVERYMILLION 500000
#define BAD -10
#define MAXCHAR 256

FILE *md_configs;
char confname[MAXCHAR];

char checkpoint[MAXCHAR];
/* MC variables */
int mc_n;
int mc_iter;
/****************/

void read_ernwin_loops(char*, int*, int**, int***, char***);

void show_usage(){
  printf("Usage: SPQR_ENERG -i <input_config> <options>\n");
  printf("<input_config> can be a .pdb file in the spqr format or a .mc binary configuration.\n");
  printf("Option values:\n");
  printf("-s  : Secondary structure\n");
  printf("-a  : Full annotations\n");
  printf("-t  : Total energy (detailed)\n");
  printf("-f  : Total energy\n");
  printf("-r  : Radius of gyration\n");
  printf("-e <ernwin_file> -T <trj_file> : Ernwin coord file containing loop definitions. Trajectory in mc format.\n");
  printf("This requires the file params.spqr to be in the current directory.\n");

}

int main(int argc, char **argv) {
  printf("SPQR annotation/energy calculation.\n");
  double *rx,*ry,*rz;
  int **LOOPINDXS;
  int NLOOPS;
  int *LOOPNNT;
  char **LOOPLABELS;
  
  int i,j,tempinit;
  int base, face, stf;
  double TOT, BB1, BB2, ST_and_BB2, WCE_and_BPH, WCE, NST, BPH;
  char targ[10];
  int rand_a;
  //  double lx, ly, lz;
  //int typ_e=-1;
  
  //parser
  int opt;
  char enefilename[MAXCHAR], ernfilename[MAXCHAR], trjfilename[MAXCHAR];
  sprintf(ernfilename, "");
  sprintf(enefilename,"");
  sprintf(trjfilename,"");

  int sflag=0,aflag=0,tflag=0,fflag=0,rflag=0,eflag=0,iflag=0,Tflag=0;
  while ((opt = getopt(argc, argv, "atfrse:i:T:")) != -1) {
    switch (opt)
      {
      case 'a':
	aflag=1;
	break;
      case 't':
	tflag=1;
	break;
      case 'f':
	fflag=1;
	break;
      case 'r':
	rflag=1;
	break;
      case 's':
	sflag=1;
	break;
      case 'e':
        sprintf(ernfilename, "%s", optarg);
	eflag=1;
        break;
      case 'i':
        sprintf(enefilename, "%s", optarg);
	iflag=1;
        break;
      case 'T':
        sprintf(trjfilename, "%s", optarg);
	Tflag=1;
        break;
      default:
	show_usage();
      }
  }

  if(aflag+tflag+fflag+rflag+sflag>1){
    fprintf(stderr, "Invalid combination of options.\n");
    show_usage();
    exit(ERR_INPUT);
  }
  
  /* initialization */
  mc_n=(int)NSOLUTE;
  int nt_n=mc_n/N_PARTS_PER_NT;
  
  int len=strlen(enefilename);
  char tmpline[1000];
  strcpy(tmpline, trjfilename);
  //now look for the pXX signature in the string
  char *init;
  init=strtok(tmpline, ".");
  init=strtok(NULL, ".");
  
  const char *last_three=&enefilename[len-3];
  if(!strcmp(last_three, ".mc")){
    MC_read_params(&mc_iter, &rand_a,-1);
    tempinit=MC_read_checkpoint(&mc_n, &rx, &ry, &rz, &rand_a, -1, enefilename, 0, NULL);
    MC_initialize_energies(mc_n, rx, ry, rz);
  }
  else if(!strcmp(last_three, "pdb")){
    MC_read_nsolute(&mc_n,-1, enefilename);
    MC_read_params(&mc_iter, &rand_a,-1);
    //MC_initialize(mc_n, &rx, &ry, &rz, lx, ly, lz, READ_MC_CONF, rand_a, mc_ini);
    MC_initialize_global(mc_n, rand_a,-1);
    MC_initialize_arrays(mc_n,&rx,&ry,&rz);
    MC_read_pdb(mc_n, &rx, &ry, &rz, -1, enefilename);
    MC_initialize_energies(mc_n, rx, ry, rz);
  }
  else{
    printf("Error opening file %s. Its extension must be pdb or mc.\n", enefilename);
    exit(ERR_INPUT);
  }
  /* before running */
  nt_n=mc_n/N_PARTS_PER_NT;
  if(tflag){
    TOT=MC_get_energy(nt_n, rx, ry, rz, 3);
    BB1=MC_get_energy(nt_n, rx, ry, rz, 0);
    WCE_and_BPH=MC_get_energy(nt_n, rx, ry, rz, 4);
    for(base=0;base<N_BASES;base++){
      for(face=0;face<WC_FACES;face++){
	nb_bp_well[base][face]=0;
      }
      nb_bp_spec_well[base][0]=0;
      nb_bp_spec_well[base][1]=0;
    }
    WCE=MC_get_energy(nt_n, rx, ry, rz, 4);
    BPH=WCE_and_BPH-WCE;
    NST=MC_get_energy(nt_n, rx, ry, rz, 2);
    ST_and_BB2=MC_get_energy(nt_n, rx, ry, rz, 1);
    for(base=0;base<N_BASES_SQ;base++){
      for(stf=0;stf<N_STFACES;stf++)
	b_st_well[base][stf]=0;
      nb_st_well[base]=0;
    }
    BB2=MC_get_energy(nt_n, rx, ry, rz, 1);
    printf("TOTAL ENERGY \t\t: %lf\nBackbone E    \t\t: %lf\nStacking(B) E     \t: %lf\nStacking(NB) E \t\t: %lf\nBase-pairing(WC) E\t: %lf\nBase-Phosphate(BPH) E  \t: %lf\n",
	   TOT,
	   BB1+BB2,
	   ST_and_BB2-BB2,
	   NST,
	   WCE,
	   BPH);
  }
  else if(sflag)
    MC_print_secondary_structure(nt_n, rx, ry, rz);
  else if(rflag)
    MC_print_radius_of_gyration(nt_n, rx, ry, rz);
  else if(aflag)
    MC_print_contact_list(nt_n, rx, ry, rz);
  else if(fflag)
    printf("Total E       : %lf\n    ",MC_get_energy(nt_n, rx, ry, rz, 3) );
  else if(eflag && Tflag){
    FILE *MCTRAJ;
    //FILE *TRTYCTCS;
    int trj=0;
    if((MCTRAJ=fopen(trjfilename, "rb"))==NULL){
      printf("Unable to read trajectory file %s\n", trjfilename);
      exit(ERR_INPUT);
    }
    //if((TRTYCTCS=fopen("interloop_ctcs.dat", "w"))==NULL){
    //  printf("Unable to open file for writing tertiary contacts.\n");
    //  exit(ERR_INPUT);
    //}
    printf("Reading Ernwin loops...\n");
    read_ernwin_loops(ernfilename, &NLOOPS, &LOOPNNT, &LOOPINDXS,&LOOPLABELS);
        
    
    //now read next config and analyze
    printf("Opening trajectory...\n");
    while(!MC_read_trajectory(mc_n, rx, ry, rz, MCTRAJ, trj)){
      MC_calc_and_write_ernwin_contacts(nt_n, rx, ry, rz, NLOOPS, LOOPNNT, LOOPINDXS, LOOPLABELS, trj, init);
      trj++;
    }
    fclose(MCTRAJ);
    //fclose(TRTYCTCS);
  }
  
  else if(!strcmp(enefilename, "") ){
    show_usage();
    exit(ERR_INPUT);
  }
  
  else{
    show_usage();
    exit(ERR_INPUT);
  }
  //if(sflag) printf("No clashes found.\n");
  printf("No clashes found.\n");
  return 0;
}


void read_ernwin_loops(char *filename, int *nl, int **nnt, int ***indexes, char ***labels){
  FILE *ernwincoord=fopen(filename, "r");
  char * line = NULL;
  int tnl=0;
  size_t len = 0;
  ssize_t read;
  char cop[6], cop2[6];
  char *token;
  int nt1,nt2,nt3,nt4,looplen;
  int ii;
  while ((read = getline(&line, &len, ernwincoord)) != -1) {
    if(read > 6){
      strncpy(cop, line, 6*sizeof(char));
      if(!strcmp(cop,"define")){
	token=strtok(line, " ");
	token=strtok(NULL, " ");
	if(token[0]=='s' || token[0]=='i' || token[0]=='h' || token[0]=='h'|| token[0]=='f' || token[0]=='t'){
	  tnl++;
	}
      }    }  }
  *nl=tnl;
  *nnt=(int *)malloc(sizeof(int)*tnl);
  *indexes=(int **)malloc(sizeof(int*)*tnl);
  *labels=(char **)malloc(sizeof(char*)*tnl);
  for(ii=0;ii<tnl;ii++) (*labels)[ii]=(char *)malloc(sizeof(char)*2);
  
  fclose(ernwincoord);

  tnl=0;
  ernwincoord=fopen(filename, "r");
  while ((read = getline(&line, &len, ernwincoord)) != -1) {
    if(read > 6){
      strncpy(cop, line, 6*sizeof(char));
      if(!strcmp(cop,"define")){
	//printf("line is %s\n", line);
	token=strtok(line, " ");
	token=strtok(NULL, " ");
	//printf("token is %s\n", token);
	if(token[0]=='s'){
	  //|| token[0]=='i'){
	  //we have a stem
	  (*labels)[tnl][0]=token[0];(*labels)[tnl][1]=token[1];
	  //token=strtok(line, " ");
	  token=strtok(NULL, " ");
	  nt1=atoi(token);
	  //token=strtok(line, " ");
	  token=strtok(NULL, " ");
	  nt2=atoi(token);
	  token=strtok(NULL, " ");
	  nt3=atoi(token);
	  token=strtok(NULL, " ");
	  nt4=atoi(token);
	  //printf("%d %d %d %d\n", nt1, nt2, nt3, nt4);
	  looplen=nt2-nt1+1 + nt4-nt3+1;
	  (*nnt)[tnl]=looplen;
	  (*indexes)[tnl]=(int *)malloc(sizeof(int)*looplen);
	  for(ii=0;ii<=nt2-nt1;ii++)
	    (*indexes)[tnl][ii]=ii+nt1-1;
	  for(ii=0;ii<=nt4-nt3;ii++)
	    (*indexes)[tnl][ii+nt2-nt1+1]=ii+nt3-1;
	  tnl++;
	}
	else if(token[0]=='i'){
	  //we have an internal loop
	  //printf("FOUND INTERNAL\n");
	  (*labels)[tnl][0]=token[0];(*labels)[tnl][1]=token[1];
	  token=strtok(NULL, " ");
	  nt1=atoi(token);
	  token=strtok(NULL, " ");
	  nt2=atoi(token);
	  looplen=nt2-nt1+1;
	  
	  token=strtok(NULL, " ");
	  if(token!=NULL){
	    nt3=atoi(token);
	    token=strtok(NULL, " ");
	    nt4=atoi(token);
	    looplen=looplen + nt4-nt3+1;
	  }
	  (*nnt)[tnl]=looplen;
	  
	  (*indexes)[tnl]=(int *)malloc(sizeof(int)*looplen);
	  for(ii=0;ii<=nt2-nt1;ii++)
	    (*indexes)[tnl][ii]=ii+nt1-1;
	  if(looplen>nt2-nt1+1){
	    for(ii=0;ii<=nt4-nt3;ii++)
	      (*indexes)[tnl][ii+nt2-nt1+1]=ii+nt3-1;
	  }
	  tnl++;
	}
	else if(token[0]=='f' || token[0]=='h' || token[0]=='t'){
	  //we have a hairpin, first or terminal loop. 
	  //printf("FOUND HAIRPIN\n");
	  (*labels)[tnl][0]=token[0];(*labels)[tnl][1]=token[1];
	  token=strtok(NULL, " ");
	  nt1=atoi(token);
	  token=strtok(NULL, " ");
	  nt2=atoi(token);
	  looplen=nt2-nt1+1;
	  (*nnt)[tnl]=looplen;
	  (*indexes)[tnl]=(int *)malloc(sizeof(int)*looplen);
	  for(ii=0;ii<=nt2-nt1;ii++)
	    (*indexes)[tnl][ii]=ii+nt1-1;
	  tnl++;
	}
	//nl++;
      }  }  }
  fclose(ernwincoord);
}
