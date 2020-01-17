#include "pt_impl.h"

void PT_read_params(int m_i, int m_tot, FILE ** ptfile){
  double temp1, temp2;
  char filename[256];
  char *pch, *cpline;
  sprintf(filename, "params.spqr");
  FILE *PT_PARAMS;
  double *temperatures;
  char tfname[256];
  temperatures=(double *)malloc(sizeof(double)*m_tot);
  int ll;
  if((PT_PARAMS=fopen(filename, "r"))==NULL){
    printf("PARALLEL TEMPERING: params.spqr file not found.\n");
    exit(ERR_INPUT);
  }
  
  else{
    int l, cnt=0;
    char s[MAXSTR], s2[MAXSTR], *line=NULL;
    static size_t st_l=0;
    while((l=getline(&line, &st_l, PT_PARAMS))>0){
      if(l>3)
	if(line[0]=='P' && line[1]=='T' && line[2]=='_'){
	  sscanf(line, "%s", s);
	  if(!strcmp(s, "PT_TEMPS")) {
	    cpline=strndup(line, MAXSTR);
	    pch=strtok(cpline, " \t");
	    for(ll=0;ll<m_tot;ll++){
	      pch=strtok(NULL," ");
	      temperatures[ll]=atof(pch);
	    }
	    free(cpline);
	  }
	}
    }
    
  }
  mc_target_temp=temperatures[m_i];
  sprintf(tfname, "pt_enetem.p%02d.dat", m_i);
  *ptfile=fopen(tfname,"w");
  free(temperatures);
}

void PT_write(int iter, double temp, double energ, FILE **ptfile){
  fprintf(*ptfile, "%d %lf %lf\n", iter, temp, energ);
}

void PT_end(FILE **ptfile){
  fclose(*ptfile);
}
