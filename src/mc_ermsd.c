#include "mc_ermsd.h"

void MC_copy_single_ermsd_g(int trial, int nt_c, int nt_n){
  
    int i,j,d;
    for(i=0;i<nt_n;i++){
      for(d=0;d<DIM+1;d++){
	G_trial[i][nt_c][d]=G_curr[i][nt_c][d];
	G_trial[nt_c][i][d]=G_curr[nt_c][i][d];
      }
    }
}

void MC_copy_ermsd_g(int trial, int nt_n){
  int i,j,d;
  //printf("copying g curr to g trial\n");
  for(i=0;i<nt_n;i++)
    for(j=0;j<nt_n;j++)
      for(d=0;d<DIM+1;d++){
	G_trial[i][j][d]=G_curr[i][j][d];
	G_trial[j][i][d]=G_curr[j][i][d];
      }
  
}

void MC_update_ermsd_g(int nt_c, int nt_n){
  int i,j,d;
  for(i=0;i<nt_n;i++){
    for(d=0;d<DIM+1;d++){
      G_curr[i][nt_c][d]=G_trial[i][nt_c][d];
      G_curr[nt_c][i][d]=G_trial[nt_c][i][d];
    }
  }
}

double MC_reduce_ermsd(int nt_c, int nt_n){
  double ret=0;
  int i, d;
  for(i=0;i<nt_n;i++)
    if(G_groups[i][nt_c]>-1)
      for(d=0;d<DIM+1;d++){
      	ret+=SQ(G_ref[i][nt_c][d]-G_curr[i][nt_c][d]);
	ret+=SQ(G_ref[nt_c][i][d]-G_curr[nt_c][i][d]);
      }
  return ret;
}

double MC_reduce_ermsd_trial(int nt_c, int nt_n){
  double ret=0;
  int i, d;
  //printf("reducing trial...");
  for(i=0;i<nt_n;i++)
    if(G_groups[i][nt_c]>-1)
      for(d=0;d<DIM+1;d++){
	ret+=SQ(G_ref[i][nt_c][d]-G_trial[i][nt_c][d]);
	ret+=SQ(G_ref[nt_c][i][d]-G_trial[nt_c][i][d]);
      }
  //printf("ret = %lf\n", ret);
  return ret;
}

double ermsd_norm(double *r){
  return sqrt(SQ(r[0]/ERMSDX)+SQ(r[1]/ERMSDY)+SQ(r[2]/ERMSDZ));
}

double get_ermsd(){
  return sqrt(ERMSD_SQ);
}

double get_first_ermsd(double **rx, double **ry, double **rz, int nt_n){
  double ret, sum=0.0;
  double vec[DIM];
  //mol can be trial or current
  //ermsd is always calculated with respect to the reference
  int i,j,d;
  //build_ermsd_g(&ermsdref_X, &ermsdref_Y, &ermsdref_Z, G_ref, nt_n);
  build_ermsd_g(rx, ry, rz, G_curr, nt_n);
  build_ermsd_g(rx, ry, rz, G_trial, nt_n); //since mc_temp_x,y,z are not completely initialized
  //  build_ermsd_g(mc_temp_x, mc_temp_y, mc_temp_z, G_trial, nt_n);
  //reduce the matrices
  for(i=0;i<nt_n;i++)
    for(j=0;j<nt_n;j++){
      if(G_groups[i][j]>-1)
	for(d=0;d<DIM+1;d++){
	  sum+=SQ(G_ref[i][j][d]-G_curr[i][j][d]);
	}
    }
  //ret=sum/nt_n;
  ret=sum/ERMSD_NNT;
  return ret;
}

double get_current_ermsd(double **rx, double **ry, double **rz, int nt_n){
  double ret, sum=0.0;
  double vec[DIM];
  //ermsd is always calculated with respect to the reference
  int i,j,d;
  build_ermsd_g(rx, ry, rz, G_curr, nt_n);
  //reduce the matrices
  for(i=0;i<nt_n;i++)
    for(j=0;j<nt_n;j++){
      if(G_groups[i][j]>-1)
	for(d=0;d<DIM+1;d++){
	  sum+=SQ(G_ref[i][j][d]-G_curr[i][j][d]);
	}
    }
  ret=sum/ERMSD_NNT;
  return ret;
}

void build_ermsd_g(double **x, double **y, double **z, double ***G_res, int nt_n){
  int i, j;
  for(i=0;i<nt_n;i++)
    for(j=0;j<nt_n;j++)
      get_ermsd_g_pair(i,j,x, y, z, G_res);
}

void get_ermsd_g_pair(int nt1, int nt2, double **x, double **y, double **z, double ***G_res){
  double pref;
  int d;
  int at1=nt1*N_PARTS_PER_NT;int at2=nt2*N_PARTS_PER_NT;
  double rdist, gammaR, ret=0;
  double rproj[DIM];
  double rvec[DIM];
  rvec[0]=(*x)[at2]-(*x)[at1];rvec[1]=(*y)[at2]-(*y)[at1];rvec[2]=(*z)[at2]-(*z)[at1];
  //project first
  rproj[0]=0.0;rproj[1]=0.0;rproj[2]=0.0;
  if(nt1!=nt2)
    proj_on_nt(rvec, *x, *y, *z, nt1, rproj);
  
  MC_assign_ermsd_g_pair(nt1, nt2, rproj, G_res);
}

void MC_assign_ermsd_g_pair(int nt1, int nt2, double *r_vec, double ***G_res){
  double edist=ermsd_norm(r_vec);
  double pref;
  if(edist<ERMSD_CUTOFF && nt1!=nt2){
    pref=sin(M_PI*edist/ERMSD_CUTOFF)/edist;
    G_res[nt1][nt2][0]=(ERMSD_CUTOFF/M_PI)*pref*r_vec[0]/ERMSDX;
    G_res[nt1][nt2][1]=(ERMSD_CUTOFF/M_PI)*pref*r_vec[1]/ERMSDY;
    G_res[nt1][nt2][2]=(ERMSD_CUTOFF/M_PI)*pref*r_vec[2]/ERMSDZ;
    G_res[nt1][nt2][3]=(ERMSD_CUTOFF/M_PI)*(1+cos(M_PI*edist/ERMSD_CUTOFF));
  }
  else{
    G_res[nt1][nt2][0]=0;
    G_res[nt1][nt2][1]=0;
    G_res[nt1][nt2][2]=0;
    G_res[nt1][nt2][3]=0;
  }
  //printf("%d %d   %lf %lf %lf %lf\n", nt1, nt2, G_res[nt1][nt2][0], G_res[nt1][nt2][1], G_res[nt1][nt2][2], G_res[nt1][nt2][3]);
}

double MC_get_pair_ermsd(double ex, double ey, double ez, double Gref_x, double Gref_y, double Gref_z, double Gref_w){
  double ret=0.0;
  double pref1, pref2;
  double edist=sqrt(SQ(ex)+SQ(ey)+SQ(ez));
  double gx=0.0, gy=0.0, gz=0.0, gw=0.0;
  if(edist<ERMSD_CUTOFF){
    //if(edist<ERMSD_CUTOFF && nt1!=nt2){
    pref1=sin(M_PI*edist/ERMSD_CUTOFF)/edist;
    pref2=pref1*(ERMSD_CUTOFF/M_PI);
    //G_res[nt1][nt2][0]=(ERMSD_CUTOFF/M_PI)*pref*ex;
    //G_res[nt1][nt2][1]=(ERMSD_CUTOFF/M_PI)*pref*ey;
    //G_res[nt1][nt2][2]=(ERMSD_CUTOFF/M_PI)*pref*ez;
    //G_res[nt1][nt2][3]=(ERMSD_CUTOFF/M_PI)*(1+cos(M_PI*edist/ERMSD_CUTOFF));
    gx=pref2*ex;
    gy=pref2*ey;
    gz=pref2*ez;
    gw=(ERMSD_CUTOFF/M_PI)*(1+cos(M_PI*edist/ERMSD_CUTOFF));
    //ret=SQ(Gref_x-pref2*ex) + SQ(Gref_y-pref2*ey) + SQ(Gref_z-pref2*ez) + SQ(Gref_w-(ERMSD_CUTOFF/M_PI)*(1+cos(M_PI*edist/ERMSD_CUTOFF)));
  }
  ret=SQ(Gref_x-gx) + SQ(Gref_y-gy) + SQ(Gref_z-gz) + SQ(Gref_w-gw);
  return ret;
}

void MC_init_ermsd_restr(int nt_n){
  FILE *ermsdfile;
  char filename[256];
  int i, mp, j, d, itemp, iN, pind, pgroup, ptype, ntemp=0, pair;
  double dtemp1, dtemp2, dx, dy, dz;
  char ctemp;
  double ervec[DIM];
  ERMSD_FLAG=-1;
  ERMSD_PREF=0;
  ERMSD_N_GROUPS=1;
  ERMSD_NNT=nt_n;
  int do_init=0;
  
  printf("Initializing ERMSD arrays.\n");
  G_ref  =(double ***)malloc(sizeof(double**)*nt_n);
  G_curr =(double ***)malloc(sizeof(double**)*nt_n);
  G_trial=(double ***)malloc(sizeof(double**)*nt_n);
  G_groups=(int **)malloc(sizeof(int *)*nt_n);
  for(i=0;i<nt_n;i++){
    G_ref[i]  =(double**)malloc(sizeof(double *)*nt_n);
    G_curr[i] =(double**)malloc(sizeof(double *)*nt_n);
    G_trial[i]=(double**)malloc(sizeof(double *)*nt_n);
    G_groups[i]=(int *) malloc(sizeof(int)*nt_n);
  }
  for(i=0;i<nt_n;i++){
    for(j=0;j<nt_n;j++){
      G_ref[i][j]  =(double*)malloc(sizeof(double *)*(DIM+1));
      G_curr[i][j] =(double*)malloc(sizeof(double *)*(DIM+1));
      G_trial[i][j]=(double*)malloc(sizeof(double *)*(DIM+1));
    }
  }
  //initialize to zero!!
  for(i=0;i<nt_n;i++)
    for(j=0;j<nt_n;j++){
      for(d=0;d<DIM+1;d++){
	G_ref[i][j][d]=0.0;
	G_curr[i][j][d]=0.0;
	G_trial[i][j][d]=0.0;
      }
      G_groups[i][j]=-1;
    }
  ermsdref_X=(double *)malloc(sizeof(double)*nt_n*N_PARTS_PER_NT);
  ermsdref_Y=(double *)malloc(sizeof(double)*nt_n*N_PARTS_PER_NT);
  ermsdref_Z=(double *)malloc(sizeof(double)*nt_n*N_PARTS_PER_NT);
  
  char *lline=NULL, *pdbrectyp, *tmp, *tmp2, *stmp, *cpline;
  int l, ll, at=0;
  char s1[MAXSTR], s2[MAXSTR], s3[MAXSTR];
  static size_t st_l=0;
  

  int group, *nnt_group, ntg, nt;//, grflag;
  int **ntind_group;
  
  sprintf(filename, "ermsd_frags.lst");
  if((ermsdfile=fopen(filename, "r"))==NULL){
    printf("ERMSD: No reference structure found for ERMSD pulling %s. No pulling.\n", filename);
  }
  else{
    //HERE WE DEAL WITH THE FRAGMENTS
    printf("ERMSD: found structures to steer!\n");
    l=getline(&lline, &st_l, ermsdfile);
    sscanf(lline,"%s %s %s %d %lf %lf", s1, s2, s3,  &itemp, &dtemp1, &dtemp2);
    if(!strcmp(s1, "REMARK") && !strcmp(s2, "ERMSD") && !strcmp(s3, "PARAMS")){
      ERMSD_N_GROUPS=itemp;
      ERMSD_PREF=dtemp1;
      ERMSD_CUTOFF=dtemp2;
      printf("ERMSD: %d fragments. K=%lf with cutoff=%lf\n", ERMSD_N_GROUPS, ERMSD_PREF, ERMSD_CUTOFF);
      ntind_group=(int **) malloc(sizeof(int *)*ERMSD_N_GROUPS);
      nnt_group=(int *)malloc(sizeof(int)*ERMSD_N_GROUPS);
      //read groups
      for(group=0;group<ERMSD_N_GROUPS;group++){
	l=getline(&lline, &st_l, ermsdfile);
	cpline=strndup(lline, MAXSTR);
	sscanf(cpline, "%s %s %s", s1, s2, s3);
	if(!(!strcmp(s1, "REMARK") && !strcmp(s2, "ERMSD") && !strcmp(s3, "GROUP"))){printf("Wrong syntax in ERMSD file ermsd_frags.lst\n");exit(ERR_INPUT);}
	nnt_group[group]=0;
	tmp=strtok(cpline, " ");while(tmp!=NULL){if(atoi(tmp)!=0 ||(atoi(tmp)==0 && tmp[0]=='0')) {nnt_group[group]++;} tmp=strtok(NULL, " ");}free(cpline);
	ntind_group[group]=(int *)malloc(sizeof(int)*nnt_group[group]);
	ntemp+=nnt_group[group];
	ntg=0;
	tmp=strtok(lline, " ");while(tmp!=NULL){if(atoi(tmp)!=0 ||(atoi(tmp)==0 && tmp[0]=='0')) {ntind_group[group][ntg]=atoi(tmp);ntg++;} tmp=strtok(NULL, " ");};
      }	
      //read coordinates after groups are defined
      for(group=0;group<ERMSD_N_GROUPS;group++){
	for(nt=0;nt<nnt_group[group];nt++){
	  for(at=0;at<N_PARTS_PER_NT;at++){
	    l=getline(&lline, &st_l, ermsdfile);
	    cpline=strndup(lline, 6);
	    if(!strcmp(cpline, "ATOM  ")){
	      stmp=strndup(lline+30, 8);dx=atof(stmp);free(stmp);
	      stmp=strndup(lline+38, 8);dy=atof(stmp);free(stmp);
	      stmp=strndup(lline+46, 8);dz=atof(stmp);free(stmp);
	      ermsdref_X[ntind_group[group][nt]*N_PARTS_PER_NT+at]=dx;
	      ermsdref_Y[ntind_group[group][nt]*N_PARTS_PER_NT+at]=dy;
	      ermsdref_Z[ntind_group[group][nt]*N_PARTS_PER_NT+at]=dz;
	      //printf("frags %d  %lf %lf %lf\n", at, dx, dy, dz);
	    } else {printf("Wrong syntax in ERMSD file ermsd_frags.lst\n");exit(ERR_INPUT);}
	    free(cpline);
	  }
	}
	for(i=0;i<nnt_group[group];i++)
	  for(j=0;j<nnt_group[group];j++){
	    G_groups[ntind_group[group][i]][ntind_group[group][j]]=group;
	    get_ermsd_g_pair(ntind_group[group][i],ntind_group[group][j],&ermsdref_X, &ermsdref_Y, &ermsdref_Z, G_ref);
	  }
      }
    }
    else{
      printf("Wrong syntax in ERMSD file ermsd_frags.lst\n");
      exit(ERR_INPUT);
    }
    fclose(ermsdfile);
    for(group=0;group<ERMSD_N_GROUPS;group++)
      free(ntind_group[group]);
    free(ntind_group);
    free(nnt_group);
  }
  ERMSD_NNT=ntemp;
  //for(group=0;group<ERMSD_N_GROUPS;group++){
  //printf("group %d\n", group);
  /* for(i=0;i<nt_n;i++) */
  /*   for(j=0;j<nt_n;j++) */
  /*     printf("%d  %d   %d   %lf %lf %lf\n", i,j,G_groups[i][j], G_ref[i][j][0], G_ref[i][j][1], G_ref[i][j][2]); */
  /* exit(1); */
}

void MC_get_ermsd_pair_type(int i, int j, double *ervec, int ptype){
  if(mc_types[i*N_PARTS_PER_NT]==TYP_ADENINE && mc_types[j*N_PARTS_PER_NT]==TYP_URACIL){
    ervec[0]=ERMSD_AU_x;ervec[1]=ERMSD_AU_y;ervec[2]=ERMSD_AU_z;}
  else if(mc_types[i*N_PARTS_PER_NT]==TYP_URACIL  && mc_types[j*N_PARTS_PER_NT]==TYP_ADENINE){
      ervec[0]=ERMSD_UA_x;ervec[1]=ERMSD_UA_y;ervec[2]=ERMSD_UA_z;}
  else if(mc_types[i*N_PARTS_PER_NT]==TYP_GUANINE  && mc_types[j*N_PARTS_PER_NT]==TYP_CYTOSINE){
    ervec[0]=ERMSD_GC_x;ervec[1]=ERMSD_GC_y;ervec[2]=ERMSD_GC_z;}
  else if(mc_types[i*N_PARTS_PER_NT]==TYP_CYTOSINE && mc_types[j*N_PARTS_PER_NT]==TYP_GUANINE){
    ervec[0]=ERMSD_CG_x;ervec[1]=ERMSD_CG_y;ervec[2]=ERMSD_CG_z;}
  else{
    printf("ERMSD: type for constraint not found!!\n");
    exit(1);
  }
}

void MC_init_ermsd_out(int proc){
  char filename[256];
  sprintf(filename, "ermsd_obs.p%02d.dat", proc);
  if((ermsd_obs=fopen(filename, "w"))==NULL){
    printf("ERMSD: can't create file %s.\n", filename);
  }
  else {
    fprintf(ermsd_obs, "# time internal_energy  ermsd total_energy\n");
  }
}

void MC_write_ermsd_obs(int step, double energ){
  double ermsd=get_ermsd();
  fprintf(ermsd_obs, "%d %lf %lf %lf\n", step, energ-0.5*ERMSD_PREF*SQ(ermsd),ermsd, energ); 
}
