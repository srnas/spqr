#include "mc_ermsd.h"


#ifdef LNKRMV

void MC_set_linked_loops(int nt_n, double *posx, double *posy, double *posz){
  //we read a file like the ermsd frags
  int i,j,d;
  int tmpout;
  char *ptmpout;
  FILE *lnkdloopsfile;
  char filename[256];
  char ltyp[MAXSTR], loopline[MAXSTR];
  int ini1,ini2,end1,end2,loo1,loo2;
  in_link=(int*)malloc(sizeof(int)*nt_n);
  for(i=0;i<nt_n;i++)in_link[i]=-1;
  
  sprintf(filename, "linked_loops.lst");
  if((lnkdloopsfile=fopen(filename, "r"))==NULL){
    printf("PULL AWAY LINKS: No reference structure found for LINK pulling %s. No pulling.\n", filename);
    //for(i=0;i<nt_n;i++)
    // printf("%d  %d\n", i, in_link[i]);
  }
  else{
    ptmpout=fgets(loopline,MAXSTR,lnkdloopsfile);
    sscanf(loopline, "%d %d %lf %lf", &mc_N_links, &mc_N_loops, &LOOP_K_cmcm, &LOOP_K_cmclp);
    
    mc_links=(int **)malloc(sizeof(int *)*mc_N_links);
    mc_loops=(int **)malloc(sizeof(int *)*mc_N_loops);
    mc_loop_type=(int *)malloc(sizeof(int)*mc_N_loops);
    mc_loop_size=(int *)malloc(sizeof(int)*mc_N_loops);
    mc_loop_CM=(double **)malloc(sizeof(double)*mc_N_loops);
    mc_loop_clpair=(double **)malloc(sizeof(double)*mc_N_loops);

    for(i=0;i<mc_N_links;i++) mc_links[i]=(int *)malloc(sizeof(int)*2);
    for(i=0;i<mc_N_loops;i++){
      ptmpout=fgets(loopline,MAXSTR,lnkdloopsfile);
      sscanf(loopline,"%s",ltyp);
      //printf("i=%d  typ = %s\n", i,ltyp);
      if(!strcmp(ltyp,"hp")) mc_loop_type[i]=LOOP_HP;
      else if(!strcmp(ltyp,"st")) mc_loop_type[i]=LOOP_ST;
      else if(!strcmp(ltyp,"il")) mc_loop_type[i]=LOOP_IL;
      else {printf("Wrong loop type at file %s!\n", filename);exit(ERR_INPUT);}
      if(mc_loop_type[i]==LOOP_HP || mc_loop_type[i]==LOOP_IL) {sscanf(loopline,"%s %d %d", ltyp,  &ini1,&end1);mc_loop_size[i]=end1-ini1+1;}
      else{sscanf(loopline,"%s %d %d %d %d", ltyp, &ini1,&end1,&ini2,&end2); mc_loop_size[i]=end1-ini1+1+end2-ini2+1;}
      mc_loops[i]=(int*)malloc(sizeof(int)*mc_loop_size[i]);
      for(j=ini1;j<end1+1;j++){
	mc_loops[i][j-ini1]=j;
	in_link[j]++;
      }
      if(mc_loop_type[i]==LOOP_ST){
	for(j=ini2;j<end2+1;j++){
	  mc_loops[i][j-ini2+(end1+1-ini1)]=j;
	  in_link[j]++;
	}
	//printf("%d %d %d %d\n", ini1,end1,ini2,end2);
      }
      mc_loop_CM[i]=(double *)malloc(sizeof(double)*DIM);
      mc_loop_clpair[i]=(double *)malloc(2*sizeof(double)*DIM);
    }
    
    for(i=0;i<mc_N_links;i++){
      ptmpout=fgets(loopline,MAXSTR,lnkdloopsfile);
      sscanf(loopline,"%d %d", &loo1,&loo2);
      mc_links[i][0]=loo1;
      mc_links[i][1]=loo2;
    }
  }
  printf("Number of links: %d\n",mc_N_links);
  printf("Number of involved loops: %d\n--------------------\n",mc_N_loops);
  for(i=0;i<mc_N_links;i++) printf("Link %d : %d %d\n", i, mc_links[i][0],mc_links[i][1]);
  for(i=0;i<mc_N_loops;i++) {printf("Loop %d (type %d, size %d) : ", i,mc_loop_type[i],mc_loop_size[i]);for(j=0;j<mc_loop_size[i];j++) printf("%d ", mc_loops[i][j]);printf("\n");}
  
}

void mc_update_loop(int iloop, int nt_c, double *rx, double *ry, double *rz){
  int d, ntind, nt1,nt2,nt3,nt4;
  double pos1[DIM],pos2[DIM],pos3[DIM],pos4[DIM];
  for(d=0;d<DIM;d++){
    mc_loop_CM[iloop][d]=0;
    mc_loop_clpair[iloop][d]=0;
    mc_loop_clpair[iloop][2*d]=0;
  }
  /*CALCULATE CENTER OF MASS*/
  for(ntind=0;ntind<mc_loop_size[iloop];ntind++){
    nt1=mc_loops[iloop][ntind];
    if(nt1!=nt_c){
      mc_loop_CM[iloop][0]+=(rx[nt1*N_PARTS_PER_NT+IPHO]+rx[nt1*N_PARTS_PER_NT+ISUG]);
      mc_loop_CM[iloop][1]+=(ry[nt1*N_PARTS_PER_NT+IPHO]+ry[nt1*N_PARTS_PER_NT+ISUG]);
      mc_loop_CM[iloop][2]+=(rz[nt1*N_PARTS_PER_NT+IPHO]+rz[nt1*N_PARTS_PER_NT+ISUG]);
    }
    else{
      mc_loop_CM[iloop][0]+=(get_unf_coo_temp_x(nt1*N_PARTS_PER_NT+IPHO)+get_unf_coo_temp_x(nt1*N_PARTS_PER_NT+ISUG));
      mc_loop_CM[iloop][1]+=(get_unf_coo_temp_y(nt1*N_PARTS_PER_NT+IPHO)+get_unf_coo_temp_y(nt1*N_PARTS_PER_NT+ISUG));
      mc_loop_CM[iloop][2]+=(get_unf_coo_temp_z(nt1*N_PARTS_PER_NT+IPHO)+get_unf_coo_temp_z(nt1*N_PARTS_PER_NT+ISUG));
    }
  }
  for(d=0;d<DIM;d++)
    mc_loop_CM[iloop][d]/=(2.0*(double)mc_loop_size[iloop]);

  /*CALCULATE CLOSING PAIR*/
  if(mc_loop_type[iloop]==LOOP_HP){
    nt1=mc_loops[iloop][0];
    nt2=mc_loops[iloop][mc_loop_size[iloop]-1];
    get_real_sugpos(nt1,nt_c,rx,ry,rz,pos1);
    get_real_sugpos(nt2,nt_c,rx,ry,rz,pos2);
    for(d=0;d<DIM;d++)
      mc_loop_clpair[iloop][d]=0.5*(pos1[d]+pos2[d]);
  } else if(mc_loop_type[iloop]==LOOP_ST){
    nt1=mc_loops[iloop][0];
    nt2=mc_loops[iloop][mc_loop_size[iloop]-1];
    get_real_sugpos(nt1,nt_c,rx,ry,rz,pos1);
    get_real_sugpos(nt2,nt_c,rx,ry,rz,pos2);
    for(d=0;d<DIM;d++){
      mc_loop_clpair[iloop][d]=pos1[d];
      mc_loop_clpair[iloop][2*d]=pos2[d];
    }
  } else if(mc_loop_type[iloop]==LOOP_IL){
    nt1=mc_loops[iloop][0];
    nt2=mc_loops[iloop][mc_loop_size[iloop]/2-1];
    nt3=mc_loops[iloop][mc_loop_size[iloop]/2];
    nt4=mc_loops[iloop][mc_loop_size[iloop]-1];
    get_real_sugpos(nt1,nt_c,rx,ry,rz,pos1);
    get_real_sugpos(nt2,nt_c,rx,ry,rz,pos2);
    get_real_sugpos(nt3,nt_c,rx,ry,rz,pos3);
    get_real_sugpos(nt4,nt_c,rx,ry,rz,pos4);
    for(d=0;d<DIM;d++){
      mc_loop_clpair[iloop][d]=0.5*(pos1[d]+pos4[d]);
      mc_loop_clpair[iloop][2*d]=0.5*(pos2[d]+pos3[d]);
    }
  } else {
    printf("Unrecognized loop type.\n"); exit(ERR_INTEG);
  }
}

int nt_is_in_link(int ilink, int nt_c){
  int ret=0;
  int loop1,loop2;
  int i;
  loop1=mc_links[ilink][0];
  loop2=mc_links[ilink][1];
  for(i=0;i<mc_loop_size[loop1];i++){
    if(mc_loops[loop1][i]==nt_c){ret=1;
      continue;}
  }
  for(i=0;i<mc_loop_size[loop2];i++){
    if(mc_loops[loop2][i]==nt_c){ret=1;
      continue;}
  }
  return ret;
}

int nt_is_in_loop(int iloop, int nt_c){
  int ret=0;
  int loop1;
  int i;
  for(i=0;i<mc_loop_size[iloop];i++){
    if(mc_loops[iloop][i]==nt_c){
      ret=1;
      continue;
    }
  }
  return ret;
}

int nts_in_same_link_but_different_loops(int nt_a, int nt_b){
  int i, ret=0;
  for(i=0;i<mc_N_links;i++){
    if((nt_is_in_loop(mc_links[i][0],nt_a)==1 &&  nt_is_in_loop(mc_links[i][1],nt_b)==1 )  || ( nt_is_in_loop(mc_links[i][0],nt_b)==1 &&  nt_is_in_loop(mc_links[i][1],nt_a)==1 )){
      ret=1;
      continue;
    }
  }
  return ret;
}

void get_real_sugpos(int nt, int nt_c, double *rx, double *ry, double *rz, double *rpos){
  if(nt !=nt_c){
    rpos[0]=rx[nt*N_PARTS_PER_NT+ISUG];
    rpos[1]=ry[nt*N_PARTS_PER_NT+ISUG];
    rpos[2]=rz[nt*N_PARTS_PER_NT+ISUG];
  } else {
    rpos[0]=get_unf_coo_temp_x(nt*N_PARTS_PER_NT+ISUG);
    rpos[1]=get_unf_coo_temp_y(nt*N_PARTS_PER_NT+ISUG);
    rpos[2]=get_unf_coo_temp_z(nt*N_PARTS_PER_NT+ISUG);
  }
}

double calc_link_energy(int nt_c, double *rx, double *ry, double *rz){
  int i,d,ilink;
  int l1,l2;
  double energ=0.0, d1sq, d2sq;
  for(ilink=0;ilink<mc_N_links;ilink++){
    if(nt_is_in_link(ilink,nt_c)==1) {
      l1=mc_links[ilink][0];
      l2=mc_links[ilink][1];
      
      mc_update_loop(l1, nt_c, rx, ry, rz);
      mc_update_loop(l2, nt_c, rx, ry, rz);
      //HERE ADD THE ENERGY BETWEEN THE LOOPS
      d1sq=0;d2sq=0;
      for(d=0;d<DIM;d++){
	d1sq+=SQ(mc_loop_CM[l1][d]-mc_loop_clpair[l2][d]);
	d2sq+=SQ(mc_loop_CM[l2][d]-mc_loop_clpair[l1][d]);
      }
      energ+=-LOOP_K_cmclp*(sqrt(d1sq)+sqrt(d2sq));
      if(mc_loop_type[l1]!=LOOP_HP){
	d2sq=0;
	for(d=0;d<DIM;d++)
	  d2sq+=SQ(mc_loop_CM[l2][2*d]-mc_loop_clpair[l1][d]);
	energ+=-LOOP_K_cmclp*sqrt(d2sq);
      }
      if(mc_loop_type[l2]!=LOOP_HP){
	d1sq=0;
	for(d=0;d<DIM;d++)
	  d1sq+=SQ(mc_loop_CM[l1][d]-mc_loop_clpair[l2][2*d]);
	energ+=-LOOP_K_cmclp*sqrt(d1sq);
      }
      if(mc_loop_type[l1]==LOOP_ST || mc_loop_type[l2]==LOOP_ST){
	d1sq=0;
	for(d=0;d<DIM;d++)
	  d1sq+=SQ(mc_loop_CM[l1][d]-mc_loop_CM[l2][d]);
	energ+=-LOOP_K_cmcm*sqrt(d1sq);
      }
    }
  }
  return energ;
}

#endif

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

double ermsd_norm(double *r){
  return sqrt(SQ(r[0]/ERMSDX)+SQ(r[1]/ERMSDY)+SQ(r[2]/ERMSDZ));
}

double get_ermsd(){
  return sqrt(ERMSD_SQ);
}

double get_ermsd_energ(){
  return ERMSD_ENERG;
}


double get_first_ermsd(double **rx, double **ry, double **rz, int nt_n, double *termsdsq, double *termsd_energ){
  double ret, sum=0.0, eesum=0.0,termsd,tpref;
  
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
	  termsd=SQ(G_ref[i][j][d]-G_curr[i][j][d]);
	  sum+=termsd;
	  eesum+=termsd*ERMSD_PREF[G_groups[i][j]];
	}
    }
  *termsdsq=sum/ERMSD_NNT;
  *termsd_energ=eesum;
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
  //*ERMSD_PREF=(double *)malloc(sizeof(double));
  //ERMSD_PREF[0]=0.0;
  ERMSD_N_GROUPS=1;
  ERMSD_NNT=nt_n;
  int do_init=0;
  int ssrflag=0;
  int ssind1=-1,ssind2=-1;
  double ERMSD_SSPREF=1.0;
  
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
  int tt;

  int group, *nnt_group, ntg, nt;//, grflag;
  int **ntind_group;
  ERMSD_SSTRUCT=(double **) malloc(sizeof(double *)*nt_n);
  for(i=0;i<nt_n;i++)
    ERMSD_SSTRUCT[i]=(double *)malloc(sizeof(double)*nt_n);
  for(i=0;i<nt_n;i++)
    for(j=0;j<nt_n;j++)
      ERMSD_SSTRUCT[i][j]=1.0;
  
  sprintf(filename, "ermsd_frags.lst");
  if((ermsdfile=fopen(filename, "r"))==NULL){
    printf("ERMSD: No reference structure found for ERMSD pulling %s. No pulling.\n", filename);
    ERMSD_PREF=(double *)malloc(sizeof(double));
    ERMSD_PREF[0]=0.0;
  }
  else{
    //HERE WE DEAL WITH THE FRAGMENTS
    printf("ERMSD: found structures to steer!\n");
    l=getline(&lline, &st_l, ermsdfile);
    sscanf(lline,"%s %s %s %d %lf", s1, s2, s3,  &itemp, &dtemp2);
    if(!strcmp(s1, "REMARK") && !strcmp(s2, "ERMSD") && !strcmp(s3, "PARAMS")){
      ERMSD_N_GROUPS=itemp;
      //ERMSD_PREF=dtemp1;
      ERMSD_PREF=(double *)malloc(sizeof(double)*ERMSD_N_GROUPS);
      //ERMSD_PREF[0]=0.0;
      ERMSD_CUTOFF=dtemp2;

            
      printf("ERMSD: %d fragments. Cutoff=%lf\n", ERMSD_N_GROUPS, ERMSD_CUTOFF);
      ntind_group=(int **) malloc(sizeof(int *)*ERMSD_N_GROUPS);
      nnt_group=(int *)malloc(sizeof(int)*ERMSD_N_GROUPS);
      //read groups

      l=getline(&lline, &st_l, ermsdfile);
      cpline=strndup(lline, MAXSTR);
      sscanf(cpline, "%s %s %s %lf", s1, s2, s3, &dtemp1);
      if(!strcmp(s1, "REMARK") && !strcmp(s2, "ERMSD") && !strcmp(s3, "SSPAIRS")){
	ERMSD_SSPREF=dtemp1;
	tt=0;
	tmp=strtok(cpline, " ");while(tmp!=NULL){if(tt>3 && (atoi(tmp)!=0 ||(atoi(tmp)==0 && tmp[0]=='0'))) {if(tt%2==0) ssind1=atoi(tmp);else {ssind2=atoi(tmp); ERMSD_SSTRUCT[ssind1][ssind2]=ERMSD_SSPREF;ERMSD_SSTRUCT[ssind2][ssind1]=ERMSD_SSPREF;ssind1=-1;ssind2=-1;}} tmp=strtok(NULL, " ");tt++;}free(cpline);
      }
      else if(!strcmp(s1, "REMARK") && !strcmp(s2, "ERMSD") && !strcmp(s3, "GROUP")){
	ssrflag=1;
      }
      else{ printf("Wrong syntax in ERMSD file ermsd_frags.lst\n");exit(ERR_INPUT);}
            
      for(group=0;group<ERMSD_N_GROUPS;group++){
	if(ssrflag==0){
	  l=getline(&lline, &st_l, ermsdfile);
	  cpline=strndup(lline, MAXSTR);
	  sscanf(cpline, "%s %s %s %lf", s1, s2, s3, &dtemp1);
	  ssrflag=0;
	}
	ssrflag=0;
	if(!(!strcmp(s1, "REMARK") && !strcmp(s2, "ERMSD") && !strcmp(s3, "GROUP"))){printf("Wrong syntax in ERMSD file ermsd_frags.lst\n");exit(ERR_INPUT);}
	
	ERMSD_PREF[group]=dtemp1;
	nnt_group[group]=0;
	tt=0;
	tmp=strtok(cpline, " ");while(tmp!=NULL){if(tt>3 && (atoi(tmp)!=0 ||(atoi(tmp)==0 && tmp[0]=='0'))) {nnt_group[group]++;} tmp=strtok(NULL, " ");tt++;}free(cpline);
	ntind_group[group]=(int *)malloc(sizeof(int)*nnt_group[group]);
	ntemp+=nnt_group[group];
	ntg=0;
	tt=0;
	tmp=strtok(lline, " ");while(tmp!=NULL){if(tt>3 && (atoi(tmp)!=0 ||(atoi(tmp)==0 && tmp[0]=='0'))) {ntind_group[group][ntg]=atoi(tmp);ntg++;} tmp=strtok(NULL, " ");tt++;};
	//printf("nnt group = %d\n", nnt_group[group]);
	//for(tt=0;tt<nnt_group[group];tt++) printf("%d ", ntind_group[group][tt]);
	//printf("\n");
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
      printf("ERMSD: ENERGY GROUPS AND PREFACTORS\t");
      for(group=0;group<ERMSD_N_GROUPS;group++){
	printf("%d: %lf   ( ", group, ERMSD_PREF[group]);
	for(nt=0;nt<nnt_group[group];nt++)
	  printf("%d ",ntind_group[group][nt]); 
	printf(" ) ");
      }
      printf("\n");
      if(ERMSD_SSPREF==1){
	printf("No constrained pairs\n");
      } else{
	printf("Constrained pairs (K=%lf):  ", ERMSD_SSPREF);
	for(i=0;i<nt_n;i++)
	  for(j=0;j<nt_n;j++)
	    if(ERMSD_SSTRUCT[i][j]!=1)
	      printf("(%d %d)  ", i,j);
	printf("\n");
      }
      //for(group=0;group<ERMSD_N_GROUPS;group++)
      //	ERMSD_PREF[group]*(0.5/((double)nnt_group[group]));
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
  //printf("SECONDARY STRUCTURE CONSTRAINTS INVOLVED\n");
  //for(i=0;i<nt_n;i++)
  // for(j=0;j<nt_n;j++)
  //   if(ERMSD_SSTRUCT[i][j]!=1) printf("%d  %d\n", i,j);
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
    fprintf(ermsd_obs, "# time internal_energy  ermsd_energy total_energy wall_energy\n");
  }
}

void MC_write_ermsd_obs(int step, double energ){
  double ermsd=get_ermsd();
  fprintf(ermsd_obs, "%d %lf %lf %lf %lf\n", step, energ-ERMSD_ENERG-WALL_ENERG,ermsd, energ, WALL_ENERG); 
}

double MC_wall_energy(double px, double py, double pz){
  double ret=0;
  int w;
  for(w=0;w<N_WALLS;w++){
    if(wall_epsilon[w]>0){
      double dist=fabs(wall_A[w]*px+wall_B[w]*py+wall_C[w]*pz+wall_D[w])/wall_MODSQ[w];
      
      double d4=1.0/(SQ(SQ(dist)));
      if(WALL_TYPE[w]==0)
	ret+=-wall_epsilon[w]*exp(-dist/wall_sigma[w])/dist+d4*d4*d4;
      else if(WALL_TYPE[w]==1)
	ret+=wall_epsilon[w]*dist+d4*d4*d4;
      else{
	printf("ERROR: Wall type not recognized!\n");
	exit(ERR_INPUT);
	
      }
      //printf("%d %lf %lf %lf %lf %lf %lf\t%lf\n",w, wall_epsilon[w],wall_A[w], wall_B[w],wall_C[w],wall_D[w],wall_MODSQ[w],ret);
    }
  }
  //printf("%lf  %lf\n", dist, ret);
  
  return ret;
}
