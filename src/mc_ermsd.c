#include "mc_ermsd.h"


#ifdef LNKRMV
/* void MC_init_def_phpull(int nt_n){ */
/*   int i; */
/*   /\* phpull_K=(double*)malloc(sizeof(double)*nt_n); *\/ */
/*   /\* phpull_p0=(double**)malloc(sizeof(double*)*nt_n); *\/ */
/*   /\* for(i=0;i<nt_n;i++){ *\/ */
/*   /\*   phpull_p0[i]=(double *)malloc(sizeof(double)*DIM); *\/ */
/*   /\* } *\/ */
/*   /\* for(i=0;i<nt_n;i++){ *\/ */
/*   /\*   phpull_K[i]=0; *\/ */
/*   /\* } *\/ */
/* } */

void MC_set_linked_loops(int nt_n, double *posx, double *posy, double *posz){
  //we read a file like the ermsd frags
  
  int i,cnt,d;
  int stloop,ini,end,p;
  double pdx, pdy, pdz,len;
  
  char *lline=NULL,  *tmp, *tmp2, *cpline;
  double KTEMP;
  int l,nl, ll, at=0;
  char s1[MAXSTR], s2[MAXSTR], s3[MAXSTR], s4[MAXSTR];
  static size_t st_l=0;
  int group, *nnt_group, ntg, nt;//, grflag;
  int N_LINKS,ntl,lind,exloops[8],inl;
  int typ1=-1, typ2=-1;
  FILE *lnkdloopsfile;
  char filename[256];
  int tempnl1=0,tempnl2=0;
  int start2;

  my_link=(int*)malloc(sizeof(int)*nt_n);
  my_loop=(int*)malloc(sizeof(int)*nt_n);
  for(i=0;i<nt_n;i++){my_link[i]=-1;my_loop[i]=-1;}
  sprintf(filename, "linked_loops.lst");
  if((lnkdloopsfile=fopen(filename, "r"))==NULL){
    printf("PULL AWAY LINKS: No reference structure found for LINK pulling %s. No pulling.\n", filename);
  }
  else{
    for(i=0;i<nt_n;i++){
      fr_is_mobile[i]=FR_MOB_FROZ;
      //HERE MAKE NUCLEOTIDES INVISIBLE
    }
    
    //HERE WE DEAL WITH THE FRAGMENTS
    printf("PULL LINKS AWAY: found structures to steer!\n");
    l=getline(&lline, &st_l, lnkdloopsfile);
    cpline=strndup(lline, MAXSTR);
    sscanf(cpline, "%s %s %s %s", s1, s2, s3, s4);
    //sscanf(lline, "%s %s %s", s1, s2, s3);
    if(!(!strcmp(s1, "REMARK") && !strcmp(s2, "LNKDLPS"))){printf("Wrong syntax in LNKDLPS file linked_loops.lst\n");exit(ERR_INPUT);}
    KTEMP=atof(s3);
    N_LINKS=atoi(s4);
    //nnt_group[group]=0;
    //tmp=strtok(cpline, " ");while(tmp!=NULL){if(atoi(tmp)!=0 ||(atoi(tmp)==0 && tmp[0]=='0')) {nnt_group[group]++;} tmp=strtok(NULL, " ");}free(cpline);
    //ntind_group[group]=(int *)malloc(sizeof(int)*nnt_group[group]);
    //ntemp+=nnt_group[group];
    //printf("tmp = %d \n", tmp);
    //printf("File open with %d links, K=%lf\n", N_LINKS,KTEMP);
    loop1=(int**)malloc(sizeof(int*)*N_LINKS);
    loop2=(int**)malloc(sizeof(int*)*N_LINKS);
    nnt_loop1=(int*)malloc(sizeof(int)*N_LINKS);
    nnt_loop2=(int*)malloc(sizeof(int)*N_LINKS);
    for(nl=0;nl<N_LINKS;nl++){
      l=getline(&lline,&st_l,lnkdloopsfile);
      cpline=strndup(lline,MAXSTR);
      for(d=0;d<8;d++) exloops[d]=-1;
      ntl=-1;
      tmp=strtok(cpline, " ");
      while(tmp!=NULL){
	if(ntl==-1){
	  //check type string
	  //HAIRPIN=0, STEM=1
	  if(tmp[0]=='h' && tmp[1]=='p') typ1=0;
	  if(tmp[0]=='s' && tmp[1]=='t') typ1=1;
	  if(tmp[2]=='h' && tmp[3]=='p') typ2=0;
	  if(tmp[2]=='s' && tmp[3]=='t') typ2=1;
	  ntl++;
	}
	if(ntl>-1){
	  if(atoi(tmp)!=0 ||(atoi(tmp)==0 && tmp[0]=='0')) {
	    exloops[ntl]=atoi(tmp);ntl++;
	  }}
	tmp=strtok(NULL, " ");
      }
      free(cpline);
      //case hairpir-hairpin or hairpin-duplex
      tempnl1=0,tempnl2=0;
      if(typ1==0) start2=2;
      if(typ1==1) start2=4;
      if(typ1==0 || typ1==1)
	for(inl=exloops[0];inl<=exloops[1];inl++)
	  tempnl1++;
      if(typ1==1)
	for(inl=exloops[2];inl<=exloops[3];inl++)
	  tempnl1++;
      nnt_loop1[nl]=tempnl1;
      /***/
      if(typ2==0 || typ2==1)
	for(inl=exloops[start2];inl<=exloops[start2+1];inl++)
	  tempnl2++;
      if(typ2==1)
	for(inl=exloops[start2+2];inl<=exloops[start2+3];inl++)
	  tempnl2++;
      nnt_loop2[nl]=tempnl2;
      /***/
      loop1[nl]=(int*)malloc(sizeof(int)*nnt_loop1[nl]);
      loop2[nl]=(int*)malloc(sizeof(int)*nnt_loop2[nl]);
      
      lind=0;
      if(typ1==0 || typ1==1)
	for(inl=exloops[0];inl<=exloops[1];inl++){
	  my_link[inl]=nl;
	  my_loop[inl]=0;
	  loop1[nl][lind]=inl;lind++;
	}
      if(typ1==1)
	for(inl=exloops[2];inl<=exloops[3];inl++){
	  my_link[inl]=nl;
	  my_loop[inl]=0;
	  loop1[nl][lind]=inl;lind++;}
      /* if(ntl==4 || ntl==6) */
      /* 	for(inl=exloops[0]; */
      /* 	    inl<=exloops[1];inl++){ */
      /* 	  my_link[inl]=nl; */
      /* 	  loop1[nl][lind]=inl;lind++;} */
      /* if(ntl==8) */
      /* 	for(inl=exloops[2];inl<=exloops[3];inl++){ */
      /* 	  my_link[inl]=nl; */
      /* 	  loop1[nl][lind]=inl;lind++;} */
      /***/
      lind=0;
      if(typ2==0 || typ2==1)
	for(inl=exloops[start2];inl<=exloops[start2+1];inl++){
	  my_link[inl]=nl;
	  my_loop[inl]=1;
	  loop2[nl][lind]=inl;lind++;}
      if(typ2==1)
	for(inl=exloops[start2+2];inl<=exloops[start2+3];inl++){
	  my_link[inl]=nl;
	  my_loop[inl]=1;
	  loop2[nl][lind]=inl;lind++;}
      
      /* if(ntl==4 || ntl==6) */
      /* 	for(inl=exloops[2];inl<=exloops[3];inl++){ */
      /* 	  my_link[inl]=nl; */
      /* 	  loop2[nl][lind]=inl;lind++;} */
      
      /* if(ntl==6 || ntl==8) */
      /* 	for(inl=exloops[4];inl<=exloops[5];inl++){ */
      /* 	  my_link[inl]=nl; */
      /* 	  loop2[nl][lind]=inl;lind++;} */
      /* if(ntl==8) */
      /* 	for(inl=exloops[6];inl<=exloops[7];inl++){ */
      /* 	  my_link[inl]=nl; */
      /* 	  loop2[nl][lind]=inl;lind++;} */
      
      //if(ntl==4 || ntl==6){
      printf("LINK %d\n", nl);
      printf("Loop 1 (%d) : ",nnt_loop1[nl]);
      for(inl=0;inl<nnt_loop1[nl];inl++){
	printf("%d (%d) ",loop1[nl][inl],my_loop[loop1[nl][inl]]);
	fr_is_mobile[loop1[nl][inl]]=FR_MOB_FULL;
      }
      printf("\t");
      printf("Loop 2 (%d) : ", nnt_loop2[nl]);
      for(inl=0;inl<nnt_loop2[nl];inl++){
	printf("%d (%d) ",loop2[nl][inl],my_loop[loop2[nl][inl]]);
	fr_is_mobile[loop2[nl][inl]]=FR_MOB_FULL;
      }
      printf("\n");
      //}
      //if(ntl==4){
      //for(inl=0;inl<nnt_loop2[nl];inl++){
      //  fr_is_mobile[loop2[nl][inl]]=FR_MOB_FULL;
      //}
      //} 
    }
  }
  phpull_K=(double*)malloc(sizeof(double)*N_LINKS);
  for(i=0;i<N_LINKS;i++){ 
    phpull_K[i]=KTEMP; 
  } 
}

double calc_link_energy(int nt_c, double *rx, double *ry, double *rz){
  //here we calculate the centers of mass of the loops
  int i,d,nt1,nt2;
  double CM1[DIM],CM2[DIM];
  int thslp=my_link[nt_c];
  double energ=0,cmdistsq;
  Dlnksq=100;
  if(my_link[nt_c]>-1){
    //we calculate the centers of mass of both loops
    for(d=0;d<DIM;d++){
      CM1[d]=0;CM2[d]=0;
    }

    for(nt1=0;nt1<nnt_loop1[thslp];nt1++){
      if(loop1[thslp][nt1]!=nt_c){
	CM1[0]+=(rx[loop1[thslp][nt1]*N_PARTS_PER_NT+IPHO]+rx[loop1[thslp][nt1]*N_PARTS_PER_NT+ISUG]);
	CM1[1]+=(ry[loop1[thslp][nt1]*N_PARTS_PER_NT+IPHO]+ry[loop1[thslp][nt1]*N_PARTS_PER_NT+ISUG]);
	CM1[2]+=(rz[loop1[thslp][nt1]*N_PARTS_PER_NT+IPHO]+rz[loop1[thslp][nt1]*N_PARTS_PER_NT+ISUG]);
      }
      else{
	//printf("making the difference nt1!\n");
	CM1[0]+=(get_unf_coo_temp_x(nt_c*N_PARTS_PER_NT+IPHO)+get_unf_coo_temp_x(nt_c*N_PARTS_PER_NT+ISUG));
	CM1[1]+=(get_unf_coo_temp_y(nt_c*N_PARTS_PER_NT+IPHO)+get_unf_coo_temp_y(nt_c*N_PARTS_PER_NT+ISUG));
	CM1[2]+=(get_unf_coo_temp_z(nt_c*N_PARTS_PER_NT+IPHO)+get_unf_coo_temp_z(nt_c*N_PARTS_PER_NT+ISUG));
      }
    }
    for(nt2=0;nt2<nnt_loop2[thslp];nt2++){
      if(loop1[thslp][nt2]!=nt_c){
	CM2[0]+=(rx[loop2[thslp][nt2]*N_PARTS_PER_NT+IPHO]+rx[loop2[thslp][nt2]*N_PARTS_PER_NT+ISUG]);
	CM2[1]+=(ry[loop2[thslp][nt2]*N_PARTS_PER_NT+IPHO]+ry[loop2[thslp][nt2]*N_PARTS_PER_NT+ISUG]);
	CM2[2]+=(rz[loop2[thslp][nt2]*N_PARTS_PER_NT+IPHO]+rz[loop2[thslp][nt2]*N_PARTS_PER_NT+ISUG]);
      }
      else{
	//printf("making the difference nt2!\n");
	CM2[0]+=(get_unf_coo_temp_x(nt_c*N_PARTS_PER_NT+IPHO)+get_unf_coo_temp_x(nt_c*N_PARTS_PER_NT+ISUG));
	CM2[1]+=(get_unf_coo_temp_y(nt_c*N_PARTS_PER_NT+IPHO)+get_unf_coo_temp_y(nt_c*N_PARTS_PER_NT+ISUG));
	CM2[2]+=(get_unf_coo_temp_z(nt_c*N_PARTS_PER_NT+IPHO)+get_unf_coo_temp_z(nt_c*N_PARTS_PER_NT+ISUG));
      }
    }
    for(d=0;d<DIM;d++){
      
      CM1[d]/=(2.0*(double)nnt_loop1[thslp]);
      CM2[d]/=(2.0*(double)nnt_loop2[thslp]);
    }
    
    cmdistsq=SQ(CM1[0]-CM2[0])+SQ(CM1[1]-CM2[1])+SQ(CM1[2]-CM2[2]);
    energ=phpull_K[thslp]*exp(-cmdistsq/Dlnksq);
  }
  //if(phpull_K[nt_c]!=0){
  //printf("%d  %lf\n", nt_c, 0.5*phpull_K[nt_c]*(SQ(phpull_p0[nt_c][0]-posx)+SQ(phpull_p0[nt_c][1]-posy)+SQ(phpull_p0[nt_c][2])-posz));
  //turn 0.5*phpull_K[nt_c]*(SQ(phpull_p0[nt_c][0]-posx)+SQ(phpull_p0[nt_c][1]-posy)+SQ(phpull_p0[nt_c][2])-posz);
  //
  //printf("%lf %lf %lf   %lf %lf %lf    %lf %lf\n",CM1[0], CM1[1], CM1[2], CM2[0], CM2[1], CM2[2],sqrt(cmdistsq), energ);
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
      for(group=0;group<ERMSD_N_GROUPS;group++){
	l=getline(&lline, &st_l, ermsdfile);
	cpline=strndup(lline, MAXSTR);
	sscanf(cpline, "%s %s %s %lf", s1, s2, s3, &dtemp1);
	if(!(!strcmp(s1, "REMARK") && !strcmp(s2, "ERMSD") && !strcmp(s3, "GROUP"))){printf("Wrong syntax in ERMSD file ermsd_frags.lst\n");exit(ERR_INPUT);}
	ERMSD_PREF[group]=dtemp1;
	//printf("%lf   \n", ERMSD_PREF[group]);
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
      printf("ERMSD: ENERGY PREFACTORS\t");
      for(group=0;group<ERMSD_N_GROUPS;group++)
	printf("%d: %lf   ", group, ERMSD_PREF[group]);
      printf("\n");
      for(group=0;group<ERMSD_N_GROUPS;group++)
	ERMSD_PREF[group]*(0.5/((double)nnt_group[group]));
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
    fprintf(ermsd_obs, "# time internal_energy  ermsd_energy total_energy wall_energy\n");
  }
}

void MC_write_ermsd_obs(int step, double energ){
  double ermsd=get_ermsd();
  fprintf(ermsd_obs, "%d %lf %lf %lf %lf\n", step, energ-ERMSD_ENERG-WALL_ENERG,ermsd, energ, WALL_ENERG); 
}

double MC_wall_energy(double px, double py, double pz){
  double ret=0;
  if(wall_epsilon>0){
    double dist=fabs(wall_A*px+wall_B*py+wall_C*pz+wall_D)/wall_MODSQ;
    double d4=1.0/(SQ(SQ(dist)));
    //ret=wall_epsilon*dist+d4*d4*d4;
    ret=-wall_epsilon*exp(-dist/wall_sigma)/dist+d4*d4*d4;
    //printf("%lf  %lf\n", dist, ret);
  }
  return ret;
}
