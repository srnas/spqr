#include "mc_energies.h"

char *ew_gtl(FILE *source) {
  static size_t st_l=0;
  int n,i;
  static char *sch_s=NULL;
  while ((n=getline(&sch_s,&st_l,source))>0) {
//    printf("%d %d %s",n,strlen(s),s);
//    assert(s[n-1]=='\n');
    if (sch_s[n-1]=='\n')
      sch_s[n-1]='\0';
    if (isalpha(sch_s[0]))
      for(i=0;sch_s[i];++i)
        if (isspace(sch_s[i]))
          sch_s[i--]='\0';
    if (sch_s[0]!='#')
      break;
    }
  return n>0?sch_s:NULL;
}

void MC_initialize_energy_parameters(int mpi_id){
  int i,j,stf;
  
  for(i=0;i<N_BASES_SQ;i++){
    for(stf=0;stf<N_STFACES;stf++)
      b_st_well[i][stf]=0;
    nb_st_well[i]=0;
    for(j=0;j<WC_FACES_SQ;j++)
      nb_wc_well[i][j]=0;
   }
  for(i=0;i<N_BASES;i++)
    for(j=0;j<WC_FACES;j++)
      nb_bp_well[i][j]=0;
  //we initialize the 'special faces' of BPh interactions: guanine(W) and cytosine(H)
  for(i=0;i<N_BASES;i++)
    for(j=0;j<3;j++)
      nb_bp_spec_well[i][j]=0;
  
  mc_lj_sig=(double **)malloc(sizeof(double*)*N_MAX_TYPES);
  mc_lj_eps=(double **)malloc(sizeof(double*)*N_MAX_TYPES);
  
  mc_harm_k=(double *)malloc(sizeof(double)*MAX_BOND_TYPES);
  mc_harm_r=(double *)malloc(sizeof(double)*MAX_BOND_TYPES);
  
  mc_ang_k=(double *)malloc(sizeof(double)*MAX_ANG_TYPES);
  mc_ang_th=(double *)malloc(sizeof(double)*MAX_ANG_TYPES);
  
  for(i=0;i<N_MAX_TYPES;i++){
    mc_lj_sig[i]=(double *)malloc(sizeof(double)*N_MAX_TYPES);
    mc_lj_eps[i]=(double *)malloc(sizeof(double)*N_MAX_TYPES);
  }

  for(i=0;i<N_MAX_TYPES;i++){
    for(j=0;j<N_MAX_TYPES;j++){
      mc_lj_sig[i][j]=0;
      mc_lj_eps[i][j]=0;
    }
  }
  
  for(i=0;i<MAX_BOND_TYPES;i++){
    mc_harm_k[i]=0;
    mc_harm_r[i]=-1;
  }
  
  for(i=0;i<MAX_ANG_TYPES;i++){
    mc_ang_k[i]=0;
    mc_ang_th[i]=0;
  }
  
  /* defaults */
  /* mc_lj_sig[0][0]=DEFAULT_LJ_SIG; */
  /* mc_lj_eps[0][0]=DEFAULT_LJ_EPS; */
  /* mc_harm_k[0]=DEFAULT_HARM_K; */
  /* mc_harm_r[0]=DEFAULT_HARM_R; */
  /* mc_ang_k[0]=DEFAULT_ANG_K; */
  /* mc_ang_th[0]=DEFAULT_ANG_TH; */
  
  mc_chk_freq=DEFAULT_CHK_FREQ;
  
  /*   dih_k[0][0]=0; */
  /*   dih_phi[0][0]=0; */
  /*   dih_n[0][0]=0; */

  mc_r_cut=DEFAULT_MC_RCUT;
  mc_wc_rcut=DEFAULT_MC_RCUT;
  mc_bph_rcut=DEFAULT_MC_RCUT;
  mc_nb_rcut=DEFAULT_MC_RCUT;
  mc_n_types=N_BASES;
  mc_n_bond_types=DEFAULT_MC_N_BOND_TYPES;
  vl_skin=DEFAULT_VL_SKIN;

  MC_read_energy_parameters();
  MC_initialize_tabulated_energies(mpi_id);

#ifdef WARMUP
  int ig, ip, idd;
  ph_pintra_ave[TYP_ADENINE][GLYC_A][PUCK_2][0]=-5.62207;  ph_pintra_ave[TYP_ADENINE][GLYC_A][PUCK_2][1]=-4.84927;  ph_pintra_ave[TYP_ADENINE][GLYC_A][PUCK_2][2]= 0.322851;
  ph_pintra_ave[TYP_ADENINE][GLYC_H][PUCK_2][0]=-4.72923;  ph_pintra_ave[TYP_ADENINE][GLYC_H][PUCK_2][1]=-5.16856;  ph_pintra_ave[TYP_ADENINE][GLYC_H][PUCK_2][2]=-2.15611;
  ph_pintra_ave[TYP_ADENINE][GLYC_S][PUCK_2][0]=3.9117;  ph_pintra_ave[TYP_ADENINE][GLYC_S][PUCK_2][1]=-3.80641;  ph_pintra_ave[TYP_ADENINE][GLYC_S][PUCK_2][2]=1.08371;
  ph_pintra_ave[TYP_ADENINE][GLYC_A][PUCK_3][0]=-6.10965;  ph_pintra_ave[TYP_ADENINE][GLYC_A][PUCK_3][1]=-4.45512;  ph_pintra_ave[TYP_ADENINE][GLYC_A][PUCK_3][2]=0.776946;
  ph_pintra_ave[TYP_ADENINE][GLYC_H][PUCK_3][0]=-4.25369;  ph_pintra_ave[TYP_ADENINE][GLYC_H][PUCK_3][1]=-4.3032;  ph_pintra_ave[TYP_ADENINE][GLYC_H][PUCK_3][2]=-3.52153;
  ph_pintra_ave[TYP_ADENINE][GLYC_S][PUCK_3][0]=2.37143;  ph_pintra_ave[TYP_ADENINE][GLYC_S][PUCK_3][1]=-3.66344;  ph_pintra_ave[TYP_ADENINE][GLYC_S][PUCK_3][2]=4.03475;
  
  ph_pintra_ave[TYP_URACIL][GLYC_A][PUCK_2][0]=-2.66998;ph_pintra_ave[TYP_URACIL][GLYC_A][PUCK_2][1]= -5.03771;ph_pintra_ave[TYP_URACIL][GLYC_A][PUCK_2][2]= 0.298707;
  ph_pintra_ave[TYP_URACIL][GLYC_H][PUCK_2][0]=-2.2033;  ph_pintra_ave[TYP_URACIL][GLYC_H][PUCK_2][1]=-5.17071;  ph_pintra_ave[TYP_URACIL][GLYC_H][PUCK_2][2]=-1.47934;
  ph_pintra_ave[TYP_URACIL][GLYC_S][PUCK_2][0]=5.89192;  ph_pintra_ave[TYP_URACIL][GLYC_S][PUCK_2][1]=-1.36619;  ph_pintra_ave[TYP_URACIL][GLYC_S][PUCK_2][2]=1.03601;
  ph_pintra_ave[TYP_URACIL][GLYC_A][PUCK_3][0]=-3.36428;  ph_pintra_ave[TYP_URACIL][GLYC_A][PUCK_3][1]=-4.5741;  ph_pintra_ave[TYP_URACIL][GLYC_A][PUCK_3][2]=0.573644;
  ph_pintra_ave[TYP_URACIL][GLYC_H][PUCK_3][0]=-1.47508;  ph_pintra_ave[TYP_URACIL][GLYC_H][PUCK_3][1]=-4.3344;  ph_pintra_ave[TYP_URACIL][GLYC_H][PUCK_3][2]=-3.66775;
  ph_pintra_ave[TYP_URACIL][GLYC_S][PUCK_3][0]=4.21648;  ph_pintra_ave[TYP_URACIL][GLYC_S][PUCK_3][1]=-1.89827;  ph_pintra_ave[TYP_URACIL][GLYC_S][PUCK_3][2]=3.97488;
  
  for(ig=0;ig<3;ig++)
    for(ip=0;ip<2;ip++)
      for(idd=0;idd<DIM;idd++){
	ph_pintra_ave[TYP_CYTOSINE][ig][ip][idd]=ph_pintra_ave[TYP_URACIL][ig][ip][idd];
	ph_pintra_ave[TYP_GUANINE][ig][ip][idd]=ph_pintra_ave[TYP_ADENINE][ig][ip][idd];
      }
  
  ph_pinter_ave[TYP_ADENINE][GLYC_A][PUCK_2][0]=-0.0403659;  ph_pinter_ave[TYP_ADENINE][GLYC_A][PUCK_2][1]= -7.75968 ;  ph_pinter_ave[TYP_ADENINE][GLYC_A][PUCK_2][2]= 1.21347;
  ph_pinter_ave[TYP_ADENINE][GLYC_H][PUCK_2][0]=-0.347973;  ph_pinter_ave[TYP_ADENINE][GLYC_H][PUCK_2][1]= -7.75747 ;  ph_pinter_ave[TYP_ADENINE][GLYC_H][PUCK_2][2]= 1.23118;
  ph_pinter_ave[TYP_ADENINE][GLYC_S][PUCK_2][0]=-0.400171;  ph_pinter_ave[TYP_ADENINE][GLYC_S][PUCK_2][1]= -7.81304;  ph_pinter_ave[TYP_ADENINE][GLYC_S][PUCK_2][2]=  -1.27329;
  ph_pinter_ave[TYP_ADENINE][GLYC_A][PUCK_3][0]=-1.65899 ;  ph_pinter_ave[TYP_ADENINE][GLYC_A][PUCK_3][1]=-5.07092 ;  ph_pinter_ave[TYP_ADENINE][GLYC_A][PUCK_3][2]= 4.47217;
  ph_pinter_ave[TYP_ADENINE][GLYC_H][PUCK_3][0]=-4.01143;  ph_pinter_ave[TYP_ADENINE][GLYC_H][PUCK_3][1]= -6.87805;  ph_pinter_ave[TYP_ADENINE][GLYC_H][PUCK_3][2]=  1.03547;
  ph_pinter_ave[TYP_ADENINE][GLYC_S][PUCK_3][0]=3.13305;  ph_pinter_ave[TYP_ADENINE][GLYC_S][PUCK_3][1]= -5.18425;  ph_pinter_ave[TYP_ADENINE][GLYC_S][PUCK_3][2]=  -1.60444;
  
  ph_pinter_ave[TYP_URACIL][GLYC_A][PUCK_2][0]=3.7066;  ph_pinter_ave[TYP_URACIL][GLYC_A][PUCK_2][1]= -6.11415;  ph_pinter_ave[TYP_URACIL][GLYC_A][PUCK_2][2]= 1.10119;
  ph_pinter_ave[TYP_URACIL][GLYC_H][PUCK_2][0]=3.36573;  ph_pinter_ave[TYP_URACIL][GLYC_H][PUCK_2][1]= -6.35114;  ph_pinter_ave[TYP_URACIL][GLYC_H][PUCK_2][2]= 1.13665;
  ph_pinter_ave[TYP_URACIL][GLYC_S][PUCK_2][0]=2.06965;  ph_pinter_ave[TYP_URACIL][GLYC_S][PUCK_2][1]= -5.9423 ;  ph_pinter_ave[TYP_URACIL][GLYC_S][PUCK_2][2]=-0.977865;
  ph_pinter_ave[TYP_URACIL][GLYC_A][PUCK_3][0]=0.948652;  ph_pinter_ave[TYP_URACIL][GLYC_A][PUCK_3][1]= -4.26526 ;  ph_pinter_ave[TYP_URACIL][GLYC_A][PUCK_3][2]=4.36462;
  ph_pinter_ave[TYP_URACIL][GLYC_H][PUCK_3][0]=-0.425209;  ph_pinter_ave[TYP_URACIL][GLYC_H][PUCK_3][1]= -6.79188;  ph_pinter_ave[TYP_URACIL][GLYC_H][PUCK_3][2]= 1.1285;
  ph_pinter_ave[TYP_URACIL][GLYC_S][PUCK_3][0]=5.75244;  ph_pinter_ave[TYP_URACIL][GLYC_S][PUCK_3][1]= -2.74303;  ph_pinter_ave[TYP_URACIL][GLYC_S][PUCK_3][2]= -1.61584;
  
  for(ig=0;ig<3;ig++)
    for(ip=0;ip<2;ip++)
      for(idd=0;idd<DIM;idd++){
	ph_pinter_ave[TYP_CYTOSINE][ig][ip][idd]=ph_pinter_ave[TYP_URACIL][ig][ip][idd];
	ph_pinter_ave[TYP_GUANINE][ig][ip][idd]=ph_pinter_ave[TYP_ADENINE][ig][ip][idd];
      }
  ss_ang_ave=100.0*M_PI/180.0;
#endif
  
}

void MC_initialize_tabulated_energies(int mpi_id){
  char tablename[256];
  int i,j, typ_ind, ener, ntot, stf;
  int ntypsq=mc_n_types*mc_n_types;
  char ctemp1, ctemp2, ctemp3, ctemp4, ctemp5;
  //int typ_ind=mc_n_types*type1+type2;
  
  /* NON BONDED POTENTIAL WELLS */
  //printf("Reading table wells\n");
  sprintf(tablename, "tab_energs/table_wells.tab");
  FILE *table;
  if((table=fopen(tablename,"r"))==NULL){
    printf("file not found!\n");
    if(mpi_id==0)
      printf("No potential wells for non bonded interactions loaded.\n");
    exit(1);
  }
  else{
     fscanf(table,"%c%c", &ctemp1, &ctemp2);
    //printf("reading st %c  %c\n", ctemp1, ctemp2);
    if(ctemp1!='b' || ctemp2 != 'b'){
      printf("Wrong syntax at table_wells.tab file (bb)!\n%c %c\n", ctemp1, ctemp2);
      exit(1);}
    fscanf(table, "%lf", &(BB_PREF));
    fscanf(table, "%lf", &(BB_PREF_A));
    
    //intf("%lf read \n", BB_PREF);
    fscanf(table,"%c%c%c", &ctemp1, &ctemp2, &ctemp3);
    //intf("read %c %c %c\n", ctemp1, ctemp2, ctemp3);
    //printf("reading st %c  %c\n", ctemp1, ctemp2);
    if(ctemp2!='g' || ctemp3 != 'l'){
      printf("Wrong syntax at table_wells.tab file (gl)!\n%c %c\n", ctemp1, ctemp2);
      exit(1);}
    for(i=0;i<N_GLYC_STATES;i++)
      fscanf(table, "%lf", &(glyc_well[i]));
    
    /* NON BONDED POTENTIAL WELLS */
    fscanf(table,"%c%c%c%c", &ctemp1, &ctemp2, &ctemp3, &ctemp4);
    //printf("reading st %c  %c\n", ctemp1, ctemp2);
    if(ctemp2!='s' || ctemp3 != 't' || ctemp4!='b'){
      printf("Wrong syntax at table_wells.tab file (st , 1)!\n%c %c %c %c\n", ctemp1, ctemp2, ctemp3, ctemp4);
      exit(1);}
    for(stf=0;stf<N_STFACES;stf++){
      for(i=0;i<N_BASES;i++)
	for(j=0;j<N_BASES;j++){
	  fscanf(table, "%lf", &(b_st_well[N_BASES*i+j][stf]));
	}
    }
    /* for(i=0;i<N_BASES;i++) */
    /*   for(j=0;j<N_BASES;j++){ */
    /* 	if(b_st_well[N_BASES*i+j]!=b_st_well[N_BASES*j+i]){ */
    /* 	  printf("Bonded matrix (stacking) is not symmetric!\n"); */
    /* 	  exit(1); */
    /* 	} */
    /*   } */
    /* for(i=0;i<N_BASES;i++){ */
    /*   for(j=0;j<N_BASES;j++) printf("%lf ", b_st_well[N_BASES*j+i]); */
    /*   printf("\n");} */
    /*************************************/
    fscanf(table,"%c%c%c",&ctemp3, &ctemp1, &ctemp2);
    //printf("reading st %c  %c\n", ctemp1, ctemp2);
    if(ctemp1!='s' || ctemp2 != 't'){
      printf("Wrong syntax at table_wells.tab file (st , 2)!\n%c %c %c", ctemp1, ctemp2, ctemp3);
      exit(1);}
    for(i=0;i<N_BASES;i++)
      for(j=0;j<N_BASES;j++){
	fscanf(table, "%lf", &(nb_st_well[N_BASES*i+j]));
      }
    for(i=0;i<N_BASES;i++)
      for(j=0;j<N_BASES;j++){
	if(nb_st_well[N_BASES*i+j]!=nb_st_well[N_BASES*j+i]){
	  printf("Non-bonded matrix (stacking) is not symmetric!\n");
	  exit(1);
	}
      }
    /* for(i=0;i<N_BASES;i++){ */
    /*   for(j=0;j<N_BASES;j++) printf("%lf ", nb_st_well[N_BASES*j+i]); */
    /*   printf("\n");} */
    fscanf(table,"%c%c%c", &ctemp3,&ctemp1, &ctemp2);
    //printf("reading wc %c  %c\n", ctemp1, ctemp2);
    if(ctemp1!='w' || ctemp2 != 'c'){
      printf("Wrong syntax at table_wells.tab file (wc)!\n %c %c %c", ctemp1, ctemp2, ctemp3);
      exit(1);}
    for(i=0;i<N_BASES;i++)
      for(j=0;j<N_BASES;j++){
	//for(k=0;k<WC_FACES;k++)
	//for(l=k;l<WC_FACES;l++){
	fscanf(table,"%lf", &(nb_wc_well[N_BASES*i+j][0*WC_FACES+0]));
	fscanf(table,"%lf", &(nb_wc_well[N_BASES*i+j][1*WC_FACES+1]));
	fscanf(table,"%lf", &(nb_wc_well[N_BASES*i+j][2*WC_FACES+2]));
	
	fscanf(table,"%lf", &(nb_wc_well[N_BASES*i+j][1*WC_FACES+2]));
	nb_wc_well[N_BASES*j+i][2*WC_FACES+1]=nb_wc_well[N_BASES*i+j][1*WC_FACES+2];
	fscanf(table,"%lf", &(nb_wc_well[N_BASES*i+j][1*WC_FACES+0]));
	nb_wc_well[N_BASES*j+i][0*WC_FACES+1]=nb_wc_well[N_BASES*i+j][1*WC_FACES+0];
	fscanf(table,"%lf", &(nb_wc_well[N_BASES*i+j][2*WC_FACES+0]));
	nb_wc_well[N_BASES*j+i][0*WC_FACES+2]=nb_wc_well[N_BASES*i+j][2*WC_FACES+0];
	//}
      }
    
    //BASE-PAIRING IN PARALLEL CONFORMATION
    fscanf(table,"%c%c%c", &ctemp3,&ctemp1, &ctemp2);
    if(ctemp1!='w' || ctemp2 != 'P' ){
      printf("Wrong syntax at table_wells.tab file (wc)!\n %c %c %c", ctemp1, ctemp2, ctemp3);
      exit(1);}
    for(i=0;i<N_BASES;i++)
      for(j=0;j<N_BASES;j++){
	fscanf(table,"%lf", &(nb_wc_well_F[N_BASES*i+j][0*WC_FACES+0]));
	fscanf(table,"%lf", &(nb_wc_well_F[N_BASES*i+j][1*WC_FACES+1]));
	fscanf(table,"%lf", &(nb_wc_well_F[N_BASES*i+j][2*WC_FACES+2]));
	
	fscanf(table,"%lf", &(nb_wc_well_F[N_BASES*i+j][1*WC_FACES+2]));
	nb_wc_well_F[N_BASES*j+i][2*WC_FACES+1]=nb_wc_well_F[N_BASES*i+j][1*WC_FACES+2];
	fscanf(table,"%lf", &(nb_wc_well_F[N_BASES*i+j][1*WC_FACES+0]));
	nb_wc_well_F[N_BASES*j+i][0*WC_FACES+1]=nb_wc_well_F[N_BASES*i+j][1*WC_FACES+0];
	fscanf(table,"%lf", &(nb_wc_well_F[N_BASES*i+j][2*WC_FACES+0]));
	nb_wc_well_F[N_BASES*j+i][0*WC_FACES+2]=nb_wc_well_F[N_BASES*i+j][2*WC_FACES+0];
	//}
      }
    
    

    
    /* int mm, nn; */
    /* for(i=0;i<N_BASES;i++){   */
    /*   for(j=0;j<N_BASES;j++)  */
    /* 	for(mm=0;mm<WC_FACES;mm++) */
    /* 	  for(nn=0;nn<WC_FACES;nn++) */
	    
    /* 	    printf("%lf ", nb_wc_well[N_BASES*j+i][mm*WC_FACES+nn]);   */
    /*   printf("\n"); */
    /* }  */
    
    
    
    
    /************** BASE PHOSPHATE ***************/
    fscanf(table,"%c%c%c", &ctemp3,&ctemp1, &ctemp2);
    //printf("reading wc %c  %c\n", ctemp1, ctemp2);
    if(ctemp1!='b' || ctemp2 != 'p'){
      printf("Wrong syntax at table_wells.tab file (wc)!\n %c %c %c", ctemp1, ctemp2, ctemp3);
      exit(1);}
    for(i=0;i<N_BASES;i++)
      for(j=0;j<WC_FACES;j++)
	fscanf(table,"%lf", &(nb_bp_well[i][j]));
    /* for(i=0;i<N_BASES;i++)  */
    /*   for(j=0;j<WC_FACES;j++)  */
    /* 	printf("%lf   \n", nb_bp_well[i][j]);  */
    
    fscanf(table, "%lf", &(nb_bp_spec_well[2][0]));
    fscanf(table, "%lf", &(nb_bp_spec_well[2][1]));
    fscanf(table, "%lf", &(nb_bp_spec_well[2][2]));
    fscanf(table, "%lf", &(nb_bp_spec_well[3][0]));
    fscanf(table, "%lf", &(nb_bp_spec_well[3][1]));
    fscanf(table, "%lf", &(nb_bp_spec_well[3][2]));
    fclose(table);
  }
  
  
  /************************ BASE-PAIR SECOND DIHEDRALS **************************/
  
  sprintf(tablename, "tab_energs/table_wc_secdih.tab");
  if((table=fopen(tablename,"r"))==NULL){
    printf("file not found!\n");
    if(mpi_id==0)
      printf("No potential wells for non bonded interactions loaded.\n");
    exit(1);
  }
  else{
    fscanf(table,"%c%c%c%c", &ctemp1,&ctemp2, &ctemp3, &ctemp4);
    //printf("reading wc %c  %c\n", ctemp1, ctemp2);
    if(ctemp1!='A' || ctemp2 != 'm' || ctemp3 != 'i' || ctemp4!='n'){
      printf("Wrong syntax at table_wc_secdih.tab file (wc)!\n %c %c %c %c", ctemp1, ctemp2, ctemp3, ctemp4);
      exit(1);}
    for(i=0;i<N_BASES;i++)
      for(j=0;j<N_BASES;j++){
	fscanf(table,"%lf", &(wc_secdih_min[N_BASES*i+j][0*WC_FACES+0]));
	fscanf(table,"%lf", &(wc_secdih_min[N_BASES*i+j][1*WC_FACES+1]));
	fscanf(table,"%lf", &(wc_secdih_min[N_BASES*i+j][2*WC_FACES+2]));
	
	fscanf(table,"%lf", &(wc_secdih_min[N_BASES*i+j][1*WC_FACES+2]));
	wc_secdih_min[N_BASES*j+i][2*WC_FACES+1]=wc_secdih_min[N_BASES*i+j][1*WC_FACES+2];
	fscanf(table,"%lf", &(wc_secdih_min[N_BASES*i+j][1*WC_FACES+0]));
	wc_secdih_min[N_BASES*j+i][0*WC_FACES+1]=wc_secdih_min[N_BASES*i+j][1*WC_FACES+0];
	fscanf(table,"%lf", &(wc_secdih_min[N_BASES*i+j][2*WC_FACES+0]));
	wc_secdih_min[N_BASES*j+i][0*WC_FACES+2]=wc_secdih_min[N_BASES*i+j][2*WC_FACES+0];
	//}
      }
    
    fscanf(table,"%c%c%c%c%c", &ctemp5,&ctemp1, &ctemp2, &ctemp3, &ctemp4);
    //printf("reading wc %c  %c\n", ctemp1, ctemp2);
    if(ctemp1!='A' || ctemp2 != 'm' || ctemp3 != 'a' || ctemp4!='x'){
      printf("Wrong syntax at table_wells.tab file (wc)!\n %c %c %c", ctemp1, ctemp2, ctemp3);
      exit(1);}
    for(i=0;i<N_BASES;i++)
      for(j=0;j<N_BASES;j++){
	fscanf(table,"%lf", &(wc_secdih_max[N_BASES*i+j][0*WC_FACES+0]));
	fscanf(table,"%lf", &(wc_secdih_max[N_BASES*i+j][1*WC_FACES+1]));
	fscanf(table,"%lf", &(wc_secdih_max[N_BASES*i+j][2*WC_FACES+2]));
	
	fscanf(table,"%lf", &(wc_secdih_max[N_BASES*i+j][1*WC_FACES+2]));
	wc_secdih_max[N_BASES*j+i][2*WC_FACES+1]=wc_secdih_max[N_BASES*i+j][1*WC_FACES+2];
	fscanf(table,"%lf", &(wc_secdih_max[N_BASES*i+j][1*WC_FACES+0]));
	wc_secdih_max[N_BASES*j+i][0*WC_FACES+1]=wc_secdih_max[N_BASES*i+j][1*WC_FACES+0];
	fscanf(table,"%lf", &(wc_secdih_max[N_BASES*i+j][2*WC_FACES+0]));
	wc_secdih_max[N_BASES*j+i][0*WC_FACES+2]=wc_secdih_max[N_BASES*i+j][2*WC_FACES+0];
	//}
      }
    
    fscanf(table,"%c%c%c%c%c", &ctemp5,&ctemp1, &ctemp2, &ctemp3, &ctemp4);
    if(ctemp1!='P' || ctemp2 != 'm' || ctemp3 != 'i' || ctemp4!='n'){
      printf("Wrong syntax at table_wells.tab file (wc)!\n %c %c %c", ctemp1, ctemp2, ctemp3);
      exit(1);}
    for(i=0;i<N_BASES;i++)
      for(j=0;j<N_BASES;j++){
	fscanf(table,"%lf", &(wc_secdih_min_F[N_BASES*i+j][0*WC_FACES+0]));
	fscanf(table,"%lf", &(wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+1]));
	fscanf(table,"%lf", &(wc_secdih_min_F[N_BASES*i+j][2*WC_FACES+2]));
	
	fscanf(table,"%lf", &(wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+2]));
	wc_secdih_min_F[N_BASES*j+i][2*WC_FACES+1]=wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+2];
	fscanf(table,"%lf", &(wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+0]));
	wc_secdih_min_F[N_BASES*j+i][0*WC_FACES+1]=wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+0];
	fscanf(table,"%lf", &(wc_secdih_min_F[N_BASES*i+j][2*WC_FACES+0]));
	wc_secdih_min_F[N_BASES*j+i][0*WC_FACES+2]=wc_secdih_min_F[N_BASES*i+j][2*WC_FACES+0];
      }
    
    fscanf(table,"%c%c%c%c%c", &ctemp5,&ctemp1, &ctemp2, &ctemp3, &ctemp4);
    if(ctemp1!='P' || ctemp2 != 'm' || ctemp3 != 'a' || ctemp4!='x'){
      printf("Wrong syntax at table_wells.tab file (wc)!\n %c %c %c", ctemp1, ctemp2, ctemp3);
      exit(1);}
    for(i=0;i<N_BASES;i++)
      for(j=0;j<N_BASES;j++){
	fscanf(table,"%lf", &(wc_secdih_max_F[N_BASES*i+j][0*WC_FACES+0]));
	fscanf(table,"%lf", &(wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+1]));
	fscanf(table,"%lf", &(wc_secdih_max_F[N_BASES*i+j][2*WC_FACES+2]));
	
	fscanf(table,"%lf", &(wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+2]));
	wc_secdih_max_F[N_BASES*j+i][2*WC_FACES+1]=wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+2];
	fscanf(table,"%lf", &(wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+0]));
	wc_secdih_max_F[N_BASES*j+i][0*WC_FACES+1]=wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+0];
	fscanf(table,"%lf", &(wc_secdih_max_F[N_BASES*i+j][2*WC_FACES+0]));
	wc_secdih_max_F[N_BASES*j+i][0*WC_FACES+2]=wc_secdih_max_F[N_BASES*i+j][2*WC_FACES+0];
      }
    fclose(table);
  }














  
  
  /*** BACKBONE ***/
  /* SUGAR - PHOSPHATE - SUGAR */
  //we have to read the nine GLYCOSIDIC states : AA (0) , AH (1) , AS (2) , HA (3) , HH (4) , HS (5) , SA (6) , SH (7) , SS (8)
  //BUT WE USE THE PUCKERS
  //TYPE 0
  sprintf(tablename, "tab_energs/table_ssB1_p33.tab");
  //FILE *table;
  if((table=fopen(tablename, "r"))==NULL){
    if(mpi_id==0)
      printf("No backbone interaction between sugars.\n");
    table_ssB1_N_33=0;
  }
  else{
    if(mpi_id==0)
      printf("Reading backbone interaction: sugar-sugar (angle), pucks 3 3\n");
    fscanf(table,"%d", &(table_ssB1_N_33));
    fscanf(table,"%lf", &(table_ssB1_params_33[0]));
    fscanf(table,"%lf", &(table_ssB1_params_33[1]));
    table_ssB1_33=(double *)malloc(sizeof(double)*table_ssB1_N_33);
    for(ener=0;ener<table_ssB1_N_33;ener++)
      fscanf(table, "%lf", &(table_ssB1_33[ener]));
    fclose(table);
  }

  //TYPE 32
  sprintf(tablename, "tab_energs/table_ssB1_p32.tab");
  //FILE *table;
  if((table=fopen(tablename, "r"))==NULL){
    if(mpi_id==0)
      printf("No backbone interaction between sugars.\n");
    table_ssB1_N_32=0;
  }
  else{
    if(mpi_id==0)
      printf("Reading backbone interaction: sugar-sugar (angle), pucks 3 3\n");
    fscanf(table,"%d", &(table_ssB1_N_32));
    fscanf(table,"%lf", &(table_ssB1_params_32[0]));
    fscanf(table,"%lf", &(table_ssB1_params_32[1]));
    table_ssB1_32=(double *)malloc(sizeof(double)*table_ssB1_N_32);
    for(ener=0;ener<table_ssB1_N_32;ener++)
      fscanf(table, "%lf", &(table_ssB1_32[ener]));
    fclose(table);
  }
  
  //TYPE 23
  sprintf(tablename, "tab_energs/table_ssB1_p23.tab");
  //FILE *table;
  if((table=fopen(tablename, "r"))==NULL){
    if(mpi_id==0)
      printf("No backbone interaction between sugars.\n");
    table_ssB1_N_23=0;
  }
  else{
    if(mpi_id==0)
      printf("Reading backbone interaction: sugar-sugar (angle), pucks 3 3\n");
    fscanf(table,"%d", &(table_ssB1_N_23));
    fscanf(table,"%lf", &(table_ssB1_params_23[0]));
    fscanf(table,"%lf", &(table_ssB1_params_23[1]));
    table_ssB1_23=(double *)malloc(sizeof(double)*table_ssB1_N_23);
    for(ener=0;ener<table_ssB1_N_23;ener++)
      fscanf(table, "%lf", &(table_ssB1_23[ener]));
    fclose(table);
  }
//TYPE 22
  sprintf(tablename, "tab_energs/table_ssB1_p22.tab");
  //FILE *table;
  if((table=fopen(tablename, "r"))==NULL){
    if(mpi_id==0)
      printf("No backbone interaction between sugars.\n");
    table_ssB1_N_22=0;
  }
  else{
    if(mpi_id==0)
      printf("Reading backbone interaction: sugar-sugar (angle), pucks 3 3\n");
    fscanf(table,"%d", &(table_ssB1_N_22));
    fscanf(table,"%lf", &(table_ssB1_params_22[0]));
    fscanf(table,"%lf", &(table_ssB1_params_22[1]));
    table_ssB1_22=(double *)malloc(sizeof(double)*table_ssB1_N_22);
    for(ener=0;ener<table_ssB1_N_22;ener++)
      fscanf(table, "%lf", &(table_ssB1_22[ener]));
    fclose(table);
  }
  
 
  
  
  /* BASE - PHOSPHATE , INTRA NT */
  //glyc ANTI
  //puck3
  table_bpI_A3=(double **)malloc(sizeof(double *)*mc_n_types);
  for(i=0;i<N_BASES;i++){
    typ_ind=i;
    sprintf(tablename, "tab_energs/table_bpI_%d_gA3.tab", i);
    if((table=fopen(tablename, "r"))==NULL){
      if(mpi_id==0)
	printf("No backbone intra BP interaction for %d, glycosidic conformation ANTI.\n", i);
      table_bpI_N_A3[typ_ind][0]=-1;
      table_bpI_N_A3[typ_ind][1]=-1;
      table_bpI_N_A3[typ_ind][2]=-1;
    }
    else{
      if(mpi_id==0)
	printf("Reading backbone intra BP for %d, ANTI.\n", i);
      fscanf(table,"%d%d%d", &(table_bpI_N_A3[typ_ind][0]), &(table_bpI_N_A3[typ_ind][1]), &(table_bpI_N_A3[typ_ind][2])); //NX, NY, NZ
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_A3[typ_ind][0][0]), &(table_bpI_params_A3[typ_ind][1][0]), &(table_bpI_params_A3[typ_ind][2][0]));//DX, DY, DZ
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_A3[typ_ind][0][1]), &(table_bpI_params_A3[typ_ind][1][1]), &(table_bpI_params_A3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      ntot=table_bpI_N_A3[typ_ind][0]*table_bpI_N_A3[typ_ind][1]*table_bpI_N_A3[typ_ind][2];
      table_bpI_A3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	fscanf(table, "%lf", &(table_bpI_A3[typ_ind][ener]));
      }
      fclose(table);
    }
  }
  //puck2
  table_bpI_A2=(double **)malloc(sizeof(double *)*mc_n_types);
  for(i=0;i<N_BASES;i++){
    typ_ind=i;
    sprintf(tablename, "tab_energs/table_bpI_%d_gA2.tab", i);
    if((table=fopen(tablename, "r"))==NULL){
      if(mpi_id==0)
	printf("No backbone intra BP interaction for %d, glycosidic conformation ANTI.\n", i);
      table_bpI_N_A2[typ_ind][0]=-1;
      table_bpI_N_A2[typ_ind][1]=-1;
      table_bpI_N_A2[typ_ind][2]=-1;
    }
    else{
      if(mpi_id==0)
	printf("Reading backbone intra BP for %d, ANTI.\n", i);
      fscanf(table,"%d%d%d", &(table_bpI_N_A2[typ_ind][0]), &(table_bpI_N_A2[typ_ind][1]), &(table_bpI_N_A2[typ_ind][2])); //NX, NY, NZ
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_A2[typ_ind][0][0]), &(table_bpI_params_A2[typ_ind][1][0]), &(table_bpI_params_A2[typ_ind][2][0]));//DX, DY, DZ
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_A2[typ_ind][0][1]), &(table_bpI_params_A2[typ_ind][1][1]), &(table_bpI_params_A2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      ntot=table_bpI_N_A2[typ_ind][0]*table_bpI_N_A2[typ_ind][1]*table_bpI_N_A2[typ_ind][2];
      table_bpI_A2[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	fscanf(table, "%lf", &(table_bpI_A2[typ_ind][ener]));
      }
      fclose(table);
    }
  }
  
  //glyc HIGH ANTI
  table_bpI_H3=(double **)malloc(sizeof(double *)*mc_n_types);
  for(i=0;i<N_BASES;i++){
    typ_ind=i;
    sprintf(tablename, "tab_energs/table_bpI_%d_gH3.tab", i);
    if((table=fopen(tablename, "r"))==NULL){
      if(mpi_id==0)
	printf("No backbone intra SP interaction for %d, glycosidic conformation HIGH ANTI.\n", i);
      table_bpI_N_H3[typ_ind][0]=-1;
      table_bpI_N_H3[typ_ind][1]=-1;
      table_bpI_N_H3[typ_ind][2]=-1;
    }
    else{
      if(mpi_id==0)
	printf("Reading backbone intra SP for %d, HIGH ANTI.\n", i);
      fscanf(table,"%d%d%d", &(table_bpI_N_H3[typ_ind][0]), &(table_bpI_N_H3[typ_ind][1]), &(table_bpI_N_H3[typ_ind][2])); //NX, NY, NZ
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_H3[typ_ind][0][0]), &(table_bpI_params_H3[typ_ind][1][0]), &(table_bpI_params_H3[typ_ind][2][0]));//DX, DY, DZ
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_H3[typ_ind][0][1]), &(table_bpI_params_H3[typ_ind][1][1]), &(table_bpI_params_H3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      ntot=table_bpI_N_H3[typ_ind][0]*table_bpI_N_H3[typ_ind][1]*table_bpI_N_H3[typ_ind][2];
      table_bpI_H3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	fscanf(table, "%lf", &(table_bpI_H3[typ_ind][ener]));
      }
      fclose(table);
    }
  }
  table_bpI_H2=(double **)malloc(sizeof(double *)*mc_n_types);
  for(i=0;i<N_BASES;i++){
    typ_ind=i;
    sprintf(tablename, "tab_energs/table_bpI_%d_gH2.tab", i);
    if((table=fopen(tablename, "r"))==NULL){
      if(mpi_id==0)
	printf("No backbone intra SP interaction for %d, glycosidic conformation HIGH ANTI.\n", i);
      table_bpI_N_H2[typ_ind][0]=-1;
      table_bpI_N_H2[typ_ind][1]=-1;
      table_bpI_N_H2[typ_ind][2]=-1;
    }
    else{
      if(mpi_id==0)
	printf("Reading backbone intra SP for %d, HIGH ANTI.\n", i);
      fscanf(table,"%d%d%d", &(table_bpI_N_H2[typ_ind][0]), &(table_bpI_N_H2[typ_ind][1]), &(table_bpI_N_H2[typ_ind][2])); //NX, NY, NZ
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_H2[typ_ind][0][0]), &(table_bpI_params_H2[typ_ind][1][0]), &(table_bpI_params_H2[typ_ind][2][0]));//DX, DY, DZ
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_H2[typ_ind][0][1]), &(table_bpI_params_H2[typ_ind][1][1]), &(table_bpI_params_H2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      ntot=table_bpI_N_H2[typ_ind][0]*table_bpI_N_H2[typ_ind][1]*table_bpI_N_H2[typ_ind][2];
      table_bpI_H2[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	fscanf(table, "%lf", &(table_bpI_H2[typ_ind][ener]));
      }
      fclose(table);
    }
  }
  //glyc SYN
  table_bpI_S3=(double **)malloc(sizeof(double *)*mc_n_types);
  for(i=0;i<N_BASES;i++){
    typ_ind=i;
    sprintf(tablename, "tab_energs/table_bpI_%d_gS3.tab", i);
    if((table=fopen(tablename, "r"))==NULL){
      if(mpi_id==0)
	printf("No backbone intra SP interaction for %d, glycosidic conformation SYN.\n", i);
      table_bpI_N_S3[typ_ind][0]=-1;
      table_bpI_N_S3[typ_ind][1]=-1;
      table_bpI_N_S3[typ_ind][2]=-1;
    }
    else{
      if(mpi_id==0)
	printf("Reading backbone intra SP for %d, SYN.\n", i);
      fscanf(table,"%d%d%d", &(table_bpI_N_S3[typ_ind][0]), &(table_bpI_N_S3[typ_ind][1]), &(table_bpI_N_S3[typ_ind][2])); //NX, NY, NZ
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_S3[typ_ind][0][0]), &(table_bpI_params_S3[typ_ind][1][0]), &(table_bpI_params_S3[typ_ind][2][0]));//DX, DY, DZ
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_S3[typ_ind][0][1]), &(table_bpI_params_S3[typ_ind][1][1]), &(table_bpI_params_S3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      ntot=table_bpI_N_S3[typ_ind][0]*table_bpI_N_S3[typ_ind][1]*table_bpI_N_S3[typ_ind][2];
      table_bpI_S3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	fscanf(table, "%lf", &(table_bpI_S3[typ_ind][ener]));
      }
      fclose(table);
    }
  }
   table_bpI_S2=(double **)malloc(sizeof(double *)*mc_n_types);
  for(i=0;i<N_BASES;i++){
    typ_ind=i;
    sprintf(tablename, "tab_energs/table_bpI_%d_gS2.tab", i);
    if((table=fopen(tablename, "r"))==NULL){
      if(mpi_id==0)
	printf("No backbone intra SP interaction for %d, glycosidic conformation SYN.\n", i);
      table_bpI_N_S2[typ_ind][0]=-1;
      table_bpI_N_S2[typ_ind][1]=-1;
      table_bpI_N_S2[typ_ind][2]=-1;
    }
    else{
      if(mpi_id==0)
	printf("Reading backbone intra SP for %d, SYN.\n", i);
      fscanf(table,"%d%d%d", &(table_bpI_N_S2[typ_ind][0]), &(table_bpI_N_S2[typ_ind][1]), &(table_bpI_N_S2[typ_ind][2])); //NX, NY, NZ
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_S2[typ_ind][0][0]), &(table_bpI_params_S2[typ_ind][1][0]), &(table_bpI_params_S2[typ_ind][2][0]));//DX, DY, DZ
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_S2[typ_ind][0][1]), &(table_bpI_params_S2[typ_ind][1][1]), &(table_bpI_params_S2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      ntot=table_bpI_N_S2[typ_ind][0]*table_bpI_N_S2[typ_ind][1]*table_bpI_N_S2[typ_ind][2];
      table_bpI_S2[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	fscanf(table, "%lf", &(table_bpI_S2[typ_ind][ener]));
      }
      fclose(table);
    }
  }
  
  
  /* SUGAR - PHOSPHATE , INTER NT */
  //GLYC ANTI
  table_bpB_A3=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      sprintf(tablename, "tab_energs/table_bpB_%d%d_gA3.tab", i, j);
      //sprintf(tablename, "tab_energs/table_bpB_%d%d_gA.tab", j,i);
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No backbone inter SP interaction for %d and %d.\n", i,j);
	table_bpB_N_A3[typ_ind][0]=-1;
	table_bpB_N_A3[typ_ind][1]=-1;
	table_bpB_N_A3[typ_ind][2]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading backbone inter SP for %d and %d.\n", i,j);
	fscanf(table,"%d%d%d", &(table_bpB_N_A3[typ_ind][0]), &(table_bpB_N_A3[typ_ind][1]), &(table_bpB_N_A3[typ_ind][2])); //NX, NY, NZ
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_A3[typ_ind][0][0]), &(table_bpB_params_A3[typ_ind][1][0]), &(table_bpB_params_A3[typ_ind][2][0]));//DX, DY, DZ
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_A3[typ_ind][0][1]), &(table_bpB_params_A3[typ_ind][1][1]), &(table_bpB_params_A3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	ntot=table_bpB_N_A3[typ_ind][0]*table_bpB_N_A3[typ_ind][1]*table_bpB_N_A3[typ_ind][2];
	table_bpB_A3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_bpB_A3[typ_ind][ener]));
	}
	fclose(table);
      }
    }
  }
  table_bpB_A2=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      sprintf(tablename, "tab_energs/table_bpB_%d%d_gA2.tab", i, j);
      //sprintf(tablename, "tab_energs/table_bpB_%d%d_gA.tab", j,i);
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No backbone inter SP interaction for %d and %d.\n", i,j);
	table_bpB_N_A2[typ_ind][0]=-1;
	table_bpB_N_A2[typ_ind][1]=-1;
	table_bpB_N_A2[typ_ind][2]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading backbone inter SP for %d and %d.\n", i,j);
	fscanf(table,"%d%d%d", &(table_bpB_N_A2[typ_ind][0]), &(table_bpB_N_A2[typ_ind][1]), &(table_bpB_N_A2[typ_ind][2])); //NX, NY, NZ
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_A2[typ_ind][0][0]), &(table_bpB_params_A2[typ_ind][1][0]), &(table_bpB_params_A2[typ_ind][2][0]));//DX, DY, DZ
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_A2[typ_ind][0][1]), &(table_bpB_params_A2[typ_ind][1][1]), &(table_bpB_params_A2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	ntot=table_bpB_N_A2[typ_ind][0]*table_bpB_N_A2[typ_ind][1]*table_bpB_N_A2[typ_ind][2];
	table_bpB_A2[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_bpB_A2[typ_ind][ener]));
	}
	fclose(table);
      }
    }
  }
  //GLYC HIGH ANTI
  table_bpB_H3=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      sprintf(tablename, "tab_energs/table_bpB_%d%d_gH3.tab", i, j);
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No backbone inter SP interaction for %d and %d.\n", i,j);
	table_bpB_N_H3[typ_ind][0]=-1;
	table_bpB_N_H3[typ_ind][1]=-1;
	table_bpB_N_H3[typ_ind][2]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading backbone inter SP for %d and %d.\n", i,j);
	fscanf(table,"%d%d%d", &(table_bpB_N_H3[typ_ind][0]), &(table_bpB_N_H3[typ_ind][1]), &(table_bpB_N_H3[typ_ind][2])); //NX, NY, NZ
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_H3[typ_ind][0][0]), &(table_bpB_params_H3[typ_ind][1][0]), &(table_bpB_params_H3[typ_ind][2][0]));//DX, DY, DZ
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_H3[typ_ind][0][1]), &(table_bpB_params_H3[typ_ind][1][1]), &(table_bpB_params_H3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	ntot=table_bpB_N_H3[typ_ind][0]*table_bpB_N_H3[typ_ind][1]*table_bpB_N_H3[typ_ind][2];
	table_bpB_H3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_bpB_H3[typ_ind][ener]));
	}
	fclose(table);
      }
    }
  }
  table_bpB_H2=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      sprintf(tablename, "tab_energs/table_bpB_%d%d_gH2.tab", i, j);
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No backbone inter SP interaction for %d and %d.\n", i,j);
	table_bpB_N_H2[typ_ind][0]=-1;
	table_bpB_N_H2[typ_ind][1]=-1;
	table_bpB_N_H2[typ_ind][2]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading backbone inter SP for %d and %d.\n", i,j);
	fscanf(table,"%d%d%d", &(table_bpB_N_H2[typ_ind][0]), &(table_bpB_N_H2[typ_ind][1]), &(table_bpB_N_H2[typ_ind][2])); //NX, NY, NZ
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_H2[typ_ind][0][0]), &(table_bpB_params_H2[typ_ind][1][0]), &(table_bpB_params_H2[typ_ind][2][0]));//DX, DY, DZ
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_H2[typ_ind][0][1]), &(table_bpB_params_H2[typ_ind][1][1]), &(table_bpB_params_H2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	ntot=table_bpB_N_H2[typ_ind][0]*table_bpB_N_H2[typ_ind][1]*table_bpB_N_H2[typ_ind][2];
	table_bpB_H2[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_bpB_H2[typ_ind][ener]));
	}
	fclose(table);
      }
    }
  }
  
  //GLYC SYN
  table_bpB_S3=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      sprintf(tablename, "tab_energs/table_bpB_%d%d_gS3.tab", i, j);
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No backbone inter SP interaction for %d and %d.\n", i,j);
	table_bpB_N_S3[typ_ind][0]=-1;
	table_bpB_N_S3[typ_ind][1]=-1;
	table_bpB_N_S3[typ_ind][2]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading backbone inter SP for %d and %d.\n", i,j);
	fscanf(table,"%d%d%d", &(table_bpB_N_S3[typ_ind][0]), &(table_bpB_N_S3[typ_ind][1]), &(table_bpB_N_S3[typ_ind][2])); //NX, NY, NZ
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_S3[typ_ind][0][0]), &(table_bpB_params_S3[typ_ind][1][0]), &(table_bpB_params_S3[typ_ind][2][0]));//DX, DY, DZ
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_S3[typ_ind][0][1]), &(table_bpB_params_S3[typ_ind][1][1]), &(table_bpB_params_S3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	ntot=table_bpB_N_S3[typ_ind][0]*table_bpB_N_S3[typ_ind][1]*table_bpB_N_S3[typ_ind][2];
	table_bpB_S3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_bpB_S3[typ_ind][ener]));
	}
	fclose(table);
      }
    }
  }
  table_bpB_S2=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      sprintf(tablename, "tab_energs/table_bpB_%d%d_gS2.tab", i, j);
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No backbone inter SP interaction for %d and %d.\n", i,j);
	table_bpB_N_S2[typ_ind][0]=-1;
	table_bpB_N_S2[typ_ind][1]=-1;
	table_bpB_N_S2[typ_ind][2]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading backbone inter SP for %d and %d.\n", i,j);
	fscanf(table,"%d%d%d", &(table_bpB_N_S2[typ_ind][0]), &(table_bpB_N_S2[typ_ind][1]), &(table_bpB_N_S2[typ_ind][2])); //NX, NY, NZ
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_S2[typ_ind][0][0]), &(table_bpB_params_S2[typ_ind][1][0]), &(table_bpB_params_S2[typ_ind][2][0]));//DX, DY, DZ
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_S2[typ_ind][0][1]), &(table_bpB_params_S2[typ_ind][1][1]), &(table_bpB_params_S2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	ntot=table_bpB_N_S2[typ_ind][0]*table_bpB_N_S2[typ_ind][1]*table_bpB_N_S2[typ_ind][2];
	table_bpB_S2[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_bpB_S2[typ_ind][ener]));
	}
	fclose(table);
      }
    }
  }
  
  
  
  
  
  /* STACKING - BONDED */
  //s35
  table_nnB_0_s35=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      sprintf(tablename, "tab_energs/table_nnB_%d%d_0_s35.tab", i,j);
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No stacking interaction between %d and %d.\n", i,j);
	table_nnB_N_0_s35[typ_ind][0]=-1;
	table_nnB_N_0_s35[typ_ind][1]=-1;
	table_nnB_N_0_s35[typ_ind][2]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading stacking interaction between %d and %d.\n", i,j);
	fscanf(table,"%d%d%d", &(table_nnB_N_0_s35[typ_ind][0]), &(table_nnB_N_0_s35[typ_ind][1]), &(table_nnB_N_0_s35[typ_ind][2])); //NX, NY, NZ
	fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s35[typ_ind][0][0]), &(table_nnB_params_0_s35[typ_ind][1][0]), &(table_nnB_params_0_s35[typ_ind][2][0]));//DX, DY, DZ
	fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s35[typ_ind][0][1]), &(table_nnB_params_0_s35[typ_ind][1][1]), &(table_nnB_params_0_s35[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	
	ntot=table_nnB_N_0_s35[typ_ind][0]*table_nnB_N_0_s35[typ_ind][1]*table_nnB_N_0_s35[typ_ind][2];
	table_nnB_0_s35[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnB_0_s35[typ_ind][ener]));
	}
     	fclose(table);
      }
    }
  }
  //S53
  table_nnB_0_s53=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      sprintf(tablename, "tab_energs/table_nnB_%d%d_0_s53.tab", i,j);
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No stacking interaction between %d and %d.\n", i,j);
	table_nnB_N_0_s53[typ_ind][0]=-1;
	table_nnB_N_0_s53[typ_ind][1]=-1;
	table_nnB_N_0_s53[typ_ind][2]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading stacking interaction between %d and %d.\n", i,j);
	fscanf(table,"%d%d%d", &(table_nnB_N_0_s53[typ_ind][0]), &(table_nnB_N_0_s53[typ_ind][1]), &(table_nnB_N_0_s53[typ_ind][2])); //NX, NY, NZ
	fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s53[typ_ind][0][0]), &(table_nnB_params_0_s53[typ_ind][1][0]), &(table_nnB_params_0_s53[typ_ind][2][0]));//DX, DY, DZ
	fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s53[typ_ind][0][1]), &(table_nnB_params_0_s53[typ_ind][1][1]), &(table_nnB_params_0_s53[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	
	ntot=table_nnB_N_0_s53[typ_ind][0]*table_nnB_N_0_s53[typ_ind][1]*table_nnB_N_0_s53[typ_ind][2];
	table_nnB_0_s53[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnB_0_s53[typ_ind][ener]));
	}
     	fclose(table);
      }
    }
  }
  //S33
  table_nnB_0_s33=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      sprintf(tablename, "tab_energs/table_nnB_%d%d_0_s33.tab", i,j);
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No stacking interaction between %d and %d.\n", i,j);
	table_nnB_N_0_s33[typ_ind][0]=-1;
	table_nnB_N_0_s33[typ_ind][1]=-1;
	table_nnB_N_0_s33[typ_ind][2]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading stacking interaction between %d and %d.\n", i,j);
	fscanf(table,"%d%d%d", &(table_nnB_N_0_s33[typ_ind][0]), &(table_nnB_N_0_s33[typ_ind][1]), &(table_nnB_N_0_s33[typ_ind][2])); //NX, NY, NZ
	fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s33[typ_ind][0][0]), &(table_nnB_params_0_s33[typ_ind][1][0]), &(table_nnB_params_0_s33[typ_ind][2][0]));//DX, DY, DZ
	fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s33[typ_ind][0][1]), &(table_nnB_params_0_s33[typ_ind][1][1]), &(table_nnB_params_0_s33[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	
	ntot=table_nnB_N_0_s33[typ_ind][0]*table_nnB_N_0_s33[typ_ind][1]*table_nnB_N_0_s33[typ_ind][2];
	table_nnB_0_s33[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnB_0_s33[typ_ind][ener]));
	}
     	fclose(table);
      }
    }
  }
  //S55
  table_nnB_0_s55=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      sprintf(tablename, "tab_energs/table_nnB_%d%d_0_s55.tab", i,j);
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No stacking interaction between %d and %d.\n", i,j);
	table_nnB_N_0_s55[typ_ind][0]=-1;
	table_nnB_N_0_s55[typ_ind][1]=-1;
	table_nnB_N_0_s55[typ_ind][2]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading stacking interaction between %d and %d.\n", i,j);
	fscanf(table,"%d%d%d", &(table_nnB_N_0_s55[typ_ind][0]), &(table_nnB_N_0_s55[typ_ind][1]), &(table_nnB_N_0_s55[typ_ind][2])); //NX, NY, NZ
	fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s55[typ_ind][0][0]), &(table_nnB_params_0_s55[typ_ind][1][0]), &(table_nnB_params_0_s55[typ_ind][2][0]));//DX, DY, DZ
	fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s55[typ_ind][0][1]), &(table_nnB_params_0_s55[typ_ind][1][1]), &(table_nnB_params_0_s55[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	
	ntot=table_nnB_N_0_s55[typ_ind][0]*table_nnB_N_0_s55[typ_ind][1]*table_nnB_N_0_s55[typ_ind][2];
	table_nnB_0_s55[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnB_0_s55[typ_ind][ener]));
	}
     	fclose(table);
      }
    }
  }
  //////////////////////////////////
  
  table_nnB_1_s33=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //typ_ind2=N_BASES*j+i;
      sprintf(tablename, "tab_energs/table_nnB_%d%d_1_s33.tab", i,j);
      //FILE *table;
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No stacking-dihedral interaction between %d and %d.\n", i,j);
	table_nnB_N_1_s33[typ_ind]=-1;
	//table_nnB_N_1[typ_ind2]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading stacking-dihedral interaction between %d and %d.\n", i,j);
	fscanf(table,"%d", &(table_nnB_N_1_s33[typ_ind])); //N_ETA
	fscanf(table,"%lf", &(table_nnB_params_1_s33[typ_ind][0]));//D_ETA
	fscanf(table,"%lf", &(table_nnB_params_1_s33[typ_ind][1])); //ETAMIN
	ntot=table_nnB_N_1_s33[typ_ind];
	table_nnB_1_s33[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnB_1_s33[typ_ind][ener]));
	}
	fclose(table);
      }
    }
  }

 table_nnB_1_s35=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //typ_ind2=N_BASES*j+i;
      sprintf(tablename, "tab_energs/table_nnB_%d%d_1_s35.tab", i,j);
      //FILE *table;
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No stacking-dihedral interaction between %d and %d.\n", i,j);
	table_nnB_N_1_s35[typ_ind]=-1;
	//table_nnB_N_1[typ_ind2]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading stacking-dihedral interaction between %d and %d.\n", i,j);
	fscanf(table,"%d", &(table_nnB_N_1_s35[typ_ind])); //N_ETA
	fscanf(table,"%lf", &(table_nnB_params_1_s35[typ_ind][0]));//D_ETA
	fscanf(table,"%lf", &(table_nnB_params_1_s35[typ_ind][1])); //ETAMIN
	ntot=table_nnB_N_1_s35[typ_ind];
	table_nnB_1_s35[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnB_1_s35[typ_ind][ener]));
	}
	fclose(table);
      }
    }
  }

 table_nnB_1_s53=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //typ_ind2=N_BASES*j+i;
      sprintf(tablename, "tab_energs/table_nnB_%d%d_1_s53.tab", i,j);
      //FILE *table;
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No stacking-dihedral interaction between %d and %d.\n", i,j);
	table_nnB_N_1_s53[typ_ind]=-1;
	//table_nnB_N_1[typ_ind2]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading stacking-dihedral interaction between %d and %d.\n", i,j);
	fscanf(table,"%d", &(table_nnB_N_1_s53[typ_ind])); //N_ETA
	fscanf(table,"%lf", &(table_nnB_params_1_s53[typ_ind][0]));//D_ETA
	fscanf(table,"%lf", &(table_nnB_params_1_s53[typ_ind][1])); //ETAMIN
	ntot=table_nnB_N_1_s53[typ_ind];
	table_nnB_1_s53[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnB_1_s53[typ_ind][ener]));
	}
	fclose(table);
      }
    }
  }

 table_nnB_1_s55=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //typ_ind2=N_BASES*j+i;
      sprintf(tablename, "tab_energs/table_nnB_%d%d_1_s55.tab", i,j);
      //FILE *table;
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No stacking-dihedral interaction between %d and %d.\n", i,j);
	table_nnB_N_1_s55[typ_ind]=-1;
	//table_nnB_N_1[typ_ind2]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading stacking-dihedral interaction between %d and %d.\n", i,j);
	fscanf(table,"%d", &(table_nnB_N_1_s55[typ_ind])); //N_ETA
	fscanf(table,"%lf", &(table_nnB_params_1_s55[typ_ind][0]));//D_ETA
	fscanf(table,"%lf", &(table_nnB_params_1_s55[typ_ind][1])); //ETAMIN
	ntot=table_nnB_N_1_s55[typ_ind];
	table_nnB_1_s55[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnB_1_s55[typ_ind][ener]));
	}
	fclose(table);
      }
    }
  }
  

  
  /* STACKING -  NON-BONDED */
  //int nbtyp=1;
  table_nnN_0s3=(double **)malloc(sizeof(double *)*ntypsq);
  //typ_ind=0;
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      /*     //typ_ind2=N_BASES*j+i; */
      sprintf(tablename, "tab_energs/table_nnN_%d%d_0s3.tab", i, j);
      //FILE *table;
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No generic non-bonded stacking interaction.\n");
	table_nnN_N_0s3[typ_ind][0]=-1;
	table_nnN_N_0s3[typ_ind][1]=-1;
	table_nnN_N_0s3[typ_ind][2]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading generic non-bonded stacking interaction between.\n");
	fscanf(table,"%d%d%d", &(table_nnN_N_0s3[typ_ind][0]), &(table_nnN_N_0s3[typ_ind][1]), &(table_nnN_N_0s3[typ_ind][2])); //NX, NY, NZ
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_0s3[typ_ind][0][0]), &(table_nnN_params_0s3[typ_ind][1][0]), &(table_nnN_params_0s3[typ_ind][2][0]));//DX, DY, DZ
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_0s3[typ_ind][0][1]), &(table_nnN_params_0s3[typ_ind][1][1]), &(table_nnN_params_0s3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	ntot=table_nnN_N_0s3[typ_ind][0]*table_nnN_N_0s3[typ_ind][1]*table_nnN_N_0s3[typ_ind][2];
	table_nnN_0s3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnN_0s3[typ_ind][ener]));
	}
	fclose(table);
      }
    }
    
  }
  
  table_nnN_0s5=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //typ_ind2=N_BASES*j+i;
      //typ_ind=0;
      sprintf(tablename, "tab_energs/table_nnN_%d%d_0s5.tab", i, j);
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No generic non-bonded stacking_inv interaction.\n");
	table_nnN_N_0s5[typ_ind][0]=-1;
	table_nnN_N_0s5[typ_ind][1]=-1;
	table_nnN_N_0s5[typ_ind][2]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading generic non-bonded stacking_inv interaction.\n");
	fscanf(table,"%d%d%d", &(table_nnN_N_0s5[typ_ind][0]), &(table_nnN_N_0s5[typ_ind][1]), &(table_nnN_N_0s5[typ_ind][2])); //NX, NY, NZ
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_0s5[typ_ind][0][0]), &(table_nnN_params_0s5[typ_ind][1][0]), &(table_nnN_params_0s5[typ_ind][2][0]));//DX, DY, DZ
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_0s5[typ_ind][0][1]), &(table_nnN_params_0s5[typ_ind][1][1]), &(table_nnN_params_0s5[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	ntot=table_nnN_N_0s5[typ_ind][0]*table_nnN_N_0s5[typ_ind][1]*table_nnN_N_0s5[typ_ind][2];
	table_nnN_0s5[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnN_0s5[typ_ind][ener]));
	}
	fclose(table);
      }
    }
  }
  
  
  /* WATSON-CRICK */
  table_nnN_2=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
	//typ_ind2=N_BASES*j+i;
	sprintf(tablename, "tab_energs/table_nnN_%d%d_2.tab", i,j);
	//FILE *table;
	if((table=fopen(tablename, "r"))==NULL){
	  if(mpi_id==0)
	    printf("No watson-crick interaction between %d and %d.\n", i,j);
	  table_nnN_N_2[typ_ind][0]=-1;
	  table_nnN_N_2[typ_ind][1]=-1;
	  table_nnN_N_2[typ_ind][2]=-1;
	  table_nnN_N_2[typ_ind][3]=-1;
	  
	  //table_nnN_N_2[typ_ind2][0]=-1;
	  //table_nnN_N_2[typ_ind2][1]=-1;
	  //table_nnN_N_2[typ_ind2][2]=-1;
	  //table_nnN_2[typ_ind]=(double *)malloc(sizeof(double)*10);
	}
      else{
	if(mpi_id==0)
	  printf("Reading watson-crick interaction between %d and %d.\n", i,j);
	fscanf(table,"%d%d%d", &(table_nnN_N_2[typ_ind][0]), &(table_nnN_N_2[typ_ind][1]), &(table_nnN_N_2[typ_ind][2])); //NX, NY, NZ
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_2[typ_ind][0][0]), &(table_nnN_params_2[typ_ind][1][0]), &(table_nnN_params_2[typ_ind][2][0]));//DX, DY, DZ
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_2[typ_ind][0][1]), &(table_nnN_params_2[typ_ind][1][1]), &(table_nnN_params_2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	
	//intf("%d %d %d\n", table_nnN_N_2[typ_ind][0],table_nnN_N_2[typ_ind][1],table_nnN_N_2[typ_ind][2]);
	ntot=table_nnN_N_2[typ_ind][0]*table_nnN_N_2[typ_ind][1]*table_nnN_N_2[typ_ind][2];
	table_nnN_N_2[typ_ind][3]=ntot;
	//printf("%d  %d   %d   %d\n",N_BASES,ntypsq,typ_ind, ntot);
	
	table_nnN_2[typ_ind]=(double *)malloc(sizeof(double)*ntot*WC_FACES);
	
	for(ener=0;ener<ntot*WC_FACES;ener++){
	  fscanf(table, "%lf", &(table_nnN_2[typ_ind][ener]));
	  //printf("%d  %lf\n", ener, table_nnN_2[typ_ind][ener]);
	  //table_nnN_2[typ_ind2][ener]=table_nnN_2[typ_ind][ener];
	}
	fclose(table);
	//exit(1);
      }
    }
  }
  table_nnN_2_inv=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //typ_ind2=N_BASES*j+i;
      sprintf(tablename, "tab_energs/table_nnN_%d%d_2i.tab", i,j);
      //FILE *table;
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No watson-crick_inv interaction between %d and %d.\n", i,j);
	table_nnN_N_2_inv[typ_ind][0]=-1;
	table_nnN_N_2_inv[typ_ind][1]=-1;
	table_nnN_N_2_inv[typ_ind][2]=-1;
	table_nnN_N_2_inv[typ_ind][3]=-1;
	//table_nnN_N_2_inv[typ_ind2][0]=-1;
	//table_nnN_N_2_inv[typ_ind2][1]=-1;
	//table_nnN_N_2_inv[typ_ind2][2]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading watson-crick_inv interaction between %d and %d.\n", i,j);
	fscanf(table,"%d%d%d", &(table_nnN_N_2_inv[typ_ind][0]), &(table_nnN_N_2_inv[typ_ind][1]), &(table_nnN_N_2_inv[typ_ind][2])); //NX, NY, NZ
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_2_inv[typ_ind][0][0]), &(table_nnN_params_2_inv[typ_ind][1][0]), &(table_nnN_params_2_inv[typ_ind][2][0]));//DX, DY, DZ
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_2_inv[typ_ind][0][1]), &(table_nnN_params_2_inv[typ_ind][1][1]), &(table_nnN_params_2_inv[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	ntot=table_nnN_N_2_inv[typ_ind][0]*table_nnN_N_2_inv[typ_ind][1]*table_nnN_N_2_inv[typ_ind][2];
	table_nnN_N_2_inv[typ_ind][3]=ntot;

	table_nnN_2_inv[typ_ind]=(double *)malloc(sizeof(double)*ntot*WC_FACES);
	
	for(ener=0;ener<ntot*WC_FACES;ener++){
	  fscanf(table, "%lf", &(table_nnN_2_inv[typ_ind][ener]));
	  //table_nnN_2_inv[typ_ind2][ener]=table_nnN_2_inv[typ_ind][ener];
	}
	fclose(table);
      }
    }
  }
  
  //printf("n max types : %d\n", N_MAX_TYPES);
  table_nnN_3=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //typ_ind2=N_BASES*j+i;
      sprintf(tablename, "tab_energs/table_nnN_%d%d_3.tab", i,j);
      //FILE *table;
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No watson-crick_dihedral interaction between %d and %d.\n", i,j);
	table_nnN_N_3[typ_ind]=-1;
	//table_nnN_N_3[typ_ind2]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading watson-crick_dihedral interaction between %d and %d.\n", i,j);
	fscanf(table,"%d", &(table_nnN_N_3[typ_ind])); //N_THETA
	fscanf(table,"%lf", &(table_nnN_params_3[typ_ind][0]));//D_THETA
	fscanf(table,"%lf", &(table_nnN_params_3[typ_ind][1])); //THETAMIN
	ntot=table_nnN_N_3[typ_ind];
	table_nnN_3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	//table_nnN_3[typ_ind2]=(double *)malloc(sizeof(double)*ntot);
	
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnN_3[typ_ind][ener]));
	  //table_nnN_3[typ_ind2][ener]=table_nnN_3[typ_ind][ener];
	}
	fclose(table);
      }
    }
  }
  
  /* for the case of parallel orientations */
  
  table_nnN_2_F=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //typ_ind2=N_BASES*j+i;
      sprintf(tablename, "tab_energs/table_nnN_%d%d_2_F.tab", i,j);
      //FILE *table;
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No watson-crick interaction between %d and %d in parallel conformation.\n", i,j);
	table_nnN_N_2_F[typ_ind][0]=-1;
	table_nnN_N_2_F[typ_ind][1]=-1;
	table_nnN_N_2_F[typ_ind][2]=-1;
	table_nnN_N_2_F[typ_ind][3]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading watson-crick interaction between %d and %d in parallel conformation.\n", i,j);
	fscanf(table,"%d%d%d", &(table_nnN_N_2_F[typ_ind][0]), &(table_nnN_N_2_F[typ_ind][1]), &(table_nnN_N_2_F[typ_ind][2])); //NX, NY, NZ
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_2_F[typ_ind][0][0]), &(table_nnN_params_2_F[typ_ind][1][0]), &(table_nnN_params_2_F[typ_ind][2][0]));//DX, DY, DZ
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_2_F[typ_ind][0][1]), &(table_nnN_params_2_F[typ_ind][1][1]), &(table_nnN_params_2_F[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	ntot=table_nnN_N_2_F[typ_ind][0]*table_nnN_N_2_F[typ_ind][1]*table_nnN_N_2_F[typ_ind][2];
	table_nnN_N_2_F[typ_ind][3]=ntot;
	table_nnN_2_F[typ_ind]=(double *)malloc(sizeof(double)*ntot*WC_FACES);
	for(ener=0;ener<ntot*WC_FACES;ener++){
	  fscanf(table, "%lf", &(table_nnN_2_F[typ_ind][ener]));
	}
	fclose(table);
      }
    }
  }




  table_nnN_2_inv_F=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //typ_ind2=N_BASES*j+i;
      sprintf(tablename, "tab_energs/table_nnN_%d%d_2i_F.tab", i,j);
      //FILE *table;
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No watson-crick_inv interaction between %d and %d in parallel conformation.\n", i,j);
	table_nnN_N_2_inv_F[typ_ind][0]=-1;
	table_nnN_N_2_inv_F[typ_ind][1]=-1;
	table_nnN_N_2_inv_F[typ_ind][2]=-1;
	table_nnN_N_2_inv_F[typ_ind][3]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading watson-crick_inv interaction between %d and %d in parallel conformation.\n", i,j);
	fscanf(table,"%d%d%d", &(table_nnN_N_2_inv_F[typ_ind][0]), &(table_nnN_N_2_inv_F[typ_ind][1]), &(table_nnN_N_2_inv_F[typ_ind][2])); //NX, NY, NZ
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_2_inv_F[typ_ind][0][0]), &(table_nnN_params_2_inv_F[typ_ind][1][0]), &(table_nnN_params_2_inv_F[typ_ind][2][0]));//DX, DY, DZ
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_2_inv_F[typ_ind][0][1]), &(table_nnN_params_2_inv_F[typ_ind][1][1]), &(table_nnN_params_2_inv_F[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	ntot=table_nnN_N_2_inv_F[typ_ind][0]*table_nnN_N_2_inv_F[typ_ind][1]*table_nnN_N_2_inv_F[typ_ind][2];
	table_nnN_N_2_inv_F[typ_ind][3]=ntot;
	table_nnN_2_inv_F[typ_ind]=(double *)malloc(sizeof(double)*ntot*WC_FACES);
	for(ener=0;ener<ntot*WC_FACES;ener++){
	  fscanf(table, "%lf", &(table_nnN_2_inv_F[typ_ind][ener]));
	}
	fclose(table);
      }
    }
  }
  

  







  table_nnN_3_F=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      sprintf(tablename, "tab_energs/table_nnN_%d%d_3_F.tab", i,j);
      if((table=fopen(tablename, "r"))==NULL){
	if(mpi_id==0)
	  printf("No watson-crick_dihedral interaction between %d and %d in parallel configuration.\n", i,j);
	table_nnN_N_3_F[typ_ind]=-1;
      }
      else{
	if(mpi_id==0)
	  printf("Reading watson-crick_dihedral interaction between %d and %d in parallel conformation.\n", i,j);
	fscanf(table,"%d", &(table_nnN_N_3_F[typ_ind])); //N_THETA
	fscanf(table,"%lf", &(table_nnN_params_3_F[typ_ind][0]));//D_THETA
	fscanf(table,"%lf", &(table_nnN_params_3_F[typ_ind][1])); //THETAMIN
	ntot=table_nnN_N_3_F[typ_ind];
	table_nnN_3_F[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnN_3_F[typ_ind][ener]));
	}
	fclose(table);
      }
    }
  }
  






  // BASE-PHOSPHATE
  table_npN_0=(double **)malloc(sizeof(double *)*N_BASES);
  for(i=0;i<N_BASES;i++){
    sprintf(tablename, "tab_energs/table_npN_%d_0.tab", i);
    if((table=fopen(tablename, "r"))==NULL){
      if(mpi_id==0)
	printf("No base-phosphate interaction for type %d.\n", i);
      table_npN_N_0[i][0]=-1;
      table_npN_N_0[i][1]=-1;
      table_npN_N_0[i][2]=-1;
    }
    else{
      if(mpi_id==0)
	printf("Reading base-phosphate interaction for %d.\n", i);
      fscanf(table,"%d%d%d", &(table_npN_N_0[i][0]), &(table_npN_N_0[i][1]), &(table_npN_N_0[i][2])); //NX, NY, NZ
      fscanf(table,"%lf%lf%lf", &(table_npN_params_0[i][0][0]), &(table_npN_params_0[i][1][0]), &(table_npN_params_0[i][2][0]));//DX, DY, DZ
      fscanf(table,"%lf%lf%lf", &(table_npN_params_0[i][0][1]), &(table_npN_params_0[i][1][1]), &(table_npN_params_0[i][2][1])); //XMIN, YMIN, ZMIN
      
      ntot=table_npN_N_0[i][0]*table_npN_N_0[i][1]*table_npN_N_0[i][2];
      
      //intf("type %d\t\t%d %d %d     %d\n", i, table_npN_N_0[i][0],table_npN_N_0[i][1],table_npN_N_0[i][2], ntot);
      table_npN_0[i]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	fscanf(table, "%lf", &(table_npN_0[i][ener]));
      }
      fclose(table);
    }
  }
}

void MC_read_energy_parameters(){
  /* using Westphal's code from mpc_io.c */
  int i, j;
  FILE *energy_params;
  char *s;
  int nb_flag=0, b_flag=0;
  
  energy_params=fopen("energies.pms", "r");
  if(energy_params==NULL){
    printf("No energy parameter file given. Using default values for energies,\n");
  } else {
    //printf("Reading energy parameters from file.\n");
    while (s=ew_gtl(energy_params)) {
      if(!strcmp(s, "MC_RCUT")) {
	assert(s=ew_gtl(energy_params));
	//assert(sscanf(s, "%lf", &mc_r_cut)==1);
	assert(sscanf(s, "%lf%lf%lf%lf", &mc_r_cut, &mc_wc_rcut, &mc_bph_rcut, &mc_nb_rcut)==4);
      }
      /* if(!strcmp(s, "MC_NTYPES")) { */
      /* 	assert(s=ew_gtl(energy_params)); */
      /* 	assert(sscanf(s, "%d\n", &mc_n_types)==1); */
      /* 	nb_flag=1; */
      /* 	if(mc_n_types>N_MAX_TYPES){ */
      /* 	  fprintf(stderr, "MC_NTYPES can not be greater than %d by default. Change value of N_MAX_TYPES in mc_energies.h.\n", N_MAX_TYPES); */
      /* 	  exit(ERR_INPUT); */
      /* 	} */
      /* } */
      if(!strcmp(s, "MC_N_BOND_TYPES")) {
	assert(s=ew_gtl(energy_params));
#ifdef DIHEDRALS
	assert(sscanf(s, "%d%d%d\n", &mc_n_bond_types, &mc_n_ang_types,&mc_n_dih_types)==3);
#else
	assert(sscanf(s, "%d%d\n", &mc_n_bond_types, &mc_n_ang_types)==2);
#endif
	b_flag=1;
      }
      if(!strcmp(s, "VL_SKIN")) {
	assert(s=ew_gtl(energy_params));
	assert(sscanf(s, "%lf\n", &vl_skin)==1);
      }
      if(!strcmp(s, "MC_BETWEEN_CHECKPOINTS")) {
	assert(s=ew_gtl(energy_params));
	assert(sscanf(s, "%d\n", &mc_chk_freq)==1);
      }
      /* if(!strcmp(s, "LJ_EPSILON")) { */
      /* 	//assert(s=ew_gtl(energy_params)); */
      /* 	if(nb_flag==0){ */
      /* 	  fprintf(stderr,"MC_NTYPES not set properly. Define it BEFORE LJ_EPSILON and LJ_SIGMA in energys.prm.\n"); */
      /* 	  exit(ERR_INPUT); */
      /* 	} */
      /* 	for(i=0;i<mc_n_types;i++){ */
      /* 	  if(!fscanf(energy_params, "%lf ", &(mc_lj_eps[i][i]))){ */
      /* 	    fprintf(stderr, "Wrong LJ_EPSILON input at energies.pms.\n"); */
      /* 	    exit(ERR_INPUT); */
      /* 	  } */
      /* 	} */
      /* 	for(i=0;i<mc_n_types;i++) */
      /* 	  for(j=i+1;j<mc_n_types;j++){ */
      /* 	    if(!fscanf(energy_params, "%lf ", &(mc_lj_eps[i][j]))){ */
      /* 	    fprintf(stderr, "Wrong LJ_EPSILON input at energies.pms.\n"); */
      /* 	    exit(ERR_INPUT); */
      /* 	    } */
      /* 	    mc_lj_eps[j][i]=mc_lj_eps[i][j]; */
      /* 	  } */
      /* 	//assert(sscanf(s, "\n")==1); */
      /* } */
      /* if(!strcmp(s, "LJ_SIGMA")) { */
      /* 	//assert(s=ew_gtl(energy_params)); */
      /* 	if(nb_flag==0){ */
      /* 	  fprintf(stderr, "MC_NTYPES not set properly. Define it BEFORE LJ_EPSILON and LJ_SIGMA in energies.pms.\n"); */
      /* 	  exit(ERR_INPUT); */
      /* 	} */
      /* 	for(i=0;i<mc_n_types;i++){ */
      /* 	  if(!fscanf(energy_params, "%lf ", &(mc_lj_sig[i][i]))){ */
      /* 	    fprintf(stderr, "Wrong LJ_SIGMA input at energies.pms.\n"); */
      /* 	    exit(ERR_INPUT); */
      /* 	  } */
      /* 	}  */
      /* 	for(i=0;i<mc_n_types;i++) */
      /* 	  for(j=i+1;j<mc_n_types;j++){ */
      /* 	    if(!fscanf(energy_params, "%lf ", &(mc_lj_sig[i][j]))){ */
      /* 	    fprintf(stderr, "Wrong LJ_SIGMA input at energies.pms.\n"); */
      /* 	    exit(ERR_INPUT); */
      /* 	    } */
      /* 	    mc_lj_sig[j][i]=mc_lj_sig[i][j]; */
      /* 	  } */
      /* } */
  
      /* if(!strcmp(s, "HARM_K")) { */
      /* 	if(b_flag==0){ */
      /* 	  fprintf(stderr, "MC_N_BOND_NTYPES not set properly. Define it BEFORE HARM_K and HARM_R in energies.pms.\n"); */
      /* 	  exit(ERR_INPUT); */
      /* 	} */
      /* 	for(i=0;i<mc_n_bond_types;i++){ */
      /* 	  //for(i=0;i<mc_n_types;i++){ */
      /* 	  if(!fscanf(energy_params, "%lf ", &(mc_harm_k[i]))){ */
      /* 	    fprintf(stderr, "Wrong HARM_K input at energies.pms.\n"); */
      /* 	    exit(ERR_INPUT); */
      /* 	  } */
      /* 	} */
      /* 	//assert(sscanf(s, "%lf ", &(mc_harm_k[i]))==1); */
      /* 	//assert(sscanf(s, "\n")==1); */
      /* } */
      /* if(!strcmp(s, "HARM_R")) { */
      /* 	//assert(s=ew_gtl(energy_params)); */
      /* 	if(b_flag==0){ */
      /* 	  fprintf(stderr,"MC_N_BOND_NTYPES not set properly. Define it BEFORE HARM_K and HARM_R in energies.pms.\n"); */
      /* 	  exit(ERR_INPUT); */
      /* 	} */
      /* 	for(i=0;i<mc_n_bond_types;i++) */
      /* 	  if(!fscanf(energy_params, "%lf ", &(mc_harm_r[i]))){ */
      /* 	    fprintf(stderr, "Wrong HARM_R input at energies.pms.\n"); */
      /* 	    exit(ERR_INPUT); */
      /* 	  } */
      /* 	//assert(sscanf(s, "%lf ", &(mc_harm_r[i]))==1); */
      /* 	//assert(sscanf(s, "\n")==1); */
      /* } */
      if(!strcmp(s, "ANG_K")) {
	//assert(s=ew_gtl(energy_params));
	if(b_flag==0){
	  fprintf(stderr,"MC_N_BOND_NTYPES not set properly. Define it BEFORE ANG_K in energies.pms.\n");
	  exit(ERR_INPUT);
	}
	for(i=0;i<mc_n_ang_types;i++)
	  if(!fscanf(energy_params, "%lf ", &(mc_ang_k[i]))){
	    fprintf(stderr, "Wrong ANG_K input at energies.pms.\n");
	    exit(ERR_INPUT);
	  }
      }
      if(!strcmp(s, "ANG_TH")) {
	//assert(s=ew_gtl(energy_params));
	if(b_flag==0){
	  fprintf(stderr,"MC_N_BOND_NTYPES not set properly. Define it BEFORE ANG_TH in energies.pms.\n");
	  exit(ERR_INPUT);
	}
	for(i=0;i<mc_n_ang_types;i++)
	  if(!fscanf(energy_params, "%lf ", &(mc_ang_th[i]))){
	    fprintf(stderr, "Wrong ANG_TH input at energies.pms.\n");
	    exit(ERR_INPUT);
	  }
      }
      
#ifdef DIHEDRALS
      if(!strcmp(s, "DIH_K")) {
	//assert(s=ew_gtl(energy_params));
	if(b_flag==0){
	  fprintf(stderr,"MC_N_BOND_NTYPES not set properly. Define it BEFORE DIH_K in energies.pms.\n");
	  exit(ERR_INPUT);
	}
	for(i=0;i<mc_n_dih_types;i++)
	  if(!fscanf(energy_params, "%lf ", &(mc_dih_k[i]))){
	    fprintf(stderr, "Wrong DIH_K input at energies.pms.\n");
	    exit(ERR_INPUT);
	  }
      }
      if(!strcmp(s, "DIH_PHI")) {
	//assert(s=ew_gtl(energy_params));
	if(b_flag==0){
	  fprintf(stderr,"MC_N_BOND_NTYPES not set properly. Define it BEFORE DIH_PHI in energies.pms.\n");
	  exit(ERR_INPUT);
	}
	for(i=0;i<mc_n_dih_types;i++)
	  if(!fscanf(energy_params, "%lf ", &(mc_dih_phi[i]))){
	    fprintf(stderr, "Wrong DIH_PHI input at energies.pms.\n");
	    exit(ERR_INPUT);
	  }
      }
      if(!strcmp(s, "DIH_N")) {
	//assert(s=ew_gtl(energy_params));
	if(b_flag==0){
	  fprintf(stderr,"MC_N_BOND_NTYPES not set properly. Define it BEFORE DIH_N in energies.pms.\n");
	  exit(ERR_INPUT);
	}
	for(i=0;i<mc_n_dih_types;i++)
	  if(!fscanf(energy_params, "%lf ", &(mc_dih_n[i]))){
	    fprintf(stderr, "Wrong DIH_N input at energies.pms.\n");
	    exit(ERR_INPUT);
	  }
      }
#endif
    }
    fclose(energy_params);
  }
}

double MC_calc_intra_energy(int nt_c, int *flag){
  double d_vec[DIM], r_vec[DIM], self_e, glyc_e;
  int at_c=N_PARTS_PER_NT*nt_c;
  
  
  d_vec[0]=get_unf_coo_temp_x(at_c+IPHO)  -  get_unf_coo_temp_x(at_c+IBAS);
  d_vec[1]=get_unf_coo_temp_y(at_c+IPHO)  -  get_unf_coo_temp_y(at_c+IBAS);
  d_vec[2]=get_unf_coo_temp_z(at_c+IPHO)  -  get_unf_coo_temp_z(at_c+IBAS);
  //d_vec[2]=get_unf_coo_temp_z(at_c+IPHO)  -  get_unf_coo_temp_z(at_c+ISUG);
  proj_on_nt(d_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
  self_e=MC_calc_BP_intra(mc_types[at_c], r_vec, mc_temp_glyc[nt_c], mc_temp_puck[nt_c], flag);
#ifdef WRMVERB
  if(self_e>=0)printf("# P %d\n# N %d\n", nt_c, nt_c);
#endif
#ifdef XYZDEBUG
  if(*flag==1){
    printf("CLASH: intra %d\t%lf %lf %lf\n", nt_c, r_vec[0], r_vec[1], r_vec[2]);
  }
#endif
  
  glyc_e = -glyc_well[mc_temp_glyc[nt_c]];
#ifdef NEW_BIA
  if(mc_target_temp>1.0)
    self_e*=(mc_target_temp);
#endif
  return self_e + glyc_e;
}



double MC_calc_bonded_energy(int nt_c, double *rx, double *ry, double *rz, int *mc_flag){
  /* we receive the nt number, so we look for its atom index */
  double mc_ev_glob_rcut_sq=EV_GLOB_RCUT*EV_GLOB_RCUT;
  int at_c=nt_c*N_PARTS_PER_NT;
  int b, at_b; 
 //double e_bond;
  int temp_flag=0;
  double dist, d_vec[DIM], d_vec1[DIM], d_vec2[DIM];
  double bonded_energy=0.0;
  double temp_ev;
  int nt_neigh;
  double cos_ang_SS;
  //double epsilon, sigma
  double bond_ev=0.0, bond_bp=0.0;
  int typ_ind, typ_ind2;
  int puckind;
  double r_vec[DIM], r_vec_inv[DIM], t_vec[DIM];
  double eta, eta2;
  double bond_st=0.0, bond_ss=0.0;
  double b_ev_temp1, b_ev_temp2;
  double ori_o[DIM], ori_n[DIM];
  //double st_theta_0=180.0/180.0*M_PI, st_theta;
  double backbone_e=0;
  int flag_ss=0;
  double phpos[DIM];
  int this_typ;
      
  *mc_flag=0;
  /* stretching (radial) potential */
  
  for(b=0;b<mc_nbonds[nt_c][0];b++){
    nt_neigh=mc_bondlist[nt_c][b];
    //first, we see if they are stacked!
    at_b=nt_neigh*N_PARTS_PER_NT;
    typ_ind = N_BASES*mc_types[at_c] + mc_types[at_b];
    typ_ind2= N_BASES*mc_types[at_b] + mc_types[at_c];
    d_vec[0]=get_unf_coo_x(rx, at_b) - get_unf_coo_temp_x(at_c);
    d_vec[1]=get_unf_coo_y(ry, at_b) - get_unf_coo_temp_y(at_c);
    d_vec[2]=get_unf_coo_z(rz, at_b) - get_unf_coo_temp_z(at_c);
    //(rx[ba]+mc_pbox[ba][0]*box_l[0]) - (mc_temp_x[0]+mc_temp_pbox[0][0]*box_l[0]);
    
    dist=sqrt(d_vec[0]*d_vec[0]+d_vec[1]*d_vec[1]+d_vec[2]*d_vec[2]);
    proj_on_nt(d_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
    proj_on_nt_inv(d_vec, rx, ry,rz, nt_neigh, r_vec_inv);
    eta=0;//calc_st_psdihedr(mc_temp_x, mc_temp_y, mc_temp_z, rx,ry,rz,nt_c,nt_neigh); 
    
    temp_flag=0;
    eta2=calc_wc_psdihedr(mc_temp_x, mc_temp_y, mc_temp_z, rx,ry,rz,nt_c,nt_neigh); 
    /* STACKING DEPENDS ON THE GLYC CONFORMATION - AA IS BONDED, ANYTHING ELSE IS NONBONDED */ 
    if(fabs(eta2)<STANG || fabs(eta2) > M_PI-STANG)
      bond_st=MC_calc_nnB_stacking(typ_ind, typ_ind2, r_vec, r_vec_inv, eta, &temp_flag);
    else {bond_st=0;temp_flag=1;}
    
    
    /* else { */
    /*   //printf("%d %d other stacking case\n", nt_c, nt_neigh); */
    /*   bond_st=MC_calc_nnN_stacking(typ_ind,typ_ind2, r_vec, r_vec_inv, &temp_flag); // NO PAIR SPECIFIC!! */
    /* } */

    
    if(temp_flag==0){
      bonded_energy+=bond_st;
    }
  
    
    /** BASE - BASE **/
    t_vec[0]=get_unf_coo_x(rx, at_b) - get_unf_coo_temp_x(at_c);
    t_vec[1]=get_unf_coo_y(ry, at_b) - get_unf_coo_temp_y(at_c);
    t_vec[2]=get_unf_coo_z(rz, at_b) - get_unf_coo_temp_z(at_c);
    if(t_vec[0]*t_vec[0]+t_vec[1]*t_vec[1]+t_vec[2]*t_vec[2] <mc_ev_glob_rcut_sq){
      proj_on_nt(t_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
      proj_on_nt_inv(t_vec, rx,ry,rz, nt_neigh, r_vec_inv);
      b_ev_temp1=energy_hardcore(r_vec,     EV_BASE1, EV_BASE2, EV_BASE3,0,0,0);
      b_ev_temp2=energy_hardcore(r_vec_inv, EV_BASE1, EV_BASE2, EV_BASE3,0,0,0);
#ifdef WRMVERB
      if(b_ev_temp1>0 || b_ev_temp2>0) printf("# N %d\n# N %d\n", nt_c, nt_neigh);
#endif
#ifndef WARMUP
      if(b_ev_temp1>MC_CAP || b_ev_temp2>MC_CAP) {*mc_flag=2;
#ifdef XYZDEBUG
	printf("NB CLASH (base-base): %d  %d\t%lf %lf %lf\n", nt_c, nt_neigh, r_vec[0], r_vec[1], r_vec[2]);
#else
	return bonded_energy;
#endif
      }
#endif
      //printf("%lf %lf %lf %lf\n", r_vec[0], r_vec[1], r_vec[2], sqrt(r_vec[0]*r_vec[0]+r_vec[1]*r_vec[1]+r_vec[2]*r_vec[2]))
      bond_ev+=(b_ev_temp1+b_ev_temp2);
    }
    /** PHOSPHATE-SUGAR **/
    t_vec[0]=get_unf_coo_x(rx, at_b+IPHO) - get_unf_coo_temp_x(at_c+ISUG);
    t_vec[1]=get_unf_coo_y(ry, at_b+IPHO) - get_unf_coo_temp_y(at_c+ISUG);
    t_vec[2]=get_unf_coo_z(rz, at_b+IPHO) - get_unf_coo_temp_z(at_c+ISUG);
    if(t_vec[0]*t_vec[0]+t_vec[1]*t_vec[1]+t_vec[2]*t_vec[2] <mc_ev_glob_rcut_sq){
      proj_on_nt(t_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
      proj_on_nt_inv(t_vec,rx,ry,rz, nt_neigh, r_vec_inv);
      
      if(nt_neigh-nt_c==1){//connected through direct bond - bonded parameters
	b_ev_temp1=energy_hardcore(r_vec    , EV_PHOSUGR, EV_PHOSUGR, EV_PHOSUGR, 0,0,0);
	b_ev_temp2=energy_hardcore(r_vec_inv, EV_SUGPHOR, EV_SUGPHOR, EV_SUGPHOR, 0,0,0);
      }
      else{
	b_ev_temp1=energy_hardcore(r_vec    , EV_PHOSUG1, EV_PHOSUG2, EV_PHOSUG3, EV_PHOSUG4, EV_PHOSUG5, EV_PHOSUG6); // phosp wrt sugar
	b_ev_temp2=energy_hardcore(r_vec_inv, EV_SUGPHO1, EV_SUGPHO2, EV_SUGPHO3, EV_SUGPHO4, EV_SUGPHO5, EV_SUGPHO6); // base wrt pho
      }
#ifdef WRMVERB
      if(b_ev_temp1>0 || b_ev_temp2>0) printf("# P %d\n# S %d\n", nt_neigh, nt_c);
#endif

#ifndef WARMUP
      if(b_ev_temp1>MC_CAP || b_ev_temp2>MC_CAP) {*mc_flag=14;
#ifdef XYZDEBUG
	//printf("B CLASH (phosp-sugar): %d  %d\n", nt_neigh, nt_c);
	printf("B CLASH (phosp-sugar  %lf %lf): %d  %d\n\t%lf %lf %lf   (%lf %lf %lf %lf %lf %lf)\n", b_ev_temp1, b_ev_temp2, nt_neigh, nt_c, r_vec[0], r_vec[1], r_vec[2],  EV_PHOSUG1, EV_PHOSUG2, EV_PHOSUG3, EV_PHOSUG4, EV_PHOSUG5, EV_PHOSUG6);
	printf("\t %lf %lf %lf   (%lf %lf %lf)\n", r_vec_inv[0], r_vec_inv[1], r_vec_inv[2],EV_SUGPHO1, EV_SUGPHO2, EV_SUGPHO3);
#endif
      }
#endif
      bond_ev+=(b_ev_temp1+b_ev_temp2);
    }
    /** SUGAR-PHOSPHATE **/
    t_vec[0]=get_unf_coo_x(rx, at_b+ISUG) - get_unf_coo_temp_x(at_c+IPHO);
    t_vec[1]=get_unf_coo_y(ry, at_b+ISUG) - get_unf_coo_temp_y(at_c+IPHO);
    t_vec[2]=get_unf_coo_z(rz, at_b+ISUG) - get_unf_coo_temp_z(at_c+IPHO);
    if(t_vec[0]*t_vec[0]+t_vec[1]*t_vec[1]+t_vec[2]*t_vec[2] <mc_ev_glob_rcut_sq){
      proj_on_nt(t_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
      proj_on_nt_inv(t_vec,rx,ry,rz, nt_neigh, r_vec_inv);
      if(nt_c-nt_neigh==1){
	b_ev_temp1=energy_hardcore(r_vec    , EV_SUGPHOR, EV_SUGPHOR, EV_SUGPHOR, 0,0,0); // base  wrt phosp
	b_ev_temp2=energy_hardcore(r_vec_inv, EV_PHOSUGR, EV_PHOSUGR, EV_PHOSUGR, 0,0,0); // phosp wrt base
      }
      else {
	b_ev_temp1=energy_hardcore(r_vec    , EV_SUGPHO1, EV_SUGPHO2, EV_SUGPHO3, EV_SUGPHO4, EV_SUGPHO5, EV_SUGPHO6); // base  wrt phosp
	b_ev_temp2=energy_hardcore(r_vec_inv, EV_PHOSUG1, EV_PHOSUG2, EV_PHOSUG3, EV_PHOSUG4, EV_PHOSUG5, EV_PHOSUG6); // phosp wrt base
      }
#ifdef WRMVERB
      if(b_ev_temp1>0 || b_ev_temp2>0) printf("# P %d\n# S %d\n", nt_c, nt_neigh);
#endif

#ifndef WARMUP
      if(b_ev_temp1>MC_CAP || b_ev_temp2>MC_CAP) {*mc_flag=15;
#ifdef XYZDEBUG
	printf("B CLASH (sugar-phosp  %lf %lf): %d  %d\n\t%lf %lf %lf   (%lf %lf %lf)\n", b_ev_temp1, b_ev_temp2, nt_c,nt_neigh, r_vec[0], r_vec[1], r_vec[2], EV_SUGPHO1, EV_SUGPHO2, EV_SUGPHO3);
	printf("\t %lf %lf %lf   (%lf %lf %lf %lf %lf %lf)\n", r_vec_inv[0], r_vec_inv[1], r_vec_inv[2], EV_PHOSUG1, EV_PHOSUG2, EV_PHOSUG3, EV_PHOSUG4, EV_PHOSUG5, EV_PHOSUG6);
#endif
      }
#endif
      bond_ev+=(b_ev_temp1+b_ev_temp2);
    }
    
    /** SUGAR-BASE **/
    t_vec[0]=get_unf_coo_x(rx, at_b+ISUG) - get_unf_coo_temp_x(at_c+IBAS);
    t_vec[1]=get_unf_coo_y(ry, at_b+ISUG) - get_unf_coo_temp_y(at_c+IBAS);
    t_vec[2]=get_unf_coo_z(rz, at_b+ISUG) - get_unf_coo_temp_z(at_c+IBAS);
    //calc_min_vec(rx[N_PARTS_PER_NT*ba+ISUG], ry[N_PARTS_PER_NT*ba+ISUG], rz[N_PARTS_PER_NT*nt_neigh+ISUG], mc_temp_x[0], mc_temp_y[0], mc_temp_z[0], t_vec, &r);
    if(t_vec[0]*t_vec[0]+t_vec[1]*t_vec[1]+t_vec[2]*t_vec[2] <mc_ev_glob_rcut_sq){
      proj_on_nt(t_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
      proj_on_nt_inv(t_vec,rx,ry,rz,nt_neigh, r_vec_inv);
      
      if(is_purine(nt_c)){
	b_ev_temp1=energy_hardcore(r_vec    , EV_SUGPUR1, EV_SUGPUR2, EV_SUGPUR3,  EV_SUGPUR4, EV_SUGPUR5, EV_SUGPUR6 ); // sugar wrt base
	b_ev_temp2=energy_hardcore(r_vec_inv, EV_PURSUG1, EV_PURSUG2, EV_PURSUG3,  EV_PURSUG4, EV_PURSUG5, EV_PURSUG6 ); //EV_BASRSUG4, EV_BASRSUG5, EV_BASRSUG6 ); // base  wrt sugar
      }
      else{//f(is_pyrimidine(nt_c)){
	b_ev_temp1=energy_hardcore(r_vec    , EV_SUGPYR1, EV_SUGPYR2, EV_SUGPYR3,  EV_SUGPYR4, EV_SUGPYR5, EV_SUGPYR6 ); // sugar wrt base
	b_ev_temp2=energy_hardcore(r_vec_inv, EV_PYRSUG1, EV_PYRSUG2, EV_PYRSUG3,  EV_PYRSUG4, EV_PYRSUG5, EV_PYRSUG6 ); //EV_BASRSUG4, EV_BASRSUG5, EV_BASRSUG6 ); // base  wrt sugar
      }
#ifdef WRMVERB
      if(b_ev_temp1>0 || b_ev_temp2>0) printf("# N %d\n# N %d\n", nt_c, nt_neigh);
#endif
#ifndef WARMUP
      if(b_ev_temp1>MC_CAP || b_ev_temp2>MC_CAP) {*mc_flag=3; 
#ifdef XYZDEBUG
	printf("NB CLASH (sugar-base): %d  %d\n", nt_c, nt_neigh);
#else
	return bonded_energy;
#endif 
      }
#endif
      
      
      bond_ev+=(b_ev_temp1+b_ev_temp2);
    }
    /** BASE-SUGAR **/    
    t_vec[0]=get_unf_coo_x(rx, at_b) - get_unf_coo_temp_x(at_c+ISUG);
    t_vec[1]=get_unf_coo_y(ry, at_b) - get_unf_coo_temp_y(at_c+ISUG);
    t_vec[2]=get_unf_coo_z(rz, at_b) - get_unf_coo_temp_z(at_c+ISUG);
    if(t_vec[0]*t_vec[0]+t_vec[1]*t_vec[1]+t_vec[2]*t_vec[2] <mc_ev_glob_rcut_sq){
      proj_on_nt(t_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
      proj_on_nt_inv(t_vec,rx,ry,rz, nt_neigh, r_vec_inv);
      if(is_purine(nt_neigh)){
	b_ev_temp1=energy_hardcore(r_vec    , EV_PURSUG1, EV_PURSUG2, EV_PURSUG3, EV_PURSUG4, EV_PURSUG5, EV_PURSUG6); //EV_PURRSUG4, EV_PURRSUG5, EV_PURRSUG6); // base  wrt sugar
	b_ev_temp2=energy_hardcore(r_vec_inv, EV_SUGPUR1, EV_SUGPUR2, EV_SUGPUR3, EV_SUGPUR4, EV_SUGPUR5, EV_SUGPUR6); // sugar wrt base
      }
      else{
	//  }    else if(is_pyrimidine(nt_neigh)){
	b_ev_temp1=energy_hardcore(r_vec    , EV_PYRSUG1, EV_PYRSUG2, EV_PYRSUG3, EV_PYRSUG4, EV_PYRSUG5, EV_PYRSUG6); //EV_PYRRSUG4, EV_PYRRSUG5, EV_PYRRSUG6); // base  wrt sugar
	b_ev_temp2=energy_hardcore(r_vec_inv, EV_SUGPYR1, EV_SUGPYR2, EV_SUGPYR3, EV_SUGPYR4, EV_SUGPYR5, EV_SUGPYR6); // sugar wrt base
      }
#ifdef WRMVERB
      if(b_ev_temp1>0 || b_ev_temp2>0) printf("# N %d\n# N %d\n", nt_c, nt_neigh);
#endif
#ifndef WARMUP
      if(b_ev_temp1>MC_CAP || b_ev_temp2>MC_CAP) {*mc_flag=4;
#ifdef XYZDEBUG
	printf("NB CLASH (base-sugar): %d  %d\n", nt_c, nt_neigh);
#else
	return bonded_energy;
#endif 
      }
#endif
      bond_ev+=(b_ev_temp1+b_ev_temp2);
    }
    //sugar - sugar case
    /* d_vec[0]=get_unf_coo_x(rx, at_b+ISUG) - get_unf_coo_temp_x(at_c+ISUG); */
    /* d_vec[1]=get_unf_coo_y(ry, at_b+ISUG) - get_unf_coo_temp_y(at_c+ISUG); */
    /* d_vec[2]=get_unf_coo_z(rz, at_b+ISUG) - get_unf_coo_temp_z(at_c+ISUG); */
    /* dist=sqrt(d_vec[0]*d_vec[0]+d_vec[1]*d_vec[1]+d_vec[2]*d_vec[2]); */
    /* bond_ss=calc_ssB_sugars(dist); */
    /* if(bond_ss>MC_CAP) {*mc_flag=5;// printf("%d %d     %lf\n", nt_c,nt_neigh, dist); */
    /*   return bonded_energy;} */
    //SUGAR-SUGAR ANGLE
    if(nt_c - nt_neigh==1){ // we use trial P
      d_vec1[0]=get_unf_coo_temp_x(at_c+ISUG) - get_unf_coo_temp_x(at_c+IPHO);
      d_vec1[1]=get_unf_coo_temp_y(at_c+ISUG) - get_unf_coo_temp_y(at_c+IPHO);
      d_vec1[2]=get_unf_coo_temp_z(at_c+ISUG) - get_unf_coo_temp_z(at_c+IPHO);
      d_vec2[0]=get_unf_coo_x(rx, at_b+ISUG)  - get_unf_coo_temp_x(at_c+IPHO);
      d_vec2[1]=get_unf_coo_y(ry, at_b+ISUG)  - get_unf_coo_temp_y(at_c+IPHO);
      d_vec2[2]=get_unf_coo_z(rz, at_b+ISUG)  - get_unf_coo_temp_z(at_c+IPHO);
      //printf("%d %d :   %lf %lf %lf       %lf %lf %lf         %lf %lf %lf\n", nt_c, nt_neigh, get_unf_coo_temp_x(at_c+ISUG),get_unf_coo_temp_y(at_c+ISUG),get_unf_coo_temp_z(at_c+ISUG),get_unf_coo_temp_x(at_c+IPHO),get_unf_coo_temp_y(at_c+IPHO),get_unf_coo_temp_z(at_c+IPHO),get_unf_coo_x(rx, at_b+ISUG),get_unf_coo_y(rx, at_b+ISUG),get_unf_coo_z(rx, at_b+ISUG));
    }
    else if (nt_c - nt_neigh == -1){
      d_vec1[0]=get_unf_coo_temp_x(at_c+ISUG) - get_unf_coo_x(rx, at_b+IPHO);
      d_vec1[1]=get_unf_coo_temp_y(at_c+ISUG) - get_unf_coo_y(ry, at_b+IPHO);
      d_vec1[2]=get_unf_coo_temp_z(at_c+ISUG) - get_unf_coo_z(rz, at_b+IPHO);
      d_vec2[0]=get_unf_coo_x(rx, at_b+ISUG)  - get_unf_coo_x(rx, at_b+IPHO);
      d_vec2[1]=get_unf_coo_y(ry, at_b+ISUG)  - get_unf_coo_y(ry, at_b+IPHO);
      d_vec2[2]=get_unf_coo_z(rz, at_b+ISUG)  - get_unf_coo_z(rz, at_b+IPHO);
    }
    cos_ang_SS=(d_vec1[0]*d_vec2[0]+d_vec1[1]*d_vec2[1]+d_vec1[2]*d_vec2[2])/sqrt((d_vec1[0]*d_vec1[0]+d_vec1[1]*d_vec1[1]+d_vec1[2]*d_vec1[2])*(d_vec2[0]*d_vec2[0]+d_vec2[1]*d_vec2[1]+d_vec2[2]*d_vec2[2]));
    
    //dist=sqrt(d_vec[0]*d_vec[0]+d_vec[1]*d_vec[1]+d_vec[2]*d_vec[2]);
    if(nt_c<nt_neigh)
      puckind = MC_get_pucks(mc_temp_puck[nt_c], mc_puck[nt_neigh]);
    else
      puckind = MC_get_pucks(mc_puck[nt_neigh],  mc_temp_puck[nt_c]);
    flag_ss=0;
    bond_ss=calc_ssB1_sugars(acos(cos_ang_SS),puckind);
    if(bond_ss>=0) flag_ss=1;
#ifdef WRMVERB
    if(bond_ss>=0) printf("# N %d\n# N %d\n", nt_c, nt_neigh);
#endif
    
#ifdef WARMUP
    if(flag_ss==1){
      flag_ss=0;
      //bond_ss=MC_CAP-1;
      bond_ss=MC_CAP*(acos(cos_ang_SS)-ss_ang_ave)*(acos(cos_ang_SS)-ss_ang_ave);
    }
#endif
    if(flag_ss==1) {*mc_flag=5;// printf("%d %d     %lf\n", nt_c,nt_neigh, dist);

#ifdef XYZDEBUG
      printf("BB CLASH (sugar angle): %d  %d\t %lf   %lf\n", nt_c, nt_neigh, acos(cos_ang_SS),  acos(cos_ang_SS)*180/M_PI);
#else
      return bonded_energy;
#endif 
    }
    
    /** PHOSPHATE-BASE - only from i to i+1 **/
    if(nt_c<nt_neigh){//nt_up=nt_c;nt_do=nt_neigh;
      // EV between P of nt_c and B of nt_neigh
      // tab between P of nt_neigh and B of nt_c
      t_vec[0]=get_unf_coo_x(rx, at_b+IBAS) - get_unf_coo_temp_x(at_c+IPHO);
      t_vec[1]=get_unf_coo_y(ry, at_b+IBAS) - get_unf_coo_temp_y(at_c+IPHO);
      t_vec[2]=get_unf_coo_z(rz, at_b+IBAS) - get_unf_coo_temp_z(at_c+IPHO);
      if(t_vec[0]*t_vec[0]+t_vec[1]*t_vec[1]+t_vec[2]*t_vec[2] <mc_ev_glob_rcut_sq){
	proj_on_nt(t_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
	proj_on_nt_inv(t_vec,rx,ry,rz, nt_neigh, r_vec_inv);
	if(is_purine(nt_neigh)){
	  b_ev_temp1=energy_hardcore(r_vec    , EV_PURPHO1, EV_PURPHO2, EV_PURPHO3,EV_PURPHO4, EV_PURPHO5, EV_PURPHO6); // base  wrt phosp
	  b_ev_temp2=energy_hardcore(r_vec_inv, EV_PHOPUR1, EV_PHOPUR2, EV_PHOPUR3,EV_PHOPUR4, EV_PHOPUR5, EV_PHOPUR6); // phosp wrt base
	}
	else{
	  b_ev_temp1=energy_hardcore(r_vec    , EV_PYRPHO1, EV_PYRPHO2, EV_PYRPHO3,EV_PYRPHO4, EV_PYRPHO5, EV_PYRPHO6); // base  wrt phosp
	  b_ev_temp2=energy_hardcore(r_vec_inv, EV_PHOPYR1, EV_PHOPYR2, EV_PHOPYR3,EV_PHOPYR4, EV_PHOPYR5, EV_PHOPYR6); // phosp wrt base
	}
	//}      else if(is_pyrimidine(nt_neigh)){
	//b_ev_temp1=energy_hardcore(r_vec    , EV_BASYPHO1, EV_BASYPHO2, EV_BASYPHO3,0,0,0); // base  wrt phosp
	//b_ev_temp2=energy_hardcore(r_vec_inv, EV_PHOBASY1, EV_PHOBASY2, EV_PHOBASY3,0,0,0); // phosp wrt base
	//}
	//printf("in pb %d  case 1, %lf %lf\t%lf %lf %lf   |  %lf %lf %lf\n", nt_c, b_ev_temp1, b_ev_temp2, r_vec[0], r_vec[1], r_vec[2], r_vec_inv[0], r_vec_inv[1], r_vec_inv[2]);
#ifdef WRMVERB
	if(b_ev_temp1>0 || b_ev_temp2>0) printf("# P %d\n# N %d\n", nt_c, nt_neigh);
#endif
#ifndef WARMUP
	if(b_ev_temp1>MC_CAP || b_ev_temp2>MC_CAP) {*mc_flag=6;
#ifdef XYZDEBUG
	  printf("NB CLASH (phosp-base): %d  %d\n", nt_c, nt_neigh);
#else
	  return bonded_energy;
#endif 	
	}//{printf("im in case 2\t %d %d \t %lf %lf %lf    %lf\n", nt_c, nt_neigh, t_vec[0], t_vec[1], t_vec[2], sqrt(t_vec[0]*t_vec[0]+t_vec[1]*t_vec[1]+t_vec[2]*t_vec[2]));*mc_flag=6;}
#endif
	bond_ev+=(b_ev_temp1+b_ev_temp2);
      }
      t_vec[0]=get_unf_coo_x(rx, at_b+IPHO) - get_unf_coo_temp_x(at_c+IBAS);
      t_vec[1]=get_unf_coo_y(ry, at_b+IPHO) - get_unf_coo_temp_y(at_c+IBAS);
      t_vec[2]=get_unf_coo_z(rz, at_b+IPHO) - get_unf_coo_temp_z(at_c+IBAS);
      proj_on_nt(t_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
      temp_flag=0;
      //printf("INTER  %d   ", nt_c);
      bond_bp=MC_calc_BP_inter(typ_ind, r_vec,mc_temp_glyc[nt_c], mc_temp_puck[nt_c], &temp_flag);
#ifdef WARMUP
      if(temp_flag==1){
	temp_flag=0;
	//energ_ret=MC_CAP-1;
	//typ_ind = N_BASES*mc_types[at_c]    + mc_types[at_ne];//N_BASES*type1+type2;
	this_typ=mc_types[at_c];
	phpos[0]=ph_pinter_ave[this_typ][mc_temp_glyc[nt_c]][mc_temp_puck[nt_c]][0];phpos[1]=ph_pinter_ave[this_typ][mc_temp_glyc[nt_c]][mc_temp_puck[nt_c]][1];phpos[2]=ph_pinter_ave[this_typ][mc_temp_glyc[nt_c]][mc_temp_puck[nt_c]][2];
	bond_bp=MC_CAP*((r_vec[0]-phpos[0])*(r_vec[0]-phpos[0])+(r_vec[1]-phpos[1])*(r_vec[1]-phpos[1])+(r_vec[2]-phpos[2])*(r_vec[2]-phpos[2]));
	//printf("   %lf %lf %lf    %lf %lf %lf\n", r_vec[0], r_vec[1], r_vec[2], phpos[0], phpos[1], phpos[2]);
	//*(r_vec[0]*r_vec[0]+r_vec[1]*r_vec[1]+r_vec[2]*r_vec[2]);typ_ind
	
      }
      
#endif
      
      
#ifdef WRMVERB
      if(bond_bp>=0) printf("# P %d\n# N %d\n", nt_c, nt_neigh);
#endif
      if(temp_flag==1){ *mc_flag=7;
#ifdef XYZDEBUG
	printf("BB CLASH (phosp-base, inter 1): %d  %d\t%lf %lf %lf  \t%lf %lf %lf\n", nt_neigh, nt_c, r_vec[0], r_vec[1], r_vec[2], r_vec_inv[0], r_vec_inv[1], r_vec_inv[2]);
#else
	return bonded_energy;
#endif 	   
      }
      //{printf("Inter PB interaction EXPLODED A\n%d %d   %lf %lf %lf\t%lf %lf %lf\n", nt_c, nt_neigh, t_vec[0], t_vec[1], t_vec[2], r_vec[0], r_vec[1], r_vec[2]);exit(-1);}
    }
    else if(nt_c>nt_neigh){//nt_up=nt_neigh;nt_do=nt_c;
      // EV between P of nt_neigh and B of nt_c
      // tab between P of nt_c     and B of nt_neigh
      t_vec[0]=get_unf_coo_temp_x(at_c+IBAS)-get_unf_coo_x(rx, at_b+IPHO);
      t_vec[1]=get_unf_coo_temp_y(at_c+IBAS)-get_unf_coo_y(ry, at_b+IPHO);
      t_vec[2]=get_unf_coo_temp_z(at_c+IBAS)-get_unf_coo_z(rz, at_b+IPHO);
      if(t_vec[0]*t_vec[0]+t_vec[1]*t_vec[1]+t_vec[2]*t_vec[2] <mc_ev_glob_rcut_sq){
	proj_on_nt(t_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
	proj_on_nt_inv(t_vec, rx,ry,rz,nt_neigh, r_vec_inv);
	if(is_purine(nt_c)){
	  b_ev_temp1=energy_hardcore(r_vec    , EV_PHOPUR1, EV_PHOPUR2, EV_PHOPUR3,EV_PHOPUR4, EV_PHOPUR5, EV_PHOPUR6); // base  wrt phosp  - we are now projecting on the SR of the base!
	  b_ev_temp2=energy_hardcore(r_vec_inv, EV_PURPHO1, EV_PURPHO2, EV_PURPHO3,EV_PURPHO4, EV_PURPHO5, EV_PURPHO6); // phosp wrt base
	}
	else{
	  b_ev_temp1=energy_hardcore(r_vec    , EV_PHOPYR1, EV_PHOPYR2, EV_PHOPYR3,EV_PHOPYR4, EV_PHOPYR5, EV_PHOPYR6); // base  wrt phosp  - we are now projecting on the SR of the base!
	  b_ev_temp2=energy_hardcore(r_vec_inv, EV_PYRPHO1, EV_PYRPHO2, EV_PYRPHO3,EV_PYRPHO4, EV_PYRPHO5, EV_PYRPHO6); // phosp wrt base
	}
	
	
	
#ifdef WRMVERB
	if(b_ev_temp1>0 || b_ev_temp2>0) printf("# N %d\n# P %d\n", nt_c, nt_neigh);
#endif
#ifndef WARMUP
	if(b_ev_temp1>MC_CAP || b_ev_temp2>MC_CAP) {*mc_flag=8;
#ifdef XYZDEBUG
	  printf("NB CLASH (base-phosp): %d  %d\n", nt_c, nt_neigh);
#else
	  return bonded_energy;
#endif 	
	  
	}//{printf("im in case 2\t %d %d \t %lf %lf %lf    %lf\n", nt_c, nt_neigh, t_vec[0], t_vec[1], t_vec[2], sqrt(t_vec[0]*t_vec[0]+t_vec[1]*t_vec[1]+t_vec[2]*t_vec[2]));*mc_flag=6;}
#endif
	bond_ev+=(b_ev_temp1+b_ev_temp2);
      }
      t_vec[0]=get_unf_coo_x(rx, at_b+IBAS) - get_unf_coo_temp_x(at_c+IPHO);
      t_vec[1]=get_unf_coo_y(ry, at_b+IBAS) - get_unf_coo_temp_y(at_c+IPHO);
      t_vec[2]=get_unf_coo_z(rz, at_b+IBAS) - get_unf_coo_temp_z(at_c+IPHO);
      proj_on_nt_inv(t_vec,rx,ry,rz, nt_neigh, r_vec_inv);
      temp_flag=0;
      //printf("INTER  %d   ", nt_c);
      bond_bp=MC_calc_BP_inter(typ_ind2, r_vec_inv, mc_glyc[nt_neigh], mc_puck[nt_neigh], &temp_flag);

#ifdef WARMUP
      if(temp_flag==1){
	temp_flag=0;
	//energ_ret=MC_CAP-1;
	//typ_ind = N_BASES*mc_types[at_c]    + mc_types[at_ne];//N_BASES*type1+type2;
	this_typ=mc_types[at_b];
	phpos[0]=ph_pinter_ave[this_typ][mc_temp_glyc[nt_neigh]][mc_temp_puck[nt_neigh]][0];phpos[1]=ph_pinter_ave[this_typ][mc_temp_glyc[nt_neigh]][mc_temp_puck[nt_neigh]][1];phpos[2]=ph_pinter_ave[this_typ][mc_temp_glyc[nt_neigh]][mc_temp_puck[nt_neigh]][2];
	bond_bp=MC_CAP*((r_vec[0]-phpos[0])*(r_vec[0]-phpos[0])+(r_vec[1]-phpos[1])*(r_vec[1]-phpos[1])+(r_vec[2]-phpos[2])*(r_vec[2]-phpos[2]));
	//printf("   %lf %lf %lf    %lf %lf %lf\n", r_vec[0], r_vec[1], r_vec[2], phpos[0], phpos[1], phpos[2]);
	//*(r_vec[0]*r_vec[0]+r_vec[1]*r_vec[1]+r_vec[2]*r_vec[2]);typ_ind
      }
      
#endif


#ifdef WRMVERB
      if(bond_bp>=0) printf("# P %d\n# N %d\n", nt_c, nt_neigh);
#endif  
      if(temp_flag==1){*mc_flag=9; 
#ifdef XYZDEBUG
	printf("BB CLASH (base-phosp, inter 2): %d  %d\t%lf %lf %lf  \t%lf %lf %lf\n", nt_neigh, nt_c, r_vec[0], r_vec[1], r_vec[2], r_vec_inv[0], r_vec_inv[1], r_vec_inv[2]);
#else
	return bonded_energy;
#endif
      }//{printf("Inter PS interaction EXPLODED B\n%d %d   %lf %lf %lf\t%lf %lf %lf  %d %d\n", nt_c, nt_neigh, t_vec[0], t_vec[1], t_vec[2], r_vec[0], r_vec[1], r_vec[2], mc_types[N_PARTS_PER_NT*nt_c]    ,mc_types[N_PARTS_PER_NT*nt_neigh]			      );exit(-1);}
      
      
      
    }
    else {printf("BOND BROKEN\n");exit(-2);}
    // if(nt_c-nt_neigh == -1 ) { // neighbor is on top of central - phosp neighbor and sugar central
    /* t_vec[0]=(rx[ba+IPHO]+mc_pbox[ba+IPHO][0]*box_l[0]) - (mc_temp_x[ISUG]+mc_temp_pbox[ISUG][0]*box_l[0]); */
    /* t_vec[1]=(ry[ba+IPHO]+mc_pbox[ba+IPHO][1]*box_l[1]) - (mc_temp_y[ISUG]+mc_temp_pbox[ISUG][1]*box_l[1]); */
    /* t_vec[2]=(rz[ba+IPHO]+mc_pbox[ba+IPHO][2]*box_l[2]) - (mc_temp_z[ISUG]+mc_temp_pbox[ISUG][2]*box_l[2]); */
    /* calc_rel_pos(t_vec, r_vec); */
    /* calc_rel_pos_inv(t_vec, nt_neigh,rx,ry,rz, r_vec_inv); */
    /* b_ev_temp1=energy_hardcore(r_vec    , SIG_PHOSUG, SIG_PHOSUG, SIG_PHOSUG); // phosp wrt sugar */
    /* b_ev_temp2=energy_hardcore(r_vec_inv, SIG_PHOSUG, SIG_PHOSUG, SIG_PHOSUG); // sugar wrt phosp */
    /* if(b_ev_temp1>MC_CAP || b_ev_temp2>MC_CAP) *mc_flag=6;// {printf("im in case 1\n"); *mc_flag=6;} */
    /* bond_ev+=(b_ev_temp1+b_ev_temp2); */
    //}
    //else if(nt_c - nt_neigh == 1){
    
    //    }
    /** PHOSPHATE-PHOSPHATE **/
    //calc_min_vec(rx[N_PARTS_PER_NT*nt_neigh+IPHO], ry[N_PARTS_PER_NT*nt_neigh+IPHO], rz[N_PARTS_PER_NT*nt_neigh+IPHO], mc_temp_x[IPHO], mc_temp_y[IPHO], mc_temp_z[IPHO], t_vec, &r);
    t_vec[0]=get_unf_coo_x(rx, at_b+IPHO) - get_unf_coo_temp_x(at_c+IPHO);
    //(rx[ba+IPHO]+mc_pbox[ba+IPHO][0]*box_l[0]) - (mc_temp_x[IPHO]+mc_temp_pbox[IPHO][0]*box_l[0]);
    t_vec[1]=get_unf_coo_y(ry, at_b+IPHO) - get_unf_coo_temp_y(at_c+IPHO);
    t_vec[2]=get_unf_coo_z(rz, at_b+IPHO) - get_unf_coo_temp_z(at_c+IPHO);
    if(t_vec[0]*t_vec[0]+t_vec[1]*t_vec[1]+t_vec[2]*t_vec[2] <mc_ev_glob_rcut_sq){
      proj_on_nt(t_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
      proj_on_nt_inv(t_vec,rx,ry,rz,  nt_neigh,r_vec_inv);
      b_ev_temp1=energy_hardcore(r_vec    , EV_PHOS1, EV_PHOS2, EV_PHOS3,EV_PHOS4, EV_PHOS5, EV_PHOS6); // sugar wrt phosp
      b_ev_temp2=energy_hardcore(r_vec_inv, EV_PHOS1, EV_PHOS2, EV_PHOS3,EV_PHOS4, EV_PHOS5, EV_PHOS6); // phosp wrt sugar
#ifdef WRMVERB
      if(b_ev_temp1>0 || b_ev_temp2>0) printf("# P %d\n# P %d\n", nt_c, nt_neigh);
#endif
#ifndef WARMUP
      if(b_ev_temp1>MC_CAP || b_ev_temp2>MC_CAP) {*mc_flag=10;
#ifdef XYZDEBUG
	printf("NB CLASH (phosp-phosp): %d  %d\n", nt_c,nt_neigh);
#else
	return bonded_energy;
#endif 	 
      }
#endif
      bond_ev+=(b_ev_temp1+b_ev_temp2);
    }
    //printf("bb: %lf   %lf\n", bond_bb, bond_ev);
    //printf("bp  %d : %lf\n", nt_c, bond_bp);
    
    //HEREFLAG
    backbone_e=bond_ss + bond_bp + bond_ev;
#ifdef NEW_BIA
    if(mc_target_temp>1.0)
      backbone_e*=(mc_target_temp);
#endif
    bonded_energy+=backbone_e;
    //} 
    //else{
    //if(bond_st<0) bonded_energy-=nb_st_well[typ_ind];
    //}
    //nbe-=nb_wc_well[typ_ind];
    //printf("%lf %lf %lf   %lf    %d\n", bond_bb, bond_ev, bond_st, bonded_energy, temp_flag);
    //THE VECTOR MUST BE THE REAL DISTANCE BETWEEN THE TWO PARTICLES, NOT THE MINIMUM.
    /* d_vec[0]=(rx[ba]+mc_pbox[ba][0]*box_l[0]) - (mc_temp_x[0]+mc_temp_pbox[0][0]*box_l[0]); */
    /* d_vec[1]=(ry[ba]+mc_pbox[ba][1]*box_l[1]) - (mc_temp_y[0]+mc_temp_pbox[0][1]*box_l[1]); */
    /* d_vec[2]=(rz[ba]+mc_pbox[ba][2]*box_l[2]) - (mc_temp_z[0]+mc_temp_pbox[0][2]*box_l[2]); */
    /* dist=sqrt(d_vec[0]*d_vec[0]+d_vec[1]*d_vec[1]+d_vec[2]*d_vec[2]); */
    /* hk=mc_harm_k[mc_bondlist[nt_c][mc_nbonds[nt_c][0]+b]]; */
    /* hr=mc_harm_r[mc_bondlist[nt_c][mc_nbonds[nt_c][0]+b]]; */
    /* e_bond=energy_harmonic(dist, hr, hk); */
    /* bonded_energy+=e_bond; */
  }

  return bonded_energy;
}

int MC_get_pucks(int p1, int p2){
  int RET;
  if(p1==PUCK_3 && p2==PUCK_3) RET= PUCKS_33;
  if(p1==PUCK_3 && p2==PUCK_2) RET= PUCKS_32;
  if(p1==PUCK_2 && p2==PUCK_3) RET= PUCKS_23;
  if(p1==PUCK_2 && p2==PUCK_2) RET= PUCKS_22;
  return RET; 
}


int MC_get_glycs(int p1, int p2){
  //int p1=mc_glyc[nt_1];
  //int p2=mc_glyc[nt_2];
  if(     p1==GLYC_A && p2==GLYC_A) return GLYCS_AA;
  else if(p1==GLYC_A && p2==GLYC_H) return GLYCS_AH;
  else if(p1==GLYC_A && p2==GLYC_S) return GLYCS_AS;
  else if(p1==GLYC_H && p2==GLYC_A) return GLYCS_HS;
  else if(p1==GLYC_H && p2==GLYC_H) return GLYCS_HH;
  else if(p1==GLYC_H && p2==GLYC_S) return GLYCS_HS;
  else if(p1==GLYC_S && p2==GLYC_A) return GLYCS_SS;
  else if(p1==GLYC_S && p2==GLYC_H) return GLYCS_SH;
  else if(p1==GLYC_S && p2==GLYC_S) return GLYCS_SS;
  
  else {printf("Undefined GLYC states at %d and %d!\n", p1, p2);exit(100);}
}

double energy_harmonic(double r, double r0, double k){
  return 0.5*k*(r-r0)*(r-r0);
}

double energy_angle(double t, double t0, double k){
  return 0.5*k*(t-t0)*(t-t0);
}



double MC_calc_non_bonded_energy(int nt_c, double *rx, double *ry, double *rz, int nt_neigh, double *d_vec, double dist, int *mc_flag){
  int at_c=N_PARTS_PER_NT*nt_c;
  int at_ne=N_PARTS_PER_NT*nt_neigh;
  double nbe=0, nbe_wc, nbe_st, nb_ev=0;
  double r_vec[DIM], r_vec_inv[DIM], t_vec[DIM], theta, r;
  double sigma[DIM];
  double nb_ev_temp1, nb_ev_temp2;
  int temp_flag_wc=0, temp_flag_st=0;
  int typ_ind, typ_ind2;
  double ori_o[DIM], ori_n[DIM];
  double eta2;
  //double  st_theta_0=180.0/180.0*M_PI, st_theta;
  *mc_flag=0;
  //base = 0, sugar = 1, phosphate = 2
  //nb_ev_temp1=energy_hardcore(r_vec,     sigma1, sigma2, sigma3);
  /** BASE-SUGAR **/
  calc_min_vec(rx[at_ne+ISUG], ry[at_ne+ISUG], rz[at_ne+ISUG], mc_temp_x[at_c], mc_temp_y[at_c], mc_temp_z[at_c], t_vec, &r);
  if(r<EV_GLOB_RCUT){
    proj_on_nt(t_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
    proj_on_nt_inv(t_vec, rx,ry,rz, nt_neigh,r_vec_inv);
    if(is_purine(nt_c)){
      nb_ev_temp1=energy_hardcore(r_vec    , EV_SUGPUR1, EV_SUGPUR2, EV_SUGPUR3, EV_SUGPUR4, EV_SUGPUR5, EV_SUGPUR6 ); // sugar wrt base
      nb_ev_temp2=energy_hardcore(r_vec_inv, EV_PURSUG1, EV_PURSUG2, EV_PURSUG3, EV_PURSUG4, EV_PURSUG5, EV_PURSUG6 );//EV_BASRSUG4, EV_BASRSUG5, EV_BASRSUG6 ); // base  wrt sugar
    }
    else{
      nb_ev_temp1=energy_hardcore(r_vec    , EV_SUGPYR1, EV_SUGPYR2, EV_SUGPYR3, EV_SUGPYR4, EV_SUGPYR5, EV_SUGPYR6 ); // sugar wrt base
      nb_ev_temp2=energy_hardcore(r_vec_inv, EV_PYRSUG1, EV_PYRSUG2, EV_PYRSUG3, EV_PYRSUG4, EV_PYRSUG5, EV_PYRSUG6 );//EV_BASRSUG4, EV_BASRSUG5, EV_BASRSUG6 ); // base  wrt sugar
    }
#ifdef WRMVERB
    if(nb_ev_temp1>0 || nb_ev_temp2>0) printf("# N %d\n# N %d\n", nt_c, nt_neigh);
#endif
    #ifndef WARMUP
    if(nb_ev_temp1>MC_CAP || nb_ev_temp2>MC_CAP) {*mc_flag=11;
#ifdef XYZDEBUG
      printf("NB CLASH (sugar-base): %d  %d\n", nt_c,nt_neigh);
#endif
    }
    #endif
    nb_ev+=(nb_ev_temp1+nb_ev_temp2);
  }
  /** SUGAR-BASE **/
  calc_min_vec(rx[at_ne], ry[at_ne], rz[at_ne], mc_temp_x[at_c+ISUG], mc_temp_y[at_c+ISUG], mc_temp_z[at_c+ISUG], t_vec, &r);
  if(r<EV_GLOB_RCUT){
    proj_on_nt(t_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
    proj_on_nt_inv(t_vec,rx,ry,rz, nt_neigh, r_vec_inv);
    if(is_purine(nt_neigh)){
      nb_ev_temp1=energy_hardcore(r_vec    , EV_PURSUG1, EV_PURSUG2, EV_PURSUG3, EV_PURSUG4, EV_PURSUG5, EV_PURSUG6); //EV_PURRSUG4, EV_PURRSUG5, EV_PURRSUG6); // base  wrt sugar
      nb_ev_temp2=energy_hardcore(r_vec_inv, EV_SUGPUR1, EV_SUGPUR2, EV_SUGPUR3, EV_SUGPUR4, EV_SUGPUR5, EV_SUGPUR6); // sugar wrt base
    }
    else{
      nb_ev_temp1=energy_hardcore(r_vec    , EV_PYRSUG1, EV_PYRSUG2, EV_PYRSUG3, EV_PYRSUG4, EV_PYRSUG5, EV_PYRSUG6); //EV_PURRSUG4, EV_PURRSUG5, EV_PURRSUG6); // base  wrt sugar
      nb_ev_temp2=energy_hardcore(r_vec_inv, EV_SUGPYR1, EV_SUGPYR2, EV_SUGPYR3, EV_SUGPYR4, EV_SUGPYR5, EV_SUGPYR6); // sugar wrt base
    }
#ifdef WRMVERB
    if(nb_ev_temp1>0 || nb_ev_temp2>0) printf("# N %d\n# N %d\n", nt_c, nt_neigh);
#endif
#ifndef WARMUP 
    if(nb_ev_temp1>MC_CAP || nb_ev_temp2>MC_CAP) {*mc_flag=12;
#ifdef XYZDEBUG
      printf("NB CLASH (base-sugar): %d  %d\n", nt_c,nt_neigh);
#endif 
    }
#endif
    nb_ev+=(nb_ev_temp1+nb_ev_temp2);
  }
  /** SUGAR-SUGAR **/
  calc_min_vec(rx[at_ne+ISUG], ry[at_ne+ISUG], rz[at_ne+ISUG], mc_temp_x[at_c+ISUG], mc_temp_y[at_c+ISUG], mc_temp_z[at_c+ISUG], t_vec, &r);
  if(r<EV_GLOB_RCUT){
    proj_on_nt(t_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
    proj_on_nt_inv(t_vec,rx,ry,rz, nt_neigh, r_vec_inv);
    nb_ev_temp1=energy_hardcore(r_vec    , EV_SUGSUG1, EV_SUGSUG2, EV_SUGSUG3,  EV_SUGSUG4, EV_SUGSUG5, EV_SUGSUG6); // base  wrt sugar
    nb_ev_temp2=energy_hardcore(r_vec_inv, EV_SUGSUG1, EV_SUGSUG2, EV_SUGSUG3,  EV_SUGSUG4, EV_SUGSUG5, EV_SUGSUG6); // sugar wrt base
#ifdef WRMVERB
    if(nb_ev_temp1>0 || nb_ev_temp2>0) printf("# N %d\n# N %d\n", nt_c, nt_neigh);
#endif
#ifndef WARMUP
    if(nb_ev_temp1>MC_CAP || nb_ev_temp2>MC_CAP) {*mc_flag=13;
#ifdef XYZDEBUG
      printf("NB CLASH (sugar-sugar): %d  %d\n", nt_c,nt_neigh);
#endif
    }
#endif
    nb_ev+=(nb_ev_temp1+nb_ev_temp2);
  }
  /** PHOSPHATE-BASE **/
  calc_min_vec(rx[at_ne], ry[at_ne], rz[at_ne], mc_temp_x[at_c+IPHO], mc_temp_y[at_c+IPHO], mc_temp_z[at_c+IPHO], t_vec, &r);
  if(r<EV_GLOB_RCUT){
    proj_on_nt(t_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
    proj_on_nt_inv(t_vec,rx,ry,rz, nt_neigh, r_vec_inv);
    if(is_purine(nt_neigh)){
      nb_ev_temp1=energy_hardcore(r_vec    , EV_PURPHO1, EV_PURPHO2, EV_PURPHO3, EV_PURPHO4, EV_PURPHO5, EV_PURPHO6 ); // base  wrt phosp
      nb_ev_temp2=energy_hardcore(r_vec_inv, EV_PHOPUR1, EV_PHOPUR2, EV_PHOPUR3, EV_PHOPUR4, EV_PHOPUR5, EV_PHOPUR6 ); // phosp wrt base
    }
    else{
      nb_ev_temp1=energy_hardcore(r_vec    , EV_PYRPHO1, EV_PYRPHO2, EV_PYRPHO3, EV_PYRPHO4, EV_PYRPHO5, EV_PYRPHO6 ); // base  wrt phosp
      nb_ev_temp2=energy_hardcore(r_vec_inv, EV_PHOPYR1, EV_PHOPYR2, EV_PHOPYR3, EV_PHOPYR4, EV_PHOPYR5, EV_PHOPYR6 ); // phosp wrt base
    }
#ifdef WRMVERB
    if(nb_ev_temp1>0 || nb_ev_temp2>0) printf("# P %d\n# N %d\n", nt_c, nt_neigh);
#endif
#ifndef WARMUP
    if(nb_ev_temp1>MC_CAP || nb_ev_temp2>MC_CAP) {*mc_flag=14;
#ifdef XYZDEBUG
      printf("NB CLASH (phosp-sugar): %d  %d\n", nt_c,nt_neigh);
#endif
    }
#endif
    nb_ev+=(nb_ev_temp1+nb_ev_temp2);
  }
  /** BASE-PHOSPHATE **/
  calc_min_vec(rx[at_ne+IPHO], ry[at_ne+IPHO], rz[at_ne+IPHO], mc_temp_x[at_c], mc_temp_y[at_c], mc_temp_z[at_c], t_vec, &r);
  if(r<EV_GLOB_RCUT){
    proj_on_nt(t_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
    proj_on_nt_inv(t_vec, rx,ry,rz,nt_neigh, r_vec_inv);
    if(is_purine(nt_c)){
      nb_ev_temp1=energy_hardcore(r_vec    , EV_PHOPUR1, EV_PHOPUR2, EV_PHOPUR3, EV_PHOPUR4, EV_PHOPUR5, EV_PHOPUR6 ); // phosp wrt base
      nb_ev_temp2=energy_hardcore(r_vec_inv, EV_PURPHO1, EV_PURPHO2, EV_PURPHO3, EV_PURPHO4, EV_PURPHO5, EV_PURPHO6 ); // base  wrt phosp
    }
    else{
      nb_ev_temp1=energy_hardcore(r_vec    , EV_PHOPYR1, EV_PHOPYR2, EV_PHOPYR3, EV_PHOPYR4, EV_PHOPYR5, EV_PHOPYR6 ); // phosp wrt base
      nb_ev_temp2=energy_hardcore(r_vec_inv, EV_PYRPHO1, EV_PYRPHO2, EV_PYRPHO3, EV_PYRPHO4, EV_PYRPHO5, EV_PYRPHO6 ); // base  wrt phosp
    }
    
#ifdef WRMVERB
    if(nb_ev_temp1>0 || nb_ev_temp2>0) printf("# P %d\n# N %d\n", nt_neigh, nt_c);
#endif

#ifndef WARMUP
    if(nb_ev_temp1>MC_CAP || nb_ev_temp2>MC_CAP) {*mc_flag=15;
#ifdef XYZDEBUG
      printf("NB CLASH (phosp-base): %d  %d\t%lf %lf\n", nt_c,nt_neigh, nb_ev_temp1, nb_ev_temp2);
#endif
    }
#endif
    nb_ev+=(nb_ev_temp1+nb_ev_temp2);
  }
  /** PHOSPHATE-SUGAR **/
  calc_min_vec(rx[at_ne+ISUG], ry[at_ne+ISUG], rz[at_ne+ISUG], mc_temp_x[at_c+IPHO], mc_temp_y[at_c+IPHO], mc_temp_z[at_c+IPHO], t_vec, &r);
  if(r<EV_GLOB_RCUT){
    proj_on_nt(t_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
    proj_on_nt_inv(t_vec, rx,ry,rz,nt_neigh, r_vec_inv);
    nb_ev_temp1=energy_hardcore(r_vec    , EV_SUGPHO1, EV_SUGPHO2, EV_SUGPHO3, EV_SUGPHO4, EV_SUGPHO5, EV_SUGPHO6); // base wrt pho
    nb_ev_temp2=energy_hardcore(r_vec_inv, EV_PHOSUG1, EV_PHOSUG2, EV_PHOSUG3, EV_PHOSUG4, EV_PHOSUG5, EV_PHOSUG6); // phosp wrt sugar
    
    //nb_ev_temp1=energy_hardcore(r_vec    , EV_PHOSUG, EV_PHOSUG, EV_PHOSUG); // sugar wrt phosp
    //nb_ev_temp2=energy_hardcore(r_vec_inv, EV_PHOSUG, EV_PHOSUG, EV_PHOSUG); // phosp wrt sugar
#ifdef WRMVERB
    if(nb_ev_temp1>0 || nb_ev_temp2>0) printf("# P %d\n# N %d\n", nt_c, nt_neigh);
#endif
#ifndef WARMUP
    if(nb_ev_temp1>MC_CAP || nb_ev_temp2>MC_CAP) {*mc_flag=16;
#ifdef XYZDEBUG
      printf("NB CLASH (phosp-sugar): %d  %d\n# ", nt_c,nt_neigh);
#endif
    }
#endif
    nb_ev+=(nb_ev_temp1+nb_ev_temp2);
  }
  /** SUGAR-PHOSPHATE **/
  calc_min_vec(rx[at_ne+IPHO], ry[at_ne+IPHO], rz[at_ne+IPHO], mc_temp_x[at_c+ISUG], mc_temp_y[at_c+ISUG], mc_temp_z[at_c+ISUG], t_vec, &r);
  if(r<EV_GLOB_RCUT){
    proj_on_nt(t_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
    proj_on_nt_inv(t_vec,rx,ry,rz, nt_neigh, r_vec_inv);
    nb_ev_temp1=energy_hardcore(r_vec    , EV_PHOSUG1, EV_PHOSUG2, EV_PHOSUG3, EV_PHOSUG4, EV_PHOSUG5, EV_PHOSUG6); // phosp wrt sugar
    nb_ev_temp2=energy_hardcore(r_vec_inv, EV_SUGPHO1, EV_SUGPHO2, EV_SUGPHO3, EV_SUGPHO4, EV_SUGPHO5, EV_SUGPHO6); // base wrt pho
    //nb_ev_temp1=energy_hardcore(r_vec    , EV_PHOSUG, EV_PHOSUG, EV_PHOSUG); // sugar wrt phosp
    //nb_ev_temp2=energy_hardcore(r_vec_inv, EV_PHOSUG, EV_PHOSUG, EV_PHOSUG); // phosp wrt sugar
#ifdef WRMVERB
    if(nb_ev_temp1>0 || nb_ev_temp2>0) printf("# N %d\n# P %d\n", nt_c, nt_neigh);
#endif
#ifndef WARMUP
    if(nb_ev_temp1>MC_CAP || nb_ev_temp2>MC_CAP) {*mc_flag=17;
#ifdef XYZDEBUG
      printf("NB CLASH (sugar-phosp): %d  %d\n", nt_c,nt_neigh);
#endif    
    }
#endif
    nb_ev+=(nb_ev_temp1+nb_ev_temp2);
  }
  /** PHOSPHATE-PHOSPHATE **/
  calc_min_vec(rx[at_ne+IPHO], ry[at_ne+IPHO], rz[at_ne+IPHO], mc_temp_x[at_c+IPHO], mc_temp_y[at_c+IPHO], mc_temp_z[at_c+IPHO], t_vec, &r);
  if(r<EV_GLOB_RCUT){
    proj_on_nt(t_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
    proj_on_nt_inv(t_vec, rx,ry,rz, nt_neigh,r_vec_inv);
    nb_ev_temp1=energy_hardcore(r_vec    , EV_PHOS1, EV_PHOS2, EV_PHOS3,EV_PHOS4, EV_PHOS5, EV_PHOS6 ); // sugar wrt phosp
    nb_ev_temp2=energy_hardcore(r_vec_inv, EV_PHOS1, EV_PHOS2, EV_PHOS3,EV_PHOS4, EV_PHOS5, EV_PHOS6 ); // phosp wrt sugar
#ifdef WRMVERB
    if(nb_ev_temp1>0 || nb_ev_temp2>0) printf("# P %d\n# P %d\n", nt_c, nt_neigh);
#endif
#ifndef WARMUP
    if(nb_ev_temp1>MC_CAP || nb_ev_temp2>MC_CAP) {*mc_flag=18;
#ifdef XYZDEBUG
      printf("NB CLASH (phosp-phosp): %d  %d\n", nt_c,nt_neigh);
#endif
    }
#endif
    nb_ev+=(nb_ev_temp1+nb_ev_temp2);
  }
#ifdef STWCBONDED
  return nbe;
#endif
  if(dist<EV_GLOB_RCUT){ 
    /** BASE-BASE - STACKING AND EV  - WC IS CALCULATED ELSEWHERE **/
    temp_flag_wc=0;
    typ_ind = N_BASES*mc_types[at_c]    + mc_types[at_ne];//N_BASES*type1+type2;
    typ_ind2=N_BASES*mc_types[at_ne]+ mc_types[at_c];
    proj_on_nt(d_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
    proj_on_nt_inv(d_vec,rx, ry,rz, nt_neigh, r_vec_inv);
    /* theta=calc_wc_psdihedr(mc_temp_x, mc_temp_y, mc_temp_z,rx,ry,rz,nt_c,nt_neigh);  */
    /* nbe_wc = MC_calc_nnN_watscric(typ_ind, r_vec, r_vec_inv, theta, &temp_flag_wc); */
    /* if(temp_flag_wc==0) */
    /*   nbe+=nbe_wc; */
    
    temp_flag_st=0;
    //printf("%d %d    %lf %lf %lf   %lf %lf %lf       %lf   \t", nt_c, nt_neigh, r_vec[0], r_vec[1], r_vec[2], r_vec_inv[0], r_vec_inv[1], r_vec_inv[2], sqrt(r_vec[0]*r_vec[0]+r_vec[1]*r_vec[1]+r_vec[2]*r_vec[2]));
    eta2=calc_wc_psdihedr(mc_temp_x, mc_temp_y, mc_temp_z, rx,ry,rz,nt_c,nt_neigh); 
    /* STACKING DEPENDS ON THE GLYC CONFORMATION - AA IS BONDED, ANYTHING ELSE IS NONBONDED */ 
    nbe_st=0;
    if(fabs(eta2)<STANG || fabs(eta2) > M_PI-STANG)
      nbe_st = MC_calc_nnN_stacking(typ_ind,typ_ind2, r_vec, r_vec_inv, &temp_flag_st); // NO PAIR SPECIFIC!!
  
  /* norex_vec_prod(dist_1d(mc_temp_x[1], mc_temp_x[0],0),dist_1d(mc_temp_y[1], mc_temp_y[0],1),dist_1d(mc_temp_z[1], mc_temp_z[0],2),dist_1d(mc_temp_x[2], mc_temp_x[0],0),dist_1d(mc_temp_y[2], mc_temp_y[0],1),dist_1d(mc_temp_z[2], mc_temp_z[0],2),ori_n); */
  /* norex_vec_prod(dist_1d(rx[nt_neigh*N_PARTS_PER_NT+1], rx[nt_neigh*N_PARTS_PER_NT+0],0),dist_1d(ry[nt_neigh*N_PARTS_PER_NT+1], ry[nt_neigh*N_PARTS_PER_NT+0],1),dist_1d(rz[nt_neigh*N_PARTS_PER_NT+1], rz[nt_neigh*N_PARTS_PER_NT+0],2),dist_1d(rx[nt_neigh*N_PARTS_PER_NT+2], rx[nt_neigh*N_PARTS_PER_NT+0],0),dist_1d(ry[nt_neigh*N_PARTS_PER_NT+2], ry[nt_neigh*N_PARTS_PER_NT+0],1),dist_1d(rz[nt_neigh*N_PARTS_PER_NT+2], rz[nt_neigh*N_PARTS_PER_NT+0],2),ori_o); */
  /* st_theta=acos(dot_prod(ori_n, ori_o)); */
  /* if(st_theta>st_theta_0) temp_flag_st=1; */
  
  if(temp_flag_st==0)
    nbe+=nbe_st;
  
  nb_ev_temp1=energy_hardcore(r_vec,     EV_BASE1, EV_BASE2, EV_BASE3,0,0,0);
  nb_ev_temp2=energy_hardcore(r_vec_inv, EV_BASE1, EV_BASE2, EV_BASE3,0,0,0);
#ifdef WRMVERB
    if(nb_ev_temp1>0 || nb_ev_temp2>0) printf("# N %d\n# N %d\n", nt_c, nt_neigh);
#endif  
  //if((nb_ev_temp1>MC_CAP ||nb_ev_temp2>MC_CAP) && temp_flag_wc != 0 && temp_flag_st != 0) *mc_flag=16;
#ifndef WARMUP
  if((nb_ev_temp1>MC_CAP ||nb_ev_temp2>MC_CAP)) {*mc_flag=19;
#ifdef XYZDEBUG
    printf("NB CLASH (base-base): %d  %d\n", nt_c,nt_neigh);
#endif
  }
#endif
  if(temp_flag_st !=0 && temp_flag_wc != 0)
    nb_ev+=(nb_ev_temp1+nb_ev_temp2);
  
}
  return nbe+nb_ev;
}									


int get_stacking_ind(int a, int b){
  int RET;
  if(a==STFACE3 && b==STFACE3)      RET=STFACES33;
  else if(a==STFACE3 && b==STFACE5) RET=STFACES35;
  else if(a==STFACE5 && b==STFACE3) RET=STFACES53;
  else if(a==STFACE5 && b==STFACE5) RET=STFACES55;
  else {printf("Stacking indexes undefined!\n");RET=0;exit(1);}
  return RET;
}


double MC_calc_nnB_stacking(int typ_ind, int typ_ind_inv, double *r_vec, double *r_vec_inv, double eta, int *flag){
  double be=0, bei=0, bed=0;
  double tot_b_st;
  
  int sface_c, sface_n;
  if(r_vec[2]>0)     sface_c=STFACE3 ; else sface_c=STFACE5;
  if(r_vec_inv[2]>0) sface_n=STFACE3 ; else sface_n=STFACE5;
  
  int stind    =get_stacking_ind(sface_c, sface_n);
  int stind_inv=get_stacking_ind(sface_n, sface_c);
  
  
  be=calc_nnB_tab_stacking(r_vec, typ_ind, stind);
  bei=calc_nnB_tab_stacking(r_vec_inv, typ_ind_inv, stind_inv);
  bed=0;//calc_nnB_tab_st_psdihedr(eta, typ_ind, stind); // stacking pseudodihedral 
  //if(be==0 || bei==0 || bed==0)
  if(be==0 || bei==0)
    *flag=1;

  /***************** THIS IS SUPER PROVISORY **************/
  tot_b_st=0;
  //if(be!=0 && bei!=0 && bed!=0 && b_st_well[typ_ind][stind]!=0){
  if(be!=0 && bei!=0 && b_st_well[typ_ind][stind]!=0){
    tot_b_st=be+bei+bed;
    tot_b_st-=b_st_well[typ_ind][stind];
#ifdef BI_ANNEALING
    tot_b_st*=bia_lambda;
    if(tot_b_st<bia_cap)
      tot_b_st=bia_cap;
#endif
  }
  return tot_b_st;
}

double MC_calc_nnN_stacking(int typ_ind, int typ_ind_2, double *r_vec, double *r_vec_inv, int *flag){
  double nbe1=0, nbe2=0, nbe1i=0, nbe2i=0;
  double nbe=0, nbei=0;
  
  if(r_vec[2]>0) nbe1=calc_nnN_tab_stacking_s3(r_vec, typ_ind);
  else nbe1=calc_nnN_tab_stacking_s5(r_vec, typ_ind);
  
  if(r_vec_inv[2]>0) nbe2=calc_nnN_tab_stacking_s3(r_vec_inv, typ_ind_2);
  else nbe2=calc_nnN_tab_stacking_s5(r_vec_inv, typ_ind_2);
  
  if(nbe1==0 || nbe2==0) *flag=1; //in this case, there is no interaction energy
    
  
  /* if(table_nnN_N_0[typ_ind][0]>0){ */
  /*   nbe1=calc_nnN_tab_stacking(r_vec,typ_ind); // stacking */
  /*   nbe2=calc_nnN_tab_stacking_inv(r_vec,typ_ind); // stacking */
  /*   if(nbe1 == 0 && nbe2==0){ */
  /*     *flag=1; */
  /*   } */
  /*   nbe=nbe1+nbe2; */
  /*   if(nbe1!=0 && nbe2!=0) {printf ("NB STACKING NOT WELL DEFINED!!\n"); exit(1);} */
  /* } */
  /* if(table_nnN_N_0_inv[typ_ind_2][0]>0){ */
  /*   nbe1i=calc_nnN_tab_stacking_inv(r_vec_inv,typ_ind_2); // stacking */
  /*   nbe2i=calc_nnN_tab_stacking(r_vec_inv,typ_ind_2); // stacking */

  /*   if(nbe1i == 0 && nbe2i == 0){ */
  /*     *flag=1; */

  /*  } */

  /*   nbei=nbe1i+nbe2i; */
  /*   if(nbe1i!=0 && nbe2i!=0) {printf ("NB STACKING NOT WELL DEFINED!!\n"); exit(1);} */
  /* } */

  double tot_nb_st=0;
  if(nbe1!=0 && nbe2!=0 && nb_st_well[typ_ind]!=0){
    tot_nb_st=nbe1+nbe2;
    tot_nb_st-=nb_st_well[typ_ind];
#ifdef BI_ANNEALING
    tot_nb_st*=bia_lambda;
    if(tot_nb_st<bia_cap)
      tot_nb_st=bia_cap;
#endif
  }
  
  return tot_nb_st;
}




/* double MC_calc_nb_watscric_no_dih(int typ_ind, double *r_vec, double *r_vec_inv, int *flag){ */
/*   double nbe=0, nbei=0; */
/*   /\* if(r_vec[2]>0){//CASE 1 *\/ */
/*   if(table_nb_N_2[typ_ind][0]>0){ */
/*     nbe=calc_nb_tab_watscric(r_vec, typ_ind); // watson-crick */
/*     if(nbe == 0){ */
/*       *flag=1; */
/*       //printf("forbidden wc! \t %lf %lf %lf     %lf %lf %lf\n", r_vec[0], r_vec[1], r_vec[2], r_vec_inv[0], r_vec_inv[1], r_vec_inv[2]); */
/*     } */
/*   } */
/*   if(table_nb_N_2_inv[typ_ind][0]>0){ */
/*     nbei=calc_nb_tab_watscric_inv(r_vec_inv, typ_ind); // watson-crick */
/*     if(nbei == 0){ */
/*       *flag=1; */
/*       //printf("forbidden wc_inv! \t %lf %lf %lf     %lf %lf %lf\n", r_vec[0], r_vec[1], r_vec[2], r_vec_inv[0], r_vec_inv[1], r_vec_inv[2]); */
/*     } */
/*   } */
/*   return (nbe+nbei); */
/* } */






//int stacking_happens(double dist, double *r_vec, double eta){
//  if(dist<5) return 1;
//  else return 0;
//}




double calc_st_psdihedr(double *rx1, double *ry1, double *rz1, double *rx2, double *ry2, double *rz2, int nt_c, int nt_neigh){
  double cent1[DIM], cent2[DIM];
  double x1[DIM], x2[DIM], y1[DIM], y2[DIM];
  double temph1[DIM], temph2[DIM], temph3[DIM], d_vec2[DIM];
  double cos_dih_st, sin_dih_st;
  int at_c=N_PARTS_PER_NT*nt_c;
  x1[0]=dist_1d(rx1[at_c+1],rx1[at_c],0);  x1[1]=dist_1d(ry1[at_c+1],ry1[at_c],1);  x1[2]=dist_1d(rz1[at_c+1],rz1[at_c],2);
  y1[0]=dist_1d(rx1[at_c+2],rx1[at_c],0);  y1[1]=dist_1d(ry1[at_c+2],ry1[at_c],1);  y1[2]=dist_1d(rz1[at_c+2],rz1[at_c],2);
  
  
  int np=nt_neigh*N_PARTS_PER_NT;
  x2[0]=dist_1d(rx2[np+1],rx2[np],0);  x2[1]=dist_1d(ry2[np+1],ry2[np],1);  x2[2]=dist_1d(rz2[np+1],rz2[np],2);
  y2[0]=dist_1d(rx2[np+2],rx2[np],0);  y2[1]=dist_1d(ry2[np+2],ry2[np],1);  y2[2]=dist_1d(rz2[np+2],rz2[np],2);
  cent1[0]=-(x1[0]+y1[0])/sqrt(2.0);  cent1[1]=-(x1[1]+y1[1])/sqrt(2.0);  cent1[2]=-(x1[2]+y1[2])/sqrt(2.0);
  cent2[0]=(x2[0]+y2[0])/sqrt(2.0);  cent2[1]=(x2[1]+y2[1])/sqrt(2.0);  cent2[2]=(x2[2]+y2[2])/sqrt(2.0);
  d_vec2[0]=dist_1d(rx2[np],rx1[at_c],0);d_vec2[1]=dist_1d(ry2[np],ry1[at_c],1);d_vec2[2]=dist_1d(rz2[np],rz1[at_c],2);
  
  vec_prod(cent1, d_vec2, temph1);
  vec_prod(d_vec2, cent2, temph2);
  cos_dih_st=dot_prod(temph1, temph2)/sqrt(dot_prod(temph1, temph1)*dot_prod(temph2, temph2));
  vec_prod(temph1, temph2, temph3);
  sin_dih_st=dot_prod(temph3, d_vec2)/sqrt(dot_prod(d_vec2, d_vec2)*dot_prod(temph1, temph1)*dot_prod(temph2,temph2));
  //printf("%lf\n", atan2(sin_dih_st, cos_dih_st));
  return atan2(sin_dih_st, cos_dih_st);
}



/* double calc_ssB0_sugars(double dist){ */
  
/*   int dist_in=(int)((dist-table_ssB0_params_0[1]+0.5*table_ssB0_params_0[0])/table_ssB0_params_0[0]); */
/*   //printf("%lf %d %d\n", dist, dist_in, table_ssB_N_0); */
/*   if(dist_in>=table_ssB0_N_0 || dist_in < 0) return MC_CAP+1.0; */
/*   return table_ssB0_0[dist_in]; */
/* } */


double calc_ssB1_sugars(double angle, int type){
  int ang_in;
  double ret=0;
  switch (type){
  case PUCKS_33:
    ang_in=(int)((angle-table_ssB1_params_33[1]+0.5*table_ssB1_params_33[0])/table_ssB1_params_33[0]);
    if(ang_in>=table_ssB1_N_33 || ang_in < 0) return MC_CAP+1.0;
    ret= table_ssB1_33[ang_in];
    break;
case PUCKS_23:
    ang_in=(int)((angle-table_ssB1_params_23[1]+0.5*table_ssB1_params_23[0])/table_ssB1_params_23[0]);
    if(ang_in>=table_ssB1_N_23 || ang_in < 0) return MC_CAP+1.0;
    ret= table_ssB1_23[ang_in];
    break;
case PUCKS_32:
    ang_in=(int)((angle-table_ssB1_params_32[1]+0.5*table_ssB1_params_32[0])/table_ssB1_params_32[0]);
    if(ang_in>=table_ssB1_N_32 || ang_in < 0) return MC_CAP+1.0;
    ret= table_ssB1_32[ang_in];
    break;
case PUCKS_22:
    ang_in=(int)((angle-table_ssB1_params_22[1]+0.5*table_ssB1_params_22[0])/table_ssB1_params_22[0]);
    if(ang_in>=table_ssB1_N_22 || ang_in < 0) return MC_CAP+1.0;
    ret= table_ssB1_22[ang_in];
    break;
  
  default: ret= 0;
    break;
  }
  
  return BB_PREF_A*ret;
  //return ret;
}

double calc_nnB_tab_st_psdihedr(double eta, int typ_ind, int stfaces){
  int eta_in;
  double RET=0;
  switch(stfaces){
  case STFACES33:
    eta_in=(int)((eta-table_nnB_params_1_s33[typ_ind][1]+0.5*table_nnB_params_1_s33[typ_ind][0])/table_nnB_params_1_s33[typ_ind][0]);
    if(eta_in>=table_nnB_N_1_s33[typ_ind] || eta_in<0) return 0;
    RET=table_nnB_1_s33[typ_ind][eta_in];
    break;
  case STFACES35:
    eta_in=(int)((eta-table_nnB_params_1_s35[typ_ind][1]+0.5*table_nnB_params_1_s35[typ_ind][0])/table_nnB_params_1_s35[typ_ind][0]);
    if(eta_in>=table_nnB_N_1_s35[typ_ind] || eta_in<0) return 0;
    RET=table_nnB_1_s35[typ_ind][eta_in];
    break;
  case STFACES53:
    eta_in=(int)((eta-table_nnB_params_1_s53[typ_ind][1]+0.5*table_nnB_params_1_s53[typ_ind][0])/table_nnB_params_1_s53[typ_ind][0]);
    if(eta_in>=table_nnB_N_1_s53[typ_ind] || eta_in<0) return 0;
    RET=table_nnB_1_s53[typ_ind][eta_in];
    break;
  case STFACES55:
    eta_in=(int)((eta-table_nnB_params_1_s55[typ_ind][1]+0.5*table_nnB_params_1_s55[typ_ind][0])/table_nnB_params_1_s55[typ_ind][0]);
    if(eta_in>=table_nnB_N_1_s55[typ_ind] || eta_in<0) return 0;
    RET=table_nnB_1_s55[typ_ind][eta_in];
    break;
  default:
    RET=0;
    printf("Unknown combination of stacking faces for dihedral energy calculation.\n");
    exit(1);
    break;
  }
  return RET;
}

//self_e=MC_calc_SP_intra(mc_types[at_c], r_vec, mc_glyc[nt_c], flag);
double MC_calc_BP_intra(int typ_ind, double *r_vec, int glyc, int puck, int *flag){
  int xin, yin, zin, pos_ind;
  int flag2=0;
  double ret=0, energ_ret;
  switch(glyc){
  case GLYC_A:
    switch(puck){
    case PUCK_3:
      xin=(int)((r_vec[0]-table_bpI_params_A3[typ_ind][0][1]+0.5*table_bpI_params_A3[typ_ind][0][0])/table_bpI_params_A3[typ_ind][0][0]);
      yin=(int)((r_vec[1]-table_bpI_params_A3[typ_ind][1][1]+0.5*table_bpI_params_A3[typ_ind][1][0])/table_bpI_params_A3[typ_ind][1][0]);
      zin=(int)((r_vec[2]-table_bpI_params_A3[typ_ind][2][1]+0.5*table_bpI_params_A3[typ_ind][2][0])/table_bpI_params_A3[typ_ind][2][0]);
      pos_ind = zin+table_bpI_N_A3[typ_ind][2]*yin+table_bpI_N_A3[typ_ind][1]*table_bpI_N_A3[typ_ind][2]*xin;
      if(xin<0 || xin>=table_bpI_N_A3[typ_ind][0]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(yin<0 || yin>=table_bpI_N_A3[typ_ind][1]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(zin<0 || zin>=table_bpI_N_A3[typ_ind][2]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(flag2==0){
	if(table_bpI_A3[typ_ind][pos_ind]==0)
	  *flag=1;
	ret=table_bpI_A3[typ_ind][pos_ind];
      }
      break;
    case PUCK_2:
      xin=(int)((r_vec[0]-table_bpI_params_A2[typ_ind][0][1]+0.5*table_bpI_params_A2[typ_ind][0][0])/table_bpI_params_A2[typ_ind][0][0]);
      yin=(int)((r_vec[1]-table_bpI_params_A2[typ_ind][1][1]+0.5*table_bpI_params_A2[typ_ind][1][0])/table_bpI_params_A2[typ_ind][1][0]);
      zin=(int)((r_vec[2]-table_bpI_params_A2[typ_ind][2][1]+0.5*table_bpI_params_A2[typ_ind][2][0])/table_bpI_params_A2[typ_ind][2][0]);
      pos_ind = zin+table_bpI_N_A2[typ_ind][2]*yin+table_bpI_N_A2[typ_ind][1]*table_bpI_N_A2[typ_ind][2]*xin;
      if(xin<0 || xin>=table_bpI_N_A2[typ_ind][0]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(yin<0 || yin>=table_bpI_N_A2[typ_ind][1]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(zin<0 || zin>=table_bpI_N_A2[typ_ind][2]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(flag2==0){
	if(table_bpI_A2[typ_ind][pos_ind]==0)
	  *flag=1;
	ret=table_bpI_A2[typ_ind][pos_ind];
      }
      break;
    }    
    break;
  case GLYC_H:
    switch(puck){
    case PUCK_3:
      xin=(int)((r_vec[0]-table_bpI_params_H3[typ_ind][0][1]+0.5*table_bpI_params_H3[typ_ind][0][0])/table_bpI_params_H3[typ_ind][0][0]);
      yin=(int)((r_vec[1]-table_bpI_params_H3[typ_ind][1][1]+0.5*table_bpI_params_H3[typ_ind][1][0])/table_bpI_params_H3[typ_ind][1][0]);
      zin=(int)((r_vec[2]-table_bpI_params_H3[typ_ind][2][1]+0.5*table_bpI_params_H3[typ_ind][2][0])/table_bpI_params_H3[typ_ind][2][0]);
      pos_ind = zin+table_bpI_N_H3[typ_ind][2]*yin+table_bpI_N_H3[typ_ind][1]*table_bpI_N_H3[typ_ind][2]*xin;
      if(xin<0 || xin>=table_bpI_N_H3[typ_ind][0]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(yin<0 || yin>=table_bpI_N_H3[typ_ind][1]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(zin<0 || zin>=table_bpI_N_H3[typ_ind][2]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(flag2==0){
	if(table_bpI_H3[typ_ind][pos_ind]==0)
	  *flag=1;
	ret=table_bpI_H3[typ_ind][pos_ind];
      }
      break;
    case PUCK_2:
      xin=(int)((r_vec[0]-table_bpI_params_H2[typ_ind][0][1]+0.5*table_bpI_params_H2[typ_ind][0][0])/table_bpI_params_H2[typ_ind][0][0]);
      yin=(int)((r_vec[1]-table_bpI_params_H2[typ_ind][1][1]+0.5*table_bpI_params_H2[typ_ind][1][0])/table_bpI_params_H2[typ_ind][1][0]);
      zin=(int)((r_vec[2]-table_bpI_params_H2[typ_ind][2][1]+0.5*table_bpI_params_H2[typ_ind][2][0])/table_bpI_params_H2[typ_ind][2][0]);
      pos_ind = zin+table_bpI_N_H2[typ_ind][2]*yin+table_bpI_N_H2[typ_ind][1]*table_bpI_N_H2[typ_ind][2]*xin;
      if(xin<0 || xin>=table_bpI_N_H2[typ_ind][0]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(yin<0 || yin>=table_bpI_N_H2[typ_ind][1]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(zin<0 || zin>=table_bpI_N_H2[typ_ind][2]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(flag2==0){
	if(table_bpI_H2[typ_ind][pos_ind]==0)
	  *flag=1;
	ret=table_bpI_H2[typ_ind][pos_ind];
      }
      break;
    }
    break;
  case GLYC_S:
    switch(puck){
    case PUCK_3:
      xin=(int)((r_vec[0]-table_bpI_params_S3[typ_ind][0][1]+0.5*table_bpI_params_S3[typ_ind][0][0])/table_bpI_params_S3[typ_ind][0][0]);
      yin=(int)((r_vec[1]-table_bpI_params_S3[typ_ind][1][1]+0.5*table_bpI_params_S3[typ_ind][1][0])/table_bpI_params_S3[typ_ind][1][0]);
      zin=(int)((r_vec[2]-table_bpI_params_S3[typ_ind][2][1]+0.5*table_bpI_params_S3[typ_ind][2][0])/table_bpI_params_S3[typ_ind][2][0]);
      pos_ind = zin+table_bpI_N_S3[typ_ind][2]*yin+table_bpI_N_S3[typ_ind][1]*table_bpI_N_S3[typ_ind][2]*xin;
      if(xin<0 || xin>=table_bpI_N_S3[typ_ind][0]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(yin<0 || yin>=table_bpI_N_S3[typ_ind][1]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(zin<0 || zin>=table_bpI_N_S3[typ_ind][2]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(flag2==0){
	if(table_bpI_S3[typ_ind][pos_ind]==0)
	  *flag=1;
	ret=table_bpI_S3[typ_ind][pos_ind];
      }
      break;
    case PUCK_2:
      xin=(int)((r_vec[0]-table_bpI_params_S2[typ_ind][0][1]+0.5*table_bpI_params_S2[typ_ind][0][0])/table_bpI_params_S2[typ_ind][0][0]);
      yin=(int)((r_vec[1]-table_bpI_params_S2[typ_ind][1][1]+0.5*table_bpI_params_S2[typ_ind][1][0])/table_bpI_params_S2[typ_ind][1][0]);
      zin=(int)((r_vec[2]-table_bpI_params_S2[typ_ind][2][1]+0.5*table_bpI_params_S2[typ_ind][2][0])/table_bpI_params_S2[typ_ind][2][0]);
      pos_ind = zin+table_bpI_N_S2[typ_ind][2]*yin+table_bpI_N_S2[typ_ind][1]*table_bpI_N_S2[typ_ind][2]*xin;
      if(xin<0 || xin>=table_bpI_N_S2[typ_ind][0]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(yin<0 || yin>=table_bpI_N_S2[typ_ind][1]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(zin<0 || zin>=table_bpI_N_S2[typ_ind][2]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(flag2==0){
	if(table_bpI_S2[typ_ind][pos_ind]==0)
	  *flag=1;
	ret=table_bpI_S2[typ_ind][pos_ind];
      }
      break;
    }
    break;
  default:
    *flag=1;
    ret= 0;
    break;
  }
  energ_ret=BB_PREF*ret;
#ifdef WARMUP
  double phpos[DIM];
  if(*flag==1){
    *flag=0;
    phpos[0]=ph_pintra_ave[typ_ind][glyc][puck][0];phpos[1]=ph_pintra_ave[typ_ind][glyc][puck][1];phpos[2]=ph_pintra_ave[typ_ind][glyc][puck][2];
    energ_ret=MC_CAP*((r_vec[0]-phpos[0])*(r_vec[0]-phpos[0])+(r_vec[1]-phpos[1])*(r_vec[1]-phpos[1])+(r_vec[2]-phpos[2])*(r_vec[2]-phpos[2]));
    //printf("   %lf %lf %lf    %lf %lf %lf     %d  %d  %d\n", r_vec[0], r_vec[1], r_vec[2], phpos[0], phpos[1], phpos[2], typ_ind, glyc, puck);
    //energ_ret=MC_CAP-1;
    //*(r_vec[0]*r_vec[0]+r_vec[1]*r_vec[1]+r_vec[2]*r_vec[2]);
  }
#endif

  return energ_ret; 
}

double MC_calc_BP_inter(int typ_ind, double *r_vec, int glyc, int puck, int *flag){
  int xin, yin, zin, pos_ind;
  int flag2=0;
  double ret=0,energ_ret;
  //typ_ind = N_BASES*mc_types[at_c]    + mc_types[at_ne];//N_BASES*type1+type2;
  switch(glyc){
  case GLYC_A:
    switch(puck){
    case PUCK_3:
      xin=(int)((r_vec[0]-table_bpB_params_A3[typ_ind][0][1]+0.5*table_bpB_params_A3[typ_ind][0][0])/table_bpB_params_A3[typ_ind][0][0]);
      yin=(int)((r_vec[1]-table_bpB_params_A3[typ_ind][1][1]+0.5*table_bpB_params_A3[typ_ind][1][0])/table_bpB_params_A3[typ_ind][1][0]);
      zin=(int)((r_vec[2]-table_bpB_params_A3[typ_ind][2][1]+0.5*table_bpB_params_A3[typ_ind][2][0])/table_bpB_params_A3[typ_ind][2][0]);
      pos_ind = zin+table_bpB_N_A3[typ_ind][2]*yin+table_bpB_N_A3[typ_ind][1]*table_bpB_N_A3[typ_ind][2]*xin;
      if(xin<0 || xin>=table_bpB_N_A3[typ_ind][0]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(yin<0 || yin>=table_bpB_N_A3[typ_ind][1]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(zin<0 || zin>=table_bpB_N_A3[typ_ind][2]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(flag2==0){
	if(table_bpB_A3[typ_ind][pos_ind]==0)
	  *flag=1;
	ret= table_bpB_A3[typ_ind][pos_ind];
	
	
      }
      break;
    case PUCK_2:
      xin=(int)((r_vec[0]-table_bpB_params_A2[typ_ind][0][1]+0.5*table_bpB_params_A2[typ_ind][0][0])/table_bpB_params_A2[typ_ind][0][0]);
      yin=(int)((r_vec[1]-table_bpB_params_A2[typ_ind][1][1]+0.5*table_bpB_params_A2[typ_ind][1][0])/table_bpB_params_A2[typ_ind][1][0]);
      zin=(int)((r_vec[2]-table_bpB_params_A2[typ_ind][2][1]+0.5*table_bpB_params_A2[typ_ind][2][0])/table_bpB_params_A2[typ_ind][2][0]);
      pos_ind = zin+table_bpB_N_A2[typ_ind][2]*yin+table_bpB_N_A2[typ_ind][1]*table_bpB_N_A2[typ_ind][2]*xin;
      if(xin<0 || xin>=table_bpB_N_A2[typ_ind][0]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(yin<0 || yin>=table_bpB_N_A2[typ_ind][1]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(zin<0 || zin>=table_bpB_N_A2[typ_ind][2]) {
	*flag=1;
	flag2=1;
	ret= 0;}
      if(flag2==0){
	if(table_bpB_A2[typ_ind][pos_ind]==0)
	  *flag=1;
	ret= table_bpB_A2[typ_ind][pos_ind];
      }
      break;
    }
    break;
  case GLYC_H:
    switch(puck){
    case PUCK_3:
    xin=(int)((r_vec[0]-table_bpB_params_H3[typ_ind][0][1]+0.5*table_bpB_params_H3[typ_ind][0][0])/table_bpB_params_H3[typ_ind][0][0]);
    yin=(int)((r_vec[1]-table_bpB_params_H3[typ_ind][1][1]+0.5*table_bpB_params_H3[typ_ind][1][0])/table_bpB_params_H3[typ_ind][1][0]);
    zin=(int)((r_vec[2]-table_bpB_params_H3[typ_ind][2][1]+0.5*table_bpB_params_H3[typ_ind][2][0])/table_bpB_params_H3[typ_ind][2][0]);
    pos_ind = zin+table_bpB_N_H3[typ_ind][2]*yin+table_bpB_N_H3[typ_ind][1]*table_bpB_N_H3[typ_ind][2]*xin;
    if(xin<0 || xin>=table_bpB_N_H3[typ_ind][0]) {
      *flag=1;
      flag2=1;
      ret= 0;}
    if(yin<0 || yin>=table_bpB_N_H3[typ_ind][1]) {
      *flag=1;
      flag2=1;
      ret= 0;}
    if(zin<0 || zin>=table_bpB_N_H3[typ_ind][2]) {
      *flag=1;
      flag2=1;
      ret= 0;}
    if(flag2==0){
      if(table_bpB_H3[typ_ind][pos_ind]==0)
	*flag=1;
      ret= table_bpB_H3[typ_ind][pos_ind];
    }
    break;
    case PUCK_2:
    xin=(int)((r_vec[0]-table_bpB_params_H2[typ_ind][0][1]+0.5*table_bpB_params_H2[typ_ind][0][0])/table_bpB_params_H2[typ_ind][0][0]);
    yin=(int)((r_vec[1]-table_bpB_params_H2[typ_ind][1][1]+0.5*table_bpB_params_H2[typ_ind][1][0])/table_bpB_params_H2[typ_ind][1][0]);
    zin=(int)((r_vec[2]-table_bpB_params_H2[typ_ind][2][1]+0.5*table_bpB_params_H2[typ_ind][2][0])/table_bpB_params_H2[typ_ind][2][0]);
    pos_ind = zin+table_bpB_N_H2[typ_ind][2]*yin+table_bpB_N_H2[typ_ind][1]*table_bpB_N_H2[typ_ind][2]*xin;
    if(xin<0 || xin>=table_bpB_N_H2[typ_ind][0]) {
      *flag=1;
      flag2=1;
      ret= 0;}
    if(yin<0 || yin>=table_bpB_N_H2[typ_ind][1]) {
      *flag=1;
      flag2=1;
      ret= 0;}
    if(zin<0 || zin>=table_bpB_N_H2[typ_ind][2]) {
      *flag=1;
      flag2=1;
      ret= 0;}
    if(flag2==0){
      if(table_bpB_H2[typ_ind][pos_ind]==0)
	*flag=1;
      ret= table_bpB_H2[typ_ind][pos_ind];
    }
    break;
    }
    break;
  case GLYC_S:
    switch(puck){
    case PUCK_3:
    xin=(int)((r_vec[0]-table_bpB_params_S3[typ_ind][0][1]+0.5*table_bpB_params_S3[typ_ind][0][0])/table_bpB_params_S3[typ_ind][0][0]);
    yin=(int)((r_vec[1]-table_bpB_params_S3[typ_ind][1][1]+0.5*table_bpB_params_S3[typ_ind][1][0])/table_bpB_params_S3[typ_ind][1][0]);
    zin=(int)((r_vec[2]-table_bpB_params_S3[typ_ind][2][1]+0.5*table_bpB_params_S3[typ_ind][2][0])/table_bpB_params_S3[typ_ind][2][0]);
    pos_ind = zin+table_bpB_N_S3[typ_ind][2]*yin+table_bpB_N_S3[typ_ind][1]*table_bpB_N_S3[typ_ind][2]*xin;
    if(xin<0 || xin>=table_bpB_N_S3[typ_ind][0]) {
      *flag=1;
      flag2=1;
      ret= 0;}
    if(yin<0 || yin>=table_bpB_N_S3[typ_ind][1]) {
      *flag=1;
      flag2=1;
      ret= 0;}
    if(zin<0 || zin>=table_bpB_N_S3[typ_ind][2]) {
      *flag=1;
      flag2=1;
      ret= 0;}
    if(flag2==0){
      if(table_bpB_S3[typ_ind][pos_ind]==0)
	*flag=1;
      ret= table_bpB_S3[typ_ind][pos_ind];
    }
    break;
    case PUCK_2:
      xin=(int)((r_vec[0]-table_bpB_params_S2[typ_ind][0][1]+0.5*table_bpB_params_S2[typ_ind][0][0])/table_bpB_params_S2[typ_ind][0][0]);
      yin=(int)((r_vec[1]-table_bpB_params_S2[typ_ind][1][1]+0.5*table_bpB_params_S2[typ_ind][1][0])/table_bpB_params_S2[typ_ind][1][0]);
      zin=(int)((r_vec[2]-table_bpB_params_S2[typ_ind][2][1]+0.5*table_bpB_params_S2[typ_ind][2][0])/table_bpB_params_S2[typ_ind][2][0]);
    pos_ind = zin+table_bpB_N_S2[typ_ind][2]*yin+table_bpB_N_S2[typ_ind][1]*table_bpB_N_S2[typ_ind][2]*xin;
    if(xin<0 || xin>=table_bpB_N_S2[typ_ind][0]) {
      *flag=1;
      flag2=1;
      ret= 0;}
    if(yin<0 || yin>=table_bpB_N_S2[typ_ind][1]) {
      *flag=1;
      flag2=1;
      ret= 0;}
    if(zin<0 || zin>=table_bpB_N_S2[typ_ind][2]) {
      *flag=1;
      flag2=1;
      ret= 0;}
    if(flag2==0){
      if(table_bpB_S2[typ_ind][pos_ind]==0)
	*flag=1;
      ret= table_bpB_S2[typ_ind][pos_ind];
    }
      
    break;
    }
    break;
    
  default: 
    *flag=1;
    ret= 0;
    break;
  }
  energ_ret=BB_PREF*ret;
  return energ_ret;
}

double calc_nnN_tab_stacking_s3(double *r_vec, int typ_ind){
  int xin, yin, zin;
  xin=(int)((r_vec[0]-table_nnN_params_0s3[typ_ind][0][1]+0.5*table_nnN_params_0s3[typ_ind][0][0])/table_nnN_params_0s3[typ_ind][0][0]);
  yin=(int)((r_vec[1]-table_nnN_params_0s3[typ_ind][1][1]+0.5*table_nnN_params_0s3[typ_ind][1][0])/table_nnN_params_0s3[typ_ind][1][0]);
  zin=(int)((r_vec[2]-table_nnN_params_0s3[typ_ind][2][1]+0.5*table_nnN_params_0s3[typ_ind][2][0])/table_nnN_params_0s3[typ_ind][2][0]);
  int pos_ind = zin+table_nnN_N_0s3[typ_ind][2]*yin+table_nnN_N_0s3[typ_ind][1]*table_nnN_N_0s3[typ_ind][2]*xin;
  if(xin<0 || xin>=table_nnN_N_0s3[typ_ind][0])
    return 0;
  if(yin<0 || yin>=table_nnN_N_0s3[typ_ind][1])
    return 0;
  if(zin<0 || zin>=table_nnN_N_0s3[typ_ind][2])
    return 0;
  return table_nnN_0s3[typ_ind][pos_ind];
}
 
double calc_nnN_tab_stacking_s5(double *r_vec, int typ_ind){
  int xin, yin, zin;
  xin=(int)((r_vec[0]-table_nnN_params_0s5[typ_ind][0][1]+0.5*table_nnN_params_0s5[typ_ind][0][0])/table_nnN_params_0s5[typ_ind][0][0]);
  yin=(int)((r_vec[1]-table_nnN_params_0s5[typ_ind][1][1]+0.5*table_nnN_params_0s5[typ_ind][1][0])/table_nnN_params_0s5[typ_ind][1][0]);
  zin=(int)((r_vec[2]-table_nnN_params_0s5[typ_ind][2][1]+0.5*table_nnN_params_0s5[typ_ind][2][0])/table_nnN_params_0s5[typ_ind][2][0]);
  int pos_ind = zin+table_nnN_N_0s5[typ_ind][2]*yin+table_nnN_N_0s5[typ_ind][1]*table_nnN_N_0s5[typ_ind][2]*xin;
  if(xin<0 || xin>=table_nnN_N_0s5[typ_ind][0])
    return 0;
  if(yin<0 || yin>=table_nnN_N_0s5[typ_ind][1])
    return 0;
  if(zin<0 || zin>=table_nnN_N_0s5[typ_ind][2])
    return 0;
  return table_nnN_0s5[typ_ind][pos_ind];
}

double calc_nnB_tab_stacking(double *r_vec, int typ_ind, int stind){
  int xin, yin, zin, nxmax, nymax, nzmax;
  double xmin, ymin, zmin, dx, dy, dz;
  double RET;
  switch(stind){
  case STFACES33:
    xmin=table_nnB_params_0_s33[typ_ind][0][1];ymin=table_nnB_params_0_s33[typ_ind][1][1];zmin=table_nnB_params_0_s33[typ_ind][2][1];
    dx=table_nnB_params_0_s33[typ_ind][0][0];dy=table_nnB_params_0_s33[typ_ind][1][0];dz=table_nnB_params_0_s33[typ_ind][2][0];
    nxmax=table_nnB_N_0_s33[typ_ind][0];nymax=table_nnB_N_0_s33[typ_ind][1];nzmax=table_nnB_N_0_s33[typ_ind][2];
    break;
  case STFACES35:
    xmin=table_nnB_params_0_s35[typ_ind][0][1];ymin=table_nnB_params_0_s35[typ_ind][1][1];zmin=table_nnB_params_0_s35[typ_ind][2][1];
    dx=table_nnB_params_0_s35[typ_ind][0][0];dy=table_nnB_params_0_s35[typ_ind][1][0];dz=table_nnB_params_0_s35[typ_ind][2][0];
    nxmax=table_nnB_N_0_s35[typ_ind][0];nymax=table_nnB_N_0_s35[typ_ind][1];nzmax=table_nnB_N_0_s35[typ_ind][2];
    break;
  case STFACES53:
    xmin=table_nnB_params_0_s53[typ_ind][0][1];ymin=table_nnB_params_0_s53[typ_ind][1][1];zmin=table_nnB_params_0_s53[typ_ind][2][1];
    dx=table_nnB_params_0_s53[typ_ind][0][0];dy=table_nnB_params_0_s53[typ_ind][1][0];dz=table_nnB_params_0_s53[typ_ind][2][0];
    nxmax=table_nnB_N_0_s53[typ_ind][0];nymax=table_nnB_N_0_s53[typ_ind][1];nzmax=table_nnB_N_0_s53[typ_ind][2];
    break;
  case STFACES55:
    xmin=table_nnB_params_0_s55[typ_ind][0][1];ymin=table_nnB_params_0_s55[typ_ind][1][1];zmin=table_nnB_params_0_s55[typ_ind][2][1];
    dx=table_nnB_params_0_s55[typ_ind][0][0];dy=table_nnB_params_0_s55[typ_ind][1][0];dz=table_nnB_params_0_s55[typ_ind][2][0];
    nxmax=table_nnB_N_0_s55[typ_ind][0];nymax=table_nnB_N_0_s55[typ_ind][1];nzmax=table_nnB_N_0_s55[typ_ind][2];
    break;
  default:
    printf("Stacking face combination not defined!!\n"); exit(1);
    break;
  }
  
  xin=(int)((r_vec[0]-xmin+dx*0.5)/dx);
  yin=(int)((r_vec[1]-ymin+dy*0.5)/dy);
  zin=(int)((r_vec[2]-zmin+dz*0.5)/dz);
  if(xin<0 || xin >=nxmax) return 0;
  if(yin<0 || yin >=nymax) return 0;
  if(zin<0 || zin >=nzmax) return 0;
  //int pos_ind = zin+table_nnB_N_0[typ_ind][2]*yin+table_nnB_N_0[typ_ind][1]*table_nnB_N_0[typ_ind][2]*xin;
  int pos_ind = zin+nzmax*yin+nymax*nzmax*xin;

  switch(stind){
  case STFACES33:
    RET=table_nnB_0_s33[typ_ind][pos_ind];
    break;
  case STFACES35:
    RET=table_nnB_0_s35[typ_ind][pos_ind];
    break;
  case STFACES53:
    RET=table_nnB_0_s53[typ_ind][pos_ind];
    break;
  case STFACES55:
    RET=table_nnB_0_s55[typ_ind][pos_ind];
    break;
  default:
    RET=0;
    break;
  }
  return RET;
}
/* double calc_nnB_tab_stacking(double *r_vec, int typ_ind){ */
/*   int xin, yin, zin; */
/*   xin=(int)((r_vec[0]-table_nnB_params_0[typ_ind][0][1]+0.5*table_nnB_params_0[typ_ind][0][0])/table_nnB_params_0[typ_ind][0][0]); */
/*   yin=(int)((r_vec[1]-table_nnB_params_0[typ_ind][1][1]+0.5*table_nnB_params_0[typ_ind][1][0])/table_nnB_params_0[typ_ind][1][0]); */
/*   zin=(int)((r_vec[2]-table_nnB_params_0[typ_ind][2][1]+0.5*table_nnB_params_0[typ_ind][2][0])/table_nnB_params_0[typ_ind][2][0]); */
/*   int pos_ind = zin+table_nnB_N_0[typ_ind][2]*yin+table_nnB_N_0[typ_ind][1]*table_nnB_N_0[typ_ind][2]*xin; */
/*   if(xin<0 || xin>=table_nnB_N_0[typ_ind][0]) { */
/*     return 0;} */
/*   if(yin<0 || yin>=table_nnB_N_0[typ_ind][1]) { */
/*     return 0;} */
/*   if(zin<0 || zin>=table_nnB_N_0[typ_ind][2]) { */
/*     return 0;} */
/*   return table_nnB_0[typ_ind][pos_ind]; */
/* } */

/* double calc_nnB_tab_stacking_inv(double *r_vec, int typ_ind){ */
/*   int xin, yin, zin; */
/*   xin=(int)((r_vec[0]-table_nnB_params_0_inv[typ_ind][0][1]+0.5*table_nnB_params_0_inv[typ_ind][0][0])/table_nnB_params_0_inv[typ_ind][0][0]); */
/*   yin=(int)((r_vec[1]-table_nnB_params_0_inv[typ_ind][1][1]+0.5*table_nnB_params_0_inv[typ_ind][1][0])/table_nnB_params_0_inv[typ_ind][1][0]); */
/*   zin=(int)((r_vec[2]-table_nnB_params_0_inv[typ_ind][2][1]+0.5*table_nnB_params_0_inv[typ_ind][2][0])/table_nnB_params_0_inv[typ_ind][2][0]); */
/*   int pos_ind = zin+table_nnB_N_0_inv[typ_ind][2]*yin+table_nnB_N_0_inv[typ_ind][1]*table_nnB_N_0_inv[typ_ind][2]*xin; */
/*   if(xin<0 || xin>=table_nnB_N_0_inv[typ_ind][0]) { */
/*     return 0;} */
/*   if(yin<0 || yin>=table_nnB_N_0_inv[typ_ind][1]) { */
/*     return 0;} */
/*   if(zin<0 || zin>=table_nnB_N_0_inv[typ_ind][2]) { */
/*     return 0;} */
/*   return table_nnB_0_inv[typ_ind][pos_ind]; */
/* } */




double energy_excludedvol(double r, double sigma){
  /* if(r<0.88) return MC_CAP+1.0; */
  /* if(r<pow(2,1.0/6.0)*sigma){ */
  /*   double rm1=sigma/r; */
  /*   double rm6=rm1*rm1*rm1*rm1*rm1*rm1; */
  /*   /\* energy is 4*epsilon*(rm6*rm6-rm6+0.25) *\/ */
  /*   double e=4*epsilon/sigma*(rm6*rm6-rm6+0.25); */
  /*   return e; */
  /* } else return 0; */
  if(r<sigma) return MC_CAP+1.0;
  else return 0;
}


/* void calc_rel_pos_inv( double *d_vec , int nt_neigh, double * rx, double *ry, double *rz, double *r_vec){ */
/*   double pbase_x[DIM], pbase_y[DIM], pbase_z[DIM];//, pbasey[N_PARTS_PER_NT], pbasez[N_PARTS_PER_NT]; */
/*   //here we use that the first auxiliar particle is x, and the second, y */
/*   int nnt=N_PARTS_PER_NT*nt_neigh; */
/*   pbase_x[0]=dist_1d(rx[nnt+1],rx[nnt],0); */
/*   pbase_x[1]=dist_1d(ry[nnt+1],ry[nnt],0); */
/*   pbase_x[2]=dist_1d(rz[nnt+1],rz[nnt],0); */
/*   pbase_y[0]=dist_1d(rx[nnt+2],rx[nnt],0); */
/*   pbase_y[1]=dist_1d(ry[nnt+2],ry[nnt],0); */
/*   pbase_y[2]=dist_1d(rz[nnt+2],rz[nnt],0); */
/*   d_vec[0]*=-1;d_vec[1]*=-1;d_vec[2]*=-1; */
/*   vec_prod(pbase_x, pbase_y, pbase_z); */
/*   r_vec[0]=(d_vec[0]*pbase_x[0]+d_vec[1]*pbase_x[1]+d_vec[2]*pbase_x[2]); */
/*   r_vec[1]=(d_vec[0]*pbase_y[0]+d_vec[1]*pbase_y[1]+d_vec[2]*pbase_y[2]); */
/*   r_vec[2]=(d_vec[0]*pbase_z[0]+d_vec[1]*pbase_z[1]+d_vec[2]*pbase_z[2]); */
/*   //printf("%lf %lf %lf   %lf  %lf\n", r_vec[0], r_vec[1], r_vec[2], pbase_x[0]*pbase_x[0]+pbase_x[1]*pbase_x[1]+pbase_x[2]*pbase_x[2],pbase_y[0]*pbase_y[0]+pbase_y[1]*pbase_y[1]+pbase_y[2]*pbase_y[2]); */
/* } */
/* void calc_rel_pos( double *d_vec, double *r_vec){ */
/*   double pbase_x[DIM], pbase_y[DIM], pbase_z[DIM];//, pbasey[N_PARTS_PER_NT], pbasez[N_PARTS_PER_NT]; */
/*   //here we use that the first auxiliar particle is x, and the second, y */
/*   pbase_x[0]=dist_1d(mc_temp_x[N_PARTS_PER_NT*nt_c+1],mc_temp_x[N_PARTS_PER_NT*nt_c],0); */
/*   pbase_x[1]=dist_1d(mc_temp_y[N_PARTS_PER_NT*nt_c+1],mc_temp_y[N_PARTS_PER_NT*nt_c],0); */
/*   pbase_x[2]=dist_1d(mc_temp_z[N_PARTS_PER_NT*nt_c+1],mc_temp_z[N_PARTS_PER_NT*nt_c],0); */
/*   pbase_y[0]=dist_1d(mc_temp_x[N_PARTS_PER_NT*nt_c+2],mc_temp_x[N_PARTS_PER_NT*nt_c],0); */
/*   pbase_y[1]=dist_1d(mc_temp_y[N_PARTS_PER_NT*nt_c+2],mc_temp_y[N_PARTS_PER_NT*nt_c],0); */
/*   pbase_y[2]=dist_1d(mc_temp_z[N_PARTS_PER_NT*nt_c+2],mc_temp_z[N_PARTS_PER_NT*nt_c],0); */
  
/*   vec_prod(pbase_x, pbase_y, pbase_z); */
/*   r_vec[0]=(d_vec[0]*pbase_x[0]+d_vec[1]*pbase_x[1]+d_vec[2]*pbase_x[2]); */
/*   r_vec[1]=(d_vec[0]*pbase_y[0]+d_vec[1]*pbase_y[1]+d_vec[2]*pbase_y[2]); */
/*   r_vec[2]=(d_vec[0]*pbase_z[0]+d_vec[1]*pbase_z[1]+d_vec[2]*pbase_z[2]); */
/* } */




void MC_free_energy_params(){
  int i;
  for(i=0;i<N_MAX_TYPES;i++){
    free(mc_lj_sig[i]);
    free(mc_lj_eps[i]);
  }
  free(mc_lj_sig);
  free(mc_lj_eps);
  free(mc_harm_k);
  free(mc_harm_r);

  free(mc_ang_k);
#ifdef DIHEDRALS
  free(mc_dih_k);
#endif
  free(mc_glyc);
  free(mc_temp_glyc);
  
  
}
