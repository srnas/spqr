#include "mc_initialization.h"

int MC_detect_initial_condition(int mpi_id){
  FILE *infile;
  char initname[256];
  int flag=-1;
  //case 1
  printf("Detecting initial condition...");
  sprintf(initname, "pdb_inits/init.p%02d.mc", mpi_id);
  if((infile=fopen(initname, "rb"))==NULL){
    sprintf(initname, "pdb_inits/init.mc");
    if((infile=fopen(initname, "rb"))==NULL){
      sprintf(initname, "pdb_inits/init.p%02d.pdb", mpi_id);
      if((infile=fopen(initname, "r"))==NULL){
	sprintf(initname, "pdb_inits/init.pdb");
	if((infile=fopen(initname, "r"))==NULL){
	  printf("No valid initial condition found for processor %d!\n", mpi_id);
	  exit(ERR_INPUT);
	}
	else{
	  fclose(infile);
	  flag=3;
	}
      }
      else{
	fclose(infile);
	flag=2;
      }
    }
    else{
      fclose(infile);
      flag=1;
    }
  }
  else{
    fclose(infile);
    flag=0;
  }
  printf("Initial condition taken from %s.\n", initname);
  return flag;
}

int MC_initialize(int *mc_n, double **rsx, double **rsy, double **rsz, int *mc_iter, int *rand_a, int mpi_id, int init_type, int read_flag, double *add_data){
  int nt_n,at,nt;
  int init=0;
  char initfile[256];
  if(init_type==0 || init_type==1){
    if(init_type==0) sprintf(initfile, "pdb_inits/init.p%02d.mc", mpi_id);
    else  sprintf(initfile, "pdb_inits/init.mc"  );
    MC_read_params(mc_iter, rand_a,mpi_id);
   
    init=MC_read_checkpoint(mc_n, rsx, rsy, rsz, rand_a, mpi_id, initfile, read_flag, add_data);
  }
  else{
    if(init_type==2) sprintf(initfile, "pdb_inits/init.p%02d.pdb", mpi_id);
    else  sprintf(initfile, "pdb_inits/init.pdb"  );
    MC_read_nsolute(mc_n, mpi_id, initfile);
    MC_read_params(mc_iter, rand_a, mpi_id);
    MC_initialize_global(*mc_n, *rand_a, mpi_id);
    MC_initialize_arrays(*mc_n,rsx,rsy,rsz);
    MC_read_pdb(*mc_n, rsx, rsy, rsz, mpi_id,initfile);
  }
  /* INITIALIZE CONFIGURATION FILE */
  MC_initialize_save_configs(*mc_n, mpi_id);
#ifdef ERMSDR
  nt_n=*mc_n/N_PARTS_PER_NT;
  MC_init_ermsd_restr(nt_n);
  MC_init_ermsd_out(mpi_id);
  //ERMSD_SQ=get_first_ermsd(rsx, rsy, rsz, nt_n);
  get_first_ermsd(rsx, rsy, rsz, nt_n,&ERMSD_SQ,&ERMSD_ENERG);
  DELTA_ERMSD_SQ=0;
  DELTA_ERMSD_ENERG=0;
  
  WALL_ENERG=0;
  DELTA_WALL_ENERG=0;
  for(nt=0;nt<nt_n;nt++){
    WALL_ENERG+=MC_wall_energy((*rsx)[nt*N_PARTS_PER_NT+IPHO],(*rsy)[nt*N_PARTS_PER_NT+IPHO],(*rsz)[nt*N_PARTS_PER_NT+IPHO]);
  }
#endif
#ifdef LNKRMV
  nt_n=*mc_n/N_PARTS_PER_NT;
  MC_set_linked_loops(nt_n,*rsx,*rsy,*rsz);
#endif  

  MC_initialize_energies(*mc_n, *rsx, *rsy, *rsz);
  return init;
}

void MC_read_pdb(int mc_n, double **rx, double **ry, double **rz, int mpi_id, char *outpdbname){
  //read pdb - read types, positions, glp, frz, and bonds
  char *pdbrectyp, resname, *fullresname, *stmp, *sresind;
  char *atname, gl, pu, mv, fr, tgl, tpu, tmv, tfr;
  char pdbname[MAXSTR];
  static size_t st_l=0;
  
  int nt=0, nt_n=mc_n/N_PARTS_PER_NT, at=0, inat=0, l, i, j, atindex;
  char chain;
  int firstat, firstnt, isfirst=1;
  int tempres, prevres;
  int *chains; chains=(int *)malloc(sizeof(int)*nt_n);
  char tempchain, prevchain; int currchain=0;
  double tx, ty, tz;
  
  FILE *pdbfile;
  int fileflag=1, fileflaggen=1;
  char *lline=NULL;
  
  if(outpdbname==NULL)
    sprintf(pdbname,"pdb_inits/mc.p%02d.pdb", mpi_id); 
  else 
    strcpy(pdbname, outpdbname); 
  if((pdbfile=fopen(pdbname, "r"))==NULL){
    if(outpdbname!=NULL){
      printf("pdb file %s not found!\n", outpdbname);
      exit(ERR_INPUT);
    }
    fileflag=0;
    sprintf(pdbname,"pdb_inits/mc.pdb"); 
    if((pdbfile=fopen(pdbname, "r"))==NULL){
      printf("No decent init files found in pdb_inits. Can not run like this.\n");
      exit(ERR_INPUT);
      fileflaggen=0;
    }
  }

  while((l=getline(&lline, &st_l, pdbfile))>6){
    //for(j=0;j<6;j++)
    //pdbrectyp[j]=lline[j];
    pdbrectyp=strndup(lline, 6);
    if(!strcmp(pdbrectyp, "ATOM  ")){
      //read    
      //stmp=strndup(lline+6, 5);atindex=strtol(stmp, NULL, 10);free(stmp);
      stmp=strndup(lline+6, 5);atindex=atoi(stmp);free(stmp);
      atname=strndup(lline+12,4);//lline[15];
      fullresname=strndup(lline+17,3);
      if(fullresname[0]==' ' && fullresname[1]==' ') resname=fullresname[2];
      else if(fullresname[1]==' ' && fullresname[2]==' ') resname=fullresname[0];
      else if(fullresname[0]==' ' && fullresname[2]==' ') resname=fullresname[1];
      else{printf("Residue not recognized in atom %d\n", at); exit(ERR_INPUT);}
      //resname=lline[19];
      //stmp=strndup(lline+21, 1);chain=atoi(stmp);free(stmp);
      //chain=atoi(lline[21]);
      if(inat==0){
	sresind=strndup(lline+22,4);tempres=atoi(sresind);free(sresind);
	tempchain=lline[21];
	if(at==0) {prevchain=tempchain;prevres=tempres-1;}
	if((tempchain!=prevchain || tempres!=(prevres+1))) currchain++;
	//chain=currchain;
	//prevchain=currchain;prevres=tempres;
      }
      //printf("%d %d %d %d -  %c %d %c\n", at, nt, tempres, prevres, tempchain, currchain, prevchain);
      stmp=strndup(lline+30, 8);tx=atof(stmp);free(stmp);
      stmp=strndup(lline+38, 8);ty=atof(stmp);free(stmp);
      stmp=strndup(lline+46, 8);tz=atof(stmp);free(stmp);
      
      gl='A';
      pu='3';
      mv='N';
      fr='A';
      if(strlen(lline)>58 && inat==0){

	tgl=lline[56];
	tpu=lline[57];
	tmv=lline[58];
	tfr=lline[59];
	if((tgl=='A' || tgl=='H' || tgl=='S') && (tpu=='3' || tpu=='2') && (tmv=='N' || tmv=='G' || tmv=='P' || tmv=='A') && (tfr=='N' || tfr=='B' || tfr=='P' || tfr=='A')){
	  gl=tgl;pu=tpu;mv=tmv;fr=tfr;
	}
	//printf("READ PDB %d  %c %c %c %c\n", nt, gl, pu, mv, fr);
      }
      //check indexes
      if(isfirst==1){
	firstat=atindex;
	//firstnt=resindex;
	isfirst=0;
      }
      
      //set position
      (*rx)[at]=tx;
      (*ry)[at]=ty;
      (*rz)[at]=tz;
      if(strcmp(atname, "BASE")==0){
	if(resname=='A') mc_types[at]=TYP_ADENINE;
	if(resname=='U') mc_types[at]=TYP_URACIL;
	if(resname=='G') mc_types[at]=TYP_GUANINE;
	if(resname=='C') mc_types[at]=TYP_CYTOSINE;
      }
      if(strcmp(atname, "XVEC")==0) mc_types[at]=-1;
      if(strcmp(atname, "YVEC")==0) mc_types[at]=-2;
      if(strcmp(atname, "SUGR")==0) mc_types[at]=-3;
      if(strcmp(atname, "PHOS")==0) mc_types[at]=-4;
      
      /*** THIS SETUP IS ONLY READ ON BASES ***/
      if(inat==0){
	chains[nt]=currchain;
	prevchain=tempchain; prevres=tempres;
	//set glp
	if(gl=='A') mc_glyc[nt]=GLYC_A;
	else if(gl=='H') mc_glyc[nt]=GLYC_H;
	else if(gl=='S') mc_glyc[nt]=GLYC_S;
	else {
	  printf("Invalid GLYCOSIDIC index for %d.\n", nt);
	  exit(ERR_INPUT);
	}
	
	if(pu=='3') mc_puck[nt]=PUCK_3;
	else if(pu=='2') mc_puck[nt]=PUCK_2;
	else {
	  printf("Invalid PUCKER index for %d in gly_pck.dat file.\n", nt);
	  exit(ERR_INPUT);
	}
	
	if(mv=='N') glp_is_flippable[nt]=GLP_FIXED;
	else if(mv=='G') glp_is_flippable[nt]=GLP_GLYC;
	else if(mv=='P') glp_is_flippable[nt]=GLP_PUCK;
	else if(mv=='A') glp_is_flippable[nt]=GLP_BOTH;
	else{
	  printf("Invalid fix flag in glyc and puck at nucleoside %d .\n", nt);
	  exit(ERR_INPUT);
	}
	
	//set frz
	if(fr=='N') fr_is_mobile[nt]=FR_MOB_FROZ;
	else if(fr=='B') fr_is_mobile[nt]=FR_MOB_BASE;
	else if(fr=='P') fr_is_mobile[nt]=FR_MOB_PHOS;
	else if(fr=='A') fr_is_mobile[nt]=FR_MOB_FULL;
	else{
	  printf("Invalid fix flag for nucleoside %d .\n", nt);
	  exit(ERR_INPUT);
	}
      }
      at++;inat++;if(inat==N_PARTS_PER_NT){inat=0;nt++;}
    }
    else if(!strcmp(pdbrectyp, "REMARK")){
      //here we just read a comment
#ifdef INPUTVERBOSE
      printf("Comment is: %s\n", pdbrectyp);
#endif
      
    }
    else{
      printf("Record in pdb file not recognized at nucleotide %d.\n", nt);
      exit(ERR_INPUT);
    }
  }
  for(i=0;i<nt_n;i++)
    MC_map_sugar(i, rx, ry, rz);
  //copy glp to temp
  for(i=0;i<nt_n;i++){
    mc_temp_glyc[i]=mc_glyc[i];
    mc_temp_puck[i]=mc_puck[i];
  }
  
  //set bonds
  int nbonds=0, currb;
  nt=0;
  if(chains[nt]==chains[nt+1]) nbonds++;
  mc_nbonds[nt][0]=nbonds;
  if(mc_nbonds[nt][0]>0){
    mc_bondlist[nt]=(int *)malloc(2*mc_nbonds[nt][0]*sizeof(int));
    currb=0;
    if(chains[nt]==chains[nt+1]) {mc_bondlist[nt][currb]=nt+1;mc_bondlist[nt][currb+mc_nbonds[nt][0]];currb++;};
  }
  nbonds=0;
  nt=nt_n-1;
  if(chains[nt]==chains[nt-1]) nbonds++;
  mc_nbonds[nt][0]=nbonds;
  if(mc_nbonds[nt][0]>0){
    mc_bondlist[nt]=(int *)malloc(2*mc_nbonds[nt][0]*sizeof(int));
    currb=0;
    if(chains[nt]==chains[nt-1]) {mc_bondlist[nt][currb]=nt-1;mc_bondlist[nt][currb+mc_nbonds[nt][0]];currb++;};
  }
  for(nt=1;nt<nt_n-1;nt++){
    nbonds=0;
    if(chains[nt]==chains[nt-1]) nbonds++;
    if(chains[nt]==chains[nt+1]) nbonds++;
    mc_nbonds[nt][0]=nbonds;
    if(nbonds>0){
      mc_bondlist[nt]=(int *)malloc(2*mc_nbonds[nt][0]*sizeof(int));
      currb=0;
      if(chains[nt]==chains[nt-1]) {mc_bondlist[nt][currb]=nt-1;mc_bondlist[nt][currb+mc_nbonds[nt][0]];currb++;};
      if(chains[nt]==chains[nt+1]) {mc_bondlist[nt][currb]=nt+1;mc_bondlist[nt][currb+mc_nbonds[nt][0]];currb++;};
    }
  }
  free(chains);
  //we check that no pyrimidines are in SYN state, but they can change their GLYC
  for(i=0;i<nt_n; i++){
    if(is_pyrimidine(i)){
      //mc_types[i*N_PARTS_PER_NT]==TYP_URACIL || mc_types[i*N_PARTS_PER_NT]==TYP_CYTOSINE){
      //if((glp_is_flippable[i]==GLP_GLYC || glp_is_flippable[i]==GLP_BOTH) || (mc_glyc[i]==GLYC_S) ){
      if(mc_glyc[i]==GLYC_S ){
	printf("Some pyrimidine (residue %d) has a SYN conformation. This is not implemented for the moment.\n", i);
	exit(1);
      }
    }
  }
}

void MC_initialize_global(int mc_n, int mcseed, int mpi_id){
  int i;
  int nt_n=mc_n/N_PARTS_PER_NT;
  /* RANDOM SEED */
  idum=(double *)malloc(sizeof(double));
  if(mpi_id>-1)
    *idum=-((double) mcseed+(double)mpi_id);
  else 
    *idum=-((double) mcseed);
  if(*idum>=0){
    printf("ERROR: Wrong sign of random seed!\n");
    exit(ERR_INPUT);
  }

  /* GENERAL INITIALIZATION */
  double lx=100.0, ly=100.0, lz=100.0;
  box_l[0]=lx;
  box_l[1]=ly;
  box_l[2]=lz;
  //gpu_mpc_get_temperature(&temp);
  //mc_target_temp=1;
  /**************************/
  
  /* PARTICLE PROPERTIES */
  mc_pbox=(int **)malloc(sizeof(int *)*mc_n);
  for(i=0;i<mc_n;i++)
    mc_pbox[i]=(int *)malloc(sizeof(int)*DIM);
  mc_types=(int *)malloc(sizeof(int)*mc_n);
  
  mc_temp_pbox=(int **)malloc(sizeof(int *)*mc_n);
  for(i=0;i<mc_n;i++)
    mc_temp_pbox[i]=(int *)malloc(sizeof(int)*DIM);
  /***********************/
  
  /* MC INITIALIZATION */
  MC_initialize_energy_parameters(mpi_id);
  mc_r_cut_sq=mc_r_cut*mc_r_cut;
  
  mc_wc_rcut_sq=mc_wc_rcut*mc_wc_rcut;
  mc_bph_rcut_sq=mc_bph_rcut*mc_bph_rcut;
  mc_nb_rcut_sq=mc_nb_rcut*mc_nb_rcut;
  
  /************************/
  
  /* ALLOCATE LINKED CELLS */
 /*  mc_n_linked_cells=1; */
  /* #ifndef MCDILUTE */
  /*   if(mc_r_cut+vl_skin>0){ */
  /*     for(i=0;i<DIM;i++){ */
  /*       mc_nc[i]=(int)floor(box_l[i]/(mc_r_cut+vl_skin)); */
  /*       mc_n_linked_cells*=mc_nc[i]; */
  /*     } */
  /*   } */
  /*   else */
  /* #endif  */
  /* for(i=0;i<DIM;i++) */
  /*       mc_nc[i]=1; */
  
  /*   mc_linked_cell_l=box_l[0]/mc_nc[0]; */
  
  /*   mc_cells=(int *)malloc((mc_n+mc_n_linked_cells)*sizeof(int)); */
  /*       /\* a somewhat arbitrary initialization: all the particles in the last cell *\/ */
  /*   for(i=0;i<mc_n-1;i++) */
  /*   mc_cells[i]=i+1; */
  /*   mc_cells[mc_n-1]=-1; */
  /*   for(i=mc_n;i<mc_n+mc_n_linked_cells-1;i++) */
  /*     mc_cells[i]=-1; */
  /*   mc_cells[mc_n+mc_n_linked_cells-1]=0; */
  /* */
}

void MC_read_nsolute(int *mc_n, int mpi_id, char *pdbname){
  FILE * pdbfile;
  int nmc=0, l, j;
  int fileflag=1, fileflaggen=1;
  char *pdbrectyp;
  static size_t st_l=0;
  static char *lline=NULL;
    
  if((pdbfile=fopen(pdbname, "r"))==NULL){
    printf("No decent init files found in pdb_inits. Can not run like this.\n");
    exit(ERR_INPUT);
  }
  
  while((l=getline(&lline, &st_l, pdbfile))>6){
    pdbrectyp=strndup(lline, 6);
    if(!strcmp(pdbrectyp, "ATOM  ")){
      nmc++;
    }
  }
  *mc_n=nmc;
}

void MC_initialize_arrays(int mc_n, double **rx, double **ry, double **rz){
  int nmc, i, j;
  int type;
  char tmp_str[MAX_BUFFER];
  double pos[DIM];
  FILE * xyzfile;
  int fileflag=1, fileflaggen=1;
  char xyzname[MAX_BUFFER];
  int nt_n=mc_n/((int)N_PARTS_PER_NT);
  /* allocate memory for rx, ry, rz */
  *rx=(double *)malloc(sizeof(double)*mc_n);
  *ry=(double *)malloc(sizeof(double)*mc_n);
  *rz=(double *)malloc(sizeof(double)*mc_n);
  
  mc_temp_x=(double *)malloc(sizeof(double)*mc_n);
  mc_temp_y=(double *)malloc(sizeof(double)*mc_n);
  mc_temp_z=(double *)malloc(sizeof(double)*mc_n);
  
  //ALLOCATE PUCKERS 
  //MC_init_glycs_and_pucks(nt_n, mpi_id); 
  mc_glyc=(int *)malloc(sizeof(int)*nt_n);
  mc_temp_glyc=(int *)malloc(sizeof(int)*nt_n);
  mc_puck=(int *)malloc(sizeof(int)*nt_n);
  mc_temp_puck=(int *)malloc(sizeof(int)*nt_n);
  glp_is_flippable=(int *)malloc(sizeof(int)*nt_n);
  //defaults
  for(i=0;i<nt_n;i++){
    mc_glyc[i]=GLYC_A;
    mc_puck[i]=PUCK_3;
    glp_is_flippable[i]=GLP_FIXED;
  }
  /* ALLOCATE BONDLISTS */
  mc_nbonds=(int **)malloc(sizeof(int *)*nt_n);
  for(i=0;i<nt_n;i++)
    mc_nbonds[i]=(int *)malloc(sizeof(int)*N_BONDED_INTERACTIONS);
  mc_bondlist=(int **)malloc(sizeof(int *)*nt_n);
  mc_anglelist=(int **)malloc(sizeof(int *)*nt_n);
  mc_anglecenter=(int **)malloc(sizeof(int *)*nt_n);
  //defaults
  for(i=0;i<nt_n;i++)
    for(j=0;j<N_BONDED_INTERACTIONS;j++)
      mc_nbonds[i][j]=0;
  
  
  /* ALLOCATE FRZ ARRAYS */
  fr_is_mobile=(int *)malloc(nt_n*sizeof(int));
  //defaults
  for(i=0;i<nt_n;i++)
    fr_is_mobile[i]=FR_MOB_FULL;
   
  for(i=0;i<N_PARTS_PER_NT;i++)
    for(j=0;j<DIM;j++)
      mc_temp_pbox[i][j]=0;
  
  MC_initialize_sugars(nt_n, rx, ry, rz);
}


void MC_initialize_sugars(int nt_n, double **rx, double **ry, double **rz){
  int i;
  MAP_SUG_GL_X[TYP_ADENINE][GLYC_A]=MAP_SUG_GL_X_A_R;
  MAP_SUG_GL_Y[TYP_ADENINE][GLYC_A]=MAP_SUG_GL_Y_A_R;
  MAP_SUG_GL_Z[TYP_ADENINE][GLYC_A]=MAP_SUG_GL_Z_A_R;
  MAP_SUG_GL_X[TYP_ADENINE][GLYC_H]=MAP_SUG_GL_X_H_R;
  MAP_SUG_GL_Y[TYP_ADENINE][GLYC_H]=MAP_SUG_GL_Y_H_R;
  MAP_SUG_GL_Z[TYP_ADENINE][GLYC_H]=MAP_SUG_GL_Z_H_R;
  MAP_SUG_GL_X[TYP_ADENINE][GLYC_S]=MAP_SUG_GL_X_S_R;
  MAP_SUG_GL_Y[TYP_ADENINE][GLYC_S]=MAP_SUG_GL_Y_S_R;
  MAP_SUG_GL_Z[TYP_ADENINE][GLYC_S]=MAP_SUG_GL_Z_S_R;
  
  MAP_SUG_GL_X[TYP_GUANINE][GLYC_A]=MAP_SUG_GL_X_A_R;
  MAP_SUG_GL_Y[TYP_GUANINE][GLYC_A]=MAP_SUG_GL_Y_A_R;
  MAP_SUG_GL_Z[TYP_GUANINE][GLYC_A]=MAP_SUG_GL_Z_A_R;
  MAP_SUG_GL_X[TYP_GUANINE][GLYC_H]=MAP_SUG_GL_X_H_R;
  MAP_SUG_GL_Y[TYP_GUANINE][GLYC_H]=MAP_SUG_GL_Y_H_R;
  MAP_SUG_GL_Z[TYP_GUANINE][GLYC_H]=MAP_SUG_GL_Z_H_R;
  MAP_SUG_GL_X[TYP_GUANINE][GLYC_S]=MAP_SUG_GL_X_S_R;
  MAP_SUG_GL_Y[TYP_GUANINE][GLYC_S]=MAP_SUG_GL_Y_S_R;
  MAP_SUG_GL_Z[TYP_GUANINE][GLYC_S]=MAP_SUG_GL_Z_S_R;

  MAP_SUG_GL_X[TYP_CYTOSINE][GLYC_A]=MAP_SUG_GL_X_A_Y;
  MAP_SUG_GL_Y[TYP_CYTOSINE][GLYC_A]=MAP_SUG_GL_Y_A_Y;
  MAP_SUG_GL_Z[TYP_CYTOSINE][GLYC_A]=MAP_SUG_GL_Z_A_Y;
  MAP_SUG_GL_X[TYP_CYTOSINE][GLYC_H]=MAP_SUG_GL_X_H_Y;
  MAP_SUG_GL_Y[TYP_CYTOSINE][GLYC_H]=MAP_SUG_GL_Y_H_Y;
  MAP_SUG_GL_Z[TYP_CYTOSINE][GLYC_H]=MAP_SUG_GL_Z_H_Y;
  MAP_SUG_GL_X[TYP_CYTOSINE][GLYC_S]=MAP_SUG_GL_X_S_Y;
  MAP_SUG_GL_Y[TYP_CYTOSINE][GLYC_S]=MAP_SUG_GL_Y_S_Y;
  MAP_SUG_GL_Z[TYP_CYTOSINE][GLYC_S]=MAP_SUG_GL_Z_S_Y;
  
  MAP_SUG_GL_X[TYP_URACIL][GLYC_A]=MAP_SUG_GL_X_A_Y;
  MAP_SUG_GL_Y[TYP_URACIL][GLYC_A]=MAP_SUG_GL_Y_A_Y;
  MAP_SUG_GL_Z[TYP_URACIL][GLYC_A]=MAP_SUG_GL_Z_A_Y;
  MAP_SUG_GL_X[TYP_URACIL][GLYC_H]=MAP_SUG_GL_X_H_Y;
  MAP_SUG_GL_Y[TYP_URACIL][GLYC_H]=MAP_SUG_GL_Y_H_Y;
  MAP_SUG_GL_Z[TYP_URACIL][GLYC_H]=MAP_SUG_GL_Z_H_Y;
  MAP_SUG_GL_X[TYP_URACIL][GLYC_S]=MAP_SUG_GL_X_S_Y;
  MAP_SUG_GL_Y[TYP_URACIL][GLYC_S]=MAP_SUG_GL_Y_S_Y;
  MAP_SUG_GL_Z[TYP_URACIL][GLYC_S]=MAP_SUG_GL_Z_S_Y;
  
  //for(i=0;i<nt_n;i++)
  //MC_map_sugar(i, rx, ry, rz);
}

void MC_initialize_energies(int mc_n, double *rx, double *ry, double *rz){
  //MC_update_linked_lists(mc_n, rx, ry, rz);
  int nt_n=mc_n/N_PARTS_PER_NT;
  MC_initialize_verlet_lists(nt_n, rx, ry, rz);
  MC_swap_local_energies(nt_n, rx,ry,rz);
  
  
  //MC_energy_calc(mc_n, rx, ry, rz, fsx, fsy, fsz);
   /* INITIALIZE WC LISTS */
  MC_init_wc_arrays(nt_n, rx, ry, rz);
  
}

void MC_swap_local_energies(int nt_n, double *rx, double *ry, double *rz){
  //here we dont move particles, but we might change their GLYC state between ANTI and HIGH-ANTI
  int nt_c, flag, glyflag, tempglp;
  double temp;
  #ifdef INPUTVERBOSE
  printf("Calculating local energies for consistency between ANTI and HIGH-ANTI conformations. ONLY ON FLIPPABLE NUCLEOTIDES\n");
  #endif
  for(nt_c=0;nt_c<nt_n;nt_c++){
    if(glp_is_flippable[nt_c]==GLP_BOTH || glp_is_flippable[nt_c]==GLP_GLYC){
      glyflag=0;
      tempglp=glp_is_flippable[nt_c];
      glp_is_flippable[nt_c]=GLP_BOTH;
      MC_copy_nt(nt_c, rx, ry, rz);
      flag=MC_calculate_local_energy(rx, ry, rz, nt_c, &temp, nt_n, -2);
      if(mc_glyc[nt_c]!=mc_temp_glyc[nt_c])
	printf("Fixing GLYC of NT %d from %d to %d\n", nt_c, mc_glyc[nt_c], mc_temp_glyc[nt_c]);
      if(mc_puck[nt_c]!=mc_temp_puck[nt_c])
	printf("Fixing PUCK of NT %d from %d to %d\n", nt_c, mc_puck[nt_c], mc_temp_puck[nt_c]);
      MC_update_glypuck(nt_c);
      glp_is_flippable[nt_c]=tempglp;
    }
  }
  #ifdef INPUTVERBOSE
  printf("GLYCS and PUCKS fixed.\n");
  #endif
}

void MC_print_positions(int mc_n, double ** rx, double **ry, double **rz){
  int i;
  for (i=0;i<mc_n;i++)
    printf("%d %lf %lf %lf\n", i, (*rx)[i], (*ry)[i], (*rz)[i]);
}

/* char *ew_gtl(FILE *source) { */
/*   static size_t st_l2=0; */
/*   int n,i; */
/*   static char *sch_s2=NULL; */
/*   while ((n=getline(&sch_s2,&st_l2,source))>0) { */
/* //    printf("%d %d %s",n,strlen(s),s); */
/* //    assert(s[n-1]=='\n'); */
/*     if (sch_s2[n-1]=='\n') */
/*       sch_s2[n-1]='\0'; */
/*     if (isalpha(sch_s2[0])) */
/*       for(i=0;sch_s2[i];++i) */
/*         if (isspace(sch_s2[i])) */
/*           sch_s2[i--]='\0'; */
/*     if (sch_s2[0]!='#') */
/*       break; */
/*     } */
/*   return n>0?sch_s2:NULL; */
/* } */

void  MC_print_params(int mc_steps, int rand_seed){
  printf("SPQR SIMULATION PARAMETERS\n");  
  printf("Temperature                 : %lf\n", mc_target_temp);
  printf("MC integration steps        : %d\n", mc_steps);
  printf("Saving trajectory every     : %d\n", mc_traj_steps);
  printf("Saving checkpoints every    : %d\n", mc_chkp_steps);
  printf("Random seed                 : %d\n", rand_seed);
  printf("MC nt position displacement : %lf\n", MC_NT_XYZ);
  printf("MC ph position displacement : %lf\n", MC_PH_XYZ);
  printf("MC nt-rotation angle        : %lf\n", MC_NT_ANGLE);
  printf("Cutoff radii                : %lf %lf (wc) %lf (bph) %lf (nb)\n", mc_r_cut, mc_wc_rcut, mc_bph_rcut, mc_nb_rcut);
  printf("Verlet skin                 : %lf\n", vl_skin);
  printf("Energy file                 : %s\n", ENERG_PATH);
  printf("\n");
  
}

void MC_read_params(int *mc_iter, int *rand_a, int mpi_id){
  FILE *mc_params;
  char PARAMS_NAME[MAXSTR];
  sprintf(PARAMS_NAME,"params.spqr");
  char stmp[MAXSTR], s[MAXSTR], s2[MAXSTR], s3[MAXSTR], *line=NULL;  static size_t st_l=0;
  int essential_flags=0, l;
  int wall_flag=0,w;double w1,w2,w3,w4,w5,w6;int w7,itemp;char *tok,*templine;
  //*lx=100.0 ;*ly=100.0;*lz=100.0;
  mc_r_cut=DEFAULT_MC_RCUT;
  mc_wc_rcut=DEFAULT_MC_WC_RCUT;
  mc_bph_rcut=DEFAULT_MC_BPH_RCUT;
  mc_nb_rcut=DEFAULT_MC_NB_RCUT;
  mc_n_types=N_BASES;
  mc_n_bond_types=DEFAULT_MC_N_BOND_TYPES;
  vl_skin=DEFAULT_VL_SKIN;
  mc_chkp_steps=DEFAULT_CHKP_STEPS;
  mc_traj_steps=DEFAULT_TRAJ_STEPS;
  PDB_OUTPUT=0;
  UMBRELLA_TYPE=-1;
  N_UMBRELLAS=0;
  strcpy(ENERG_PATH,".");
  
  MC_NT_XYZ=(double)MC_NT_XYZ_DEF;
  MC_PH_XYZ=(double)MC_PH_XYZ_DEF;
  MC_NT_ANGLE=(double)MC_NT_ANGLE_DEF;
  mc_target_temp=(double)MC_TEMP_DEF;
  *rand_a=(double)MC_RAND_SEED_DEF;
  *mc_iter=(double)MC_ITER_DEF;
  KRG=0;
  RG_target=0;
  mc_params=fopen(PARAMS_NAME, "r");
  
  if(mc_params==NULL){
    printf("No MC parameters file found.\n");
    exit(ERR_INPUT);
  } else {
    //if(mpi_id==0)
    printf("Reading MC parameters from file %s ...\n", PARAMS_NAME);
    //while (s=ew_gtl(mc_params)) {
    while((l=getline(&line, &st_l, mc_params))>0){
      //printf("ret %d  '%s'\n", sscanf(line, "%s", s), s);
      //if(line[0]!='#' && sscanf(line, "%s", stmp)>0 && (line[0]!='S' && (line[1]!='A' && line[1]!='T' ) && (line[2]!='_'))){
      if(line[0]!='#' && sscanf(line, "%s", stmp)>0){
	  if(((line[0]!='S' || (line[1]!='A' && line[1]!='T' ) || (line[2]!='_')) && (line[0]!='P' || line[1]!='T') )){
	    sscanf(line, "%s", s);
	    //printf("read %s - line = %s\n", s2, line);
	    if(!strcmp(s, "TEMPERATURE")) {
	      if(sscanf(line, "%s %lf", s2, &mc_target_temp)!=2){printf("Invalid value of TEMPERATURE in %s\n", PARAMS_NAME);exit(ERR_INPUT);}
	      essential_flags++;
	    }
	    else if(!strcmp(s, "RG_COUPL")) {
	      if(sscanf(line, "%s %lf %lf", s2, &KRG, &RG_target)!=3){printf("Invalid value of RG_COUPL in %s\n", PARAMS_NAME);exit(ERR_INPUT);}
	      essential_flags++;
	    }
	    else if(!strcmp(s, "WALL_COUPL")) {
	      if(sscanf(line, "%s %d ", s2, &N_WALLS)!=2){
		printf("Invalid value of WALL_COUPL in %s\n", PARAMS_NAME);exit(ERR_INPUT);}
	      if(N_WALLS>0){
		wall_epsilon=(double*)malloc(sizeof(double)*N_WALLS);
		wall_sigma  =(double*)malloc(sizeof(double)*N_WALLS);
		wall_A=(double*)malloc(sizeof(double)*N_WALLS);
		wall_B=(double*)malloc(sizeof(double)*N_WALLS);
		wall_C=(double*)malloc(sizeof(double)*N_WALLS);
		wall_D=(double*)malloc(sizeof(double)*N_WALLS);
		WALL_TYPE=(int*)malloc(sizeof(int)*N_WALLS);
		wall_MODSQ=(double *)malloc(sizeof(double)*N_WALLS);
		templine=strndup(line, MAXSTR);
		tok=strtok(templine," ");tok=strtok(NULL," ");
		
		for(w=0;w<N_WALLS;w++){
		  tok=strtok(NULL," ");wall_epsilon[w]=strtof(tok,NULL);
		  tok=strtok(NULL," ");wall_sigma[w]=strtof(tok,NULL);
		  tok=strtok(NULL," ");wall_A[w]=strtof(tok,NULL);
		  tok=strtok(NULL," ");wall_B[w]=strtof(tok,NULL);
		  tok=strtok(NULL," ");wall_C[w]=strtof(tok,NULL);
		  tok=strtok(NULL," ");wall_D[w]=strtof(tok,NULL);
		  tok=strtok(NULL," ");WALL_TYPE[w]=atoi(tok);
		  wall_MODSQ[w]=sqrt(SQ(wall_A[w])+SQ(wall_B[w])+SQ(wall_C[w]));
		  wall_epsilon[w]*=4.0;
		}
		free(templine);
	      }
	      wall_flag++;
	    }
	    else if(!strcmp(s, "UMBRELLA_SAMPLING")) {
	      if(sscanf(line, "%s %d ", s2, &N_UMBRELLAS)!=2){
		printf("Invalid value of WALL_COUPL in %s\n", PARAMS_NAME);exit(ERR_INPUT);}
	      if(N_UMBRELLAS>0){
		templine=strndup(line, MAXSTR);
		tok=strtok(templine," ");tok=strtok(NULL," ");
		tok=strtok(NULL," ");UMBRELLA_TYPE=atoi(tok);
		tok=strtok(NULL," ");UCM[0]=strtof(tok,NULL);
		tok=strtok(NULL," ");UCM[1]=strtof(tok,NULL);
		tok=strtok(NULL," ");UCM[2]=strtof(tok,NULL);
		tok=strtok(NULL," ");UCMK0=strtof(tok,NULL);
		tok=strtok(NULL," ");UCMK1=strtof(tok,NULL);
		tok=strtok(NULL," ");UCMK2=strtof(tok,NULL);
		printf("UMBRELLA SAMPLING constraint : at ( %lf %lf %lf )  with constants ( %lf %lf %lf ) \n",UCM[0], UCM[1], UCM[2],UCMK0,UCMK1,UCMK2);
		free(templine);
	      }
	    }
	    else if(!strcmp(s, "PDB_OUTPUT")) {
	      if(sscanf(line, "%s %d", s2, &PDB_OUTPUT)!=2){printf("Invalid value of PDB_OUTPUT in %s\n", PARAMS_NAME);exit(ERR_INPUT);}
	      essential_flags++;
	    }
	    else if(!strcmp(s, "ENERGS_PATH")) {
	      if(sscanf(line, "%s %s", s2, s3)!=2){printf("Invalid value of ENERGS_PATH in %s\n", PARAMS_NAME);exit(ERR_INPUT);}
	      strcpy(ENERG_PATH, s3);
	      essential_flags++;
	    }
	    else if(!strcmp(s, "MC_NT_XYZ")) {
	      if(sscanf(line, "%s %lf", s2, &MC_NT_XYZ)!=2){printf("Invalid value of MC_NT_XYZ in %s\n", PARAMS_NAME);exit(ERR_INPUT);}
	      essential_flags++;
	    }
	    else if(!strcmp(s, "MC_PH_XYZ")) {
	      if(sscanf(line, "%s %lf", s2, &MC_PH_XYZ)!=2){printf("Invalid value of MC_PH_XYZ in %s\n", PARAMS_NAME);exit(ERR_INPUT);}
	      essential_flags++;
	    }
	    else if(!strcmp(s, "MC_STEPS")) {
	      if(sscanf(line, "%s %d", s2, mc_iter)!=2){printf("Invalid value of MC_STEPS in %s\n", PARAMS_NAME);exit(ERR_INPUT);}
	      essential_flags++;
	    }
	    else if(!strcmp(s, "RANDOM_SEED")) {
	      if(sscanf(line, "%s %d", s2, rand_a)!=2){printf("Invalid value of RANDOM_SEED in %s\n", PARAMS_NAME);exit(ERR_INPUT);}
	      essential_flags++;
	    }
	    else if(!strcmp(s, "MC_NT_ANGLE")) {
	      if(sscanf(line, "%s %lf %lf", s2, &MC_NT_ANGLE, &MC_BB_ANGLE)!=3){printf("Invalid values of MC_NT_ANGLE in %s\n", PARAMS_NAME);exit(ERR_INPUT);}
	      MC_NT_ANGLE_COS=cos(MC_NT_ANGLE);
	      MC_NT_ANGLE_SIN=sin(MC_NT_ANGLE);
	      MC_BB_ANGLE_COS=cos(MC_BB_ANGLE);
	      MC_BB_ANGLE_SIN=sin(MC_BB_ANGLE);
	      essential_flags++;
	    }
	    /* else if(!strcmp(s, "MC_RCUT")) { */
	    /*   if(sscanf(line, "%s %lf %lf %lf %lf", s2, &mc_r_cut, &mc_wc_rcut, &mc_bph_rcut, &mc_nb_rcut)!=5){printf("Invalid values of MC_RCUT in %s\n", PARAMS_NAME);exit(ERR_INPUT);} */
	    /*   essential_flags++; */
	    /* } */
	    /* else if(!strcmp(s, "VL_SKIN")) { */
	    /*   if(sscanf(line, "%s %lf", s2, &vl_skin)!=2){printf("Invalid value of VL_SKIN in %s\n", PARAMS_NAME);exit(ERR_INPUT);} */
	    /*   essential_flags++; */
	    /* } */
	    else if(!strcmp(s, "MC_TRAJ_STEPS")) {
	      if(sscanf(line, "%s %d", s2, &mc_traj_steps)!=2){printf("Invalid value of MC_TRAJ_STEPS in %s\n", PARAMS_NAME);exit(ERR_INPUT);}
	      essential_flags++;
	    }
	    else if(!strcmp(s, "MC_CHKP_STEPS")) {
	      if(sscanf(line, "%s %d", s2, &mc_chkp_steps)!=2){printf("Invalid value of MC_CHKP_STEPS in %s\n", PARAMS_NAME);exit(ERR_INPUT);}
	      essential_flags++;
	    }
	    else {printf("Parameter %s not recognized in %s\n", s, PARAMS_NAME);exit(ERR_INPUT);}
	  }
	}
    }
    if(essential_flags!=11){printf("Missing parameters in %s!\n", PARAMS_NAME); exit(ERR_INPUT);}
    if(wall_flag==0){N_WALLS=0;
      wall_epsilon=(double *)malloc(sizeof(double));wall_epsilon[0]=0;
      wall_sigma  =(double *)malloc(sizeof(double));wall_sigma[0]=0;
      wall_A=(double *)malloc(sizeof(double));wall_A[0]=0;
      wall_B=(double *)malloc(sizeof(double));wall_B[0]=0;
      wall_C=(double *)malloc(sizeof(double));wall_C[0]=0;
      wall_D=(double *)malloc(sizeof(double));wall_D[0]=0;
      WALL_TYPE=(int *)malloc(sizeof(int));WALL_TYPE[0]=0;
      wall_MODSQ=(double *)malloc(sizeof(double));wall_MODSQ[0]=0;
      //wall_epsilon=0; wall_sigma=0, wall_A=0; wall_B=0; wall_C=0;wall_D=0;WALL_TYPE=0; // defaults}
    }
    if(UMBRELLA_TYPE==-1){
      UCM[0]=0;
      UCM[1]=0;
      UCM[2]=0;
      UCMK0=0.0;
      UCMK1=0.0;
      UCMK2=0.0;
    }
  }
  fclose(mc_params);
  /***** check some parameters *****/
  
  /*********************************/
  //if(mpi_id==0)
  if(mpi_id>-1)
  MC_print_params(*mc_iter, *rand_a);
}
