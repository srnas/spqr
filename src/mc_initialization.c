#include "mc_initialization.h"

void MC_initialize(int mc_n, double **rsx, double **rsy, double **rsz, double lx, double ly, double lz, int mc_read_conf_flag, int mcseed, int ini, int mpi_id){
  MC_initialize_global(mc_n, lx, ly, lz, mcseed, mpi_id);
  //if(mc_read_conf_flag==READ_MC_CONF){
  MC_initialize_positions(mc_n,rsx,rsy,rsz,mc_read_conf_flag, mpi_id);
  /* INITIALIZE CONFIGURATION FILE */
  MC_initialize_save_configs(mc_n, ini, mpi_id);
  MC_initialize_energies(mc_n, *rsx, *rsy, *rsz);
}

void MC_initialize_global(int mc_n, double lx, double ly, double lz, int mcseed, int mpi_id){
  int i;
  int nt_n=mc_n/N_PARTS_PER_NT;
  /* RANDOM SEED */
  idum=(double *)malloc(sizeof(double));
  *idum=-((double) mcseed+(double)mpi_id);
  
  /* GENERAL INITIALIZATION */
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
  mc_n_linked_cells=1;
#ifndef MCDILUTE
  if(mc_r_cut+vl_skin>0){
    for(i=0;i<DIM;i++){
      mc_nc[i]=(int)floor(box_l[i]/(mc_r_cut+vl_skin));
      mc_n_linked_cells*=mc_nc[i];
    }
  }
  else
#endif 
   for(i=0;i<DIM;i++)
      mc_nc[i]=1;
  
  mc_linked_cell_l=box_l[0]/mc_nc[0];
  
  mc_cells=(int *)malloc((mc_n+mc_n_linked_cells)*sizeof(int));
  /* a somewhat arbitrary initialization: all the particles in the last cell */
  for(i=0;i<mc_n-1;i++)
    mc_cells[i]=i+1;
  mc_cells[mc_n-1]=-1;
  for(i=mc_n;i<mc_n+mc_n_linked_cells-1;i++)
    mc_cells[i]=-1;
  mc_cells[mc_n+mc_n_linked_cells-1]=0;

  /*************************/
  
  /* INITIALIZE BONDLISTS */
  mc_nbonds=(int **)malloc(sizeof(int *)*nt_n);
  for(i=0;i<nt_n;i++)
    mc_nbonds[i]=(int *)malloc(sizeof(int)*N_BONDED_INTERACTIONS);
  mc_bondlist=(int **)malloc(sizeof(int *)*nt_n);
  mc_anglelist=(int **)malloc(sizeof(int *)*nt_n);
  mc_anglecenter=(int **)malloc(sizeof(int *)*nt_n);

#ifdef STWCBONDED
  mc_stbondlist=(int **)malloc(sizeof(int *)*nt_n);
  mc_wcbondlist=(int **)malloc(sizeof(int *)*nt_n);
#endif
#ifdef DIHEDRALS
  mc_dihedrallist=(int **)malloc(sizeof(int *)*nt_n);
#endif
  if(MC_open_bondlist_file()){
    //printf("Open bonds.dat with %d types of bonded interactions\n", (int) N_BONDED_INTERACTIONS);
    for(i=0;i<nt_n;i++){
      MC_read_bondlist_from_file(i, mpi_id);
    }
    MC_close_bondlist_file();
  }else{
    printf("Using default bonds.\n");
    MC_default_bonds(nt_n);
  }

  /* int j;  */
  /* for(i=0;i<nt_n;i++){ */
  /*   printf("%d :  ", i); */
  /*   for(j=0;j<mc_nbonds[i][2];j++) */
  /*     printf("st %d ", mc_stbondlist[i][j]); */
  /*   for(j=0;j<mc_nbonds[i][3];j++) */
  /*     printf("wc %d ", mc_wcbondlist[i][j]); */
    
    
  /*   printf("\n"); */
  /* } */
 
  /************************/
  
  
#ifdef FROZEN
  MC_init_frozen_file(nt_n);
#endif
  if(mpi_id==0)
    MC_print_parameters(mc_n);
  //INITIALIZE PUCKERS
  MC_init_glycs_and_pucks(nt_n);

#ifdef SECONDSTC
  MC_init_sec_str_c_lists(nt_n);
#endif

}

#ifdef SECONDSTC
void MC_init_sec_str_c_lists(int nt_n){
  int i, p;
  int temp1, temp2, temp3, temp4;
  char filename[256];
  sprintf(filename, "sec_strc.lst");
  FILE *sscfile;
  
  wc_sstruct_N=(int *)malloc(sizeof(int)*nt_n);
  wc_sstruct_neigh=(int **)malloc(sizeof(int *)*nt_n);
  for(i=0;i<nt_n;i++){
    wc_sstruct_N[i]=0;
  }
  if((sscfile=fopen(filename,"r"))==NULL){
    printf("List of secondary_structure constrain not found - ignoring it\n");
  }
  else{
    printf("SECONDARY STRUCTURE IS CONSTRAINED!!\n");
    for(i=0;i<nt_n;i++){
      fscanf(sscfile, "%d%d", &temp1, &temp2);
      if(temp1!=i){
	printf("Invalid index in sec_strc.lst file at line %d\n", i);
	exit(ERR_INPUT);
      }
      wc_sstruct_N[i]=temp2; //THIS IS THE NUMBER OF NT CONSTRAINED TO INTERACT WITH THE CURRENT NT. 0 MEANS NO RESTRICTION, -1 MEANS FULL RESTRICTION
      wc_sstruct_neigh[i]=(int *)malloc(sizeof(int)*2*wc_sstruct_N[i]);
      for(p=0;p<wc_sstruct_N[i];p++){
	fscanf(sscfile, "%d%d", &temp3, &temp4);
	wc_sstruct_neigh[i][2*p]=temp3; //THIS IS THE INDEX OF THE NEIGHBOR
	wc_sstruct_neigh[i][2*p+1]  =temp4; //THIS MUST BE THE TYPE OF CONSTRAINT - 0 , ONLY BASE-PAIR, 1 - ALL, 2 - THIS BASE + PHOS , 3 - THIS PHOS + BASE
	
      }
      
    }
    /* for(p=0;p<nt_n;p++){ */
    /*   printf("%d %d  ", p, wc_sstruct_N[p]); */
    /*   for(i=0;i<wc_sstruct_N[p];i++){ */
    /* 	printf("%d %d  ", wc_sstruct_neigh[p][1], wc_sstruct_neigh[p][0]); */
    /*   } */
    /*   printf("\n"); */
    /* } */
    fclose(sscfile);
  }
}
#endif
#ifdef FROZEN
void MC_init_frozen_file(int nt_n){
  int i, temp1, temp2;
  fr_is_mobile=(int *)malloc(nt_n*sizeof(int));
  for(i=0;i<nt_n;i++)
    fr_is_mobile[i]=FR_MOB_FULL;
  char filename[256];
  sprintf(filename, "mobile.lst");
  FILE *frfile;
  if((frfile=fopen(filename,"r"))==NULL){
    printf("List of mobile particles not found. All atoms are mobile by default.\n");
  }
  else{
    printf("List of mobile particles found.\n");
    for(i=0;i<nt_n;i++){
      fscanf(frfile, "%d %d", &temp1, &temp2);
      if(temp1!=i){
	printf("Invalid index in mobile.lst file at line %d\n", i);
	exit(ERR_INPUT);
      }
      if(temp2<FR_MOB_FROZ || temp2 > FR_MOB_FULL){
	printf("Invalid value in mobile.lst file at line %d\n", i);
	exit(ERR_INPUT);
      }
      fr_is_mobile[temp1]=temp2;
    }
    fclose(frfile);
  }
}
#endif

void MC_init_glycs_and_pucks(int nt_n){ 
  int i, temp1, temp3;
  char temp2;
  char filename[256];
  mc_glyc=(int *)malloc(sizeof(int)*nt_n);
  mc_temp_glyc=(int *)malloc(sizeof(int)*nt_n);
  
  mc_puck=(int *)malloc(sizeof(int)*nt_n);
  mc_temp_puck=(int *)malloc(sizeof(int)*nt_n);
  
  
  sprintf(filename, "gly_pck.dat");
  FILE *glycfile;
  if((glycfile=fopen(filename, "r"))==NULL){
    printf("No glycosidic indexes file found. Setting all glycosidic states to ANTI and all puckers to C3' endo\n");
    for(i=0;i<nt_n;i++){
      mc_glyc[i]=GLYC_A;
      mc_puck[i]=PUCK_3;
    }
  }
  else{
    for(i=0;i<nt_n;i++){
      fscanf(glycfile, "%d %c %d", &temp1, &temp2, &temp3);
      if(temp1!=i){
	printf("Wrong index for %d in gly_pck.dat file.\n", i);
	exit(ERR_INPUT);
      }
      if(temp2=='A') mc_glyc[i]=GLYC_A;
      else if(temp2=='H') mc_glyc[i]=GLYC_H;
      else if(temp2=='S') mc_glyc[i]=GLYC_S;
      else {
	printf("Invalid GLYCOSIDIC index for %d in gly_pck.dat file.\n", i);
	exit(ERR_INPUT);
      }
      
      if(temp3==3) mc_puck[i]=PUCK_3;
      else if(temp3==2) mc_puck[i]=PUCK_2;
      else {
	printf("Invalid PUCKER index for %d in gly_pck.dat file.\n", i);
	exit(ERR_INPUT);
      }
      
      
    }
    fclose(glycfile);
  }
  for(i=0;i<nt_n;i++){
    mc_temp_glyc[i]=mc_glyc[i];
    mc_temp_puck[i]=mc_puck[i];
  }
}

void MC_print_parameters(int mc_n){
  int i,j;
  printf("\nMC parameters:\n");
  printf("Number of particles                                     : %d\n", mc_n);
  printf("MC linked cell length                                   : %lf\n",mc_linked_cell_l);
  printf("Number of MC linked cells                               : %d\n", mc_n_linked_cells);
  printf("Maximum cutoff radius between nucleotides               : %lf\n", mc_r_cut);
  printf("Number of particle types (for non-bonded interactions)  : %d\n", mc_n_types);
  printf("Number of bond types                                    : %d  %d", mc_n_bond_types, mc_n_ang_types);
#ifdef DIHEDRALS
  printf("  %d", mc_n_dih_types);
#endif
  printf("\n");
  printf("Skin for Verlet list                                    : %lf\n", vl_skin);
  /* printf("LJ epsilon (for non-bonded interactions)                : "); */
  /* for(i=0;i<mc_n_types;i++) */
  /*   printf(" ( %d %d %lf ) ", i, i, mc_lj_eps[i][i]); */
  /* for(i=0;i<mc_n_types;i++) */
  /*   for(j=i+1;j<mc_n_types;j++) */
  /*     printf(" ( %d %d %lf ) ", i,j,mc_lj_eps[i][j]); */
  /* printf("\n"); */
  
  /* printf("LJ sigma (for non-bonded interactions)                  : "); */
  /* for(i=0;i<mc_n_types;i++) */
  /*   printf("%lf ", mc_lj_sig[i][i]); */
  /* for(i=0;i<mc_n_types;i++) */
  /*   for(j=i+1;j<mc_n_types;j++) */
  /*     printf("%lf ", mc_lj_sig[i][j]); */
  /* printf("\n"); */
  
  /* printf("Harmonic bond K                                         : "); */
  /* for(i=0;i<mc_n_bond_types;i++) */
  /*   printf("%lf ", mc_harm_k[i]); */
  /* printf("\n"); */
  
  /* printf("Harmonic bond R                                         : "); */
  /* for(i=0;i<mc_n_bond_types;i++) */
  /*   printf("%lf ", mc_harm_r[i]); */
  /* printf("\n"); */
  
  /* printf("Angle constant K                                        : "); */
  /* for(i=0;i<mc_n_ang_types;i++) */
  /*   printf("%lf ", mc_ang_k[i]); */
  /* printf("\n"); */
  
  /* printf("Angle constant THETA                                    : "); */
  /* for(i=0;i<mc_n_ang_types;i++) */
  /*   printf("%lf ", mc_ang_th[i]); */
  /* printf("\n"); */

#ifdef FROZEN
  printf("Some particles are MOBILE: ");
  for(i=0;i<mc_n/N_PARTS_PER_NT;i++){
    if(fr_is_mobile[i]==FR_MOB_FULL) printf(" %d FULL ", i);
    if(fr_is_mobile[i]==FR_MOB_BASE) printf(" %d (N) ", i);
    if(fr_is_mobile[i]==FR_MOB_PHOS) printf(" %d (P) ", i);
    if(fr_is_mobile[i]==FR_MOB_FROZ) printf(" %d FROZEN ", i);
  }
  printf("\n");
#endif
#ifdef DIHEDRALS
  printf("Dihedral constant K                                     : ");
  for(i=0;i<mc_n_dih_types;i++)
    printf("%lf ", mc_dih_k[i]);
  printf("\n");
  
  printf("Dihedral constant THETA                                 : ");
  for(i=0;i<mc_n_dih_types;i++)
    printf("%lf ", mc_dih_phi[i]);
  printf("\n");
  
  printf("Dihedral constant N                                     : ");
  for(i=0;i<mc_n_dih_types;i++)
    printf("%lf ", mc_dih_n[i]);
  printf("\n");
#endif
  
  printf("MC random seed                                          : %lf\n",*idum); 
  printf("MC configurations saved every %d MPC steps.\n", mc_chk_freq);
}

void MC_read_nsolute(int *mc_n, int mpi_id){
  FILE * xyzfile;
  FILE * xyzfile_gen;
  
  char xyzname[256];
  int nmc;
  int fileflag=1, fileflaggen=1;
  sprintf(xyzname,"xyz_inits/mc.p%02d.xyz", mpi_id);
  if((xyzfile=fopen(xyzname, "r"))==NULL){
    fileflag=0;
    sprintf(xyzname,"xyz_inits/mc.xyz");
    if((xyzfile=fopen(xyzname, "r"))==NULL){
      printf("No file  mc.p%02d.xyz nor mc.xyz found. Number of MC particles set to default %d.\n", mpi_id, *mc_n);
      fileflaggen=0;
    }
  }
  //printf("fileflags %d and %d\t%s\n", fileflag, fileflaggen, xyzname);
  if(fileflag==1 || fileflaggen==1){
    if(!(fscanf(xyzfile, "%d", &nmc))){
      fprintf(stderr, "Invalid number of MC particles in xyz file.\n");
      exit(ERR_INPUT);
    }
    else{
      *mc_n=nmc;
      printf("Number of MC particles %d read from file %s\n", nmc,xyzname);
      fclose(xyzfile);
    }
  }
}

void MC_initialize_positions(int mc_n, double **rx, double **ry, double **rz, int mc_read_conf_flag, int mpi_id){
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
  
  switch(mc_read_conf_flag){
  case READ_MC_CONF:
    //printf("Positions to be read from file mc.%02d.xyz.\n", mpi_id);
    sprintf(xyzname,"xyz_inits/mc.p%02d.xyz", mpi_id);
    if((xyzfile=fopen(xyzname, "r"))==NULL){
      fileflag=0;
      sprintf(xyzname,"xyz_inits/mc.xyz");
      if((xyzfile=fopen(xyzname, "r"))==NULL){
	//printf("No file  mc.%02d.xyz nor mc.xyz found. Number of MC particles set to default %d.\n", mpc_id, xyzname, *mc_n);
	fileflaggen=0;
	
	//printf("No file %s. Therefore, the MC positions will be left as they are!\n", xyzname, (int) DEFAULT_MC_MASS);
	printf("No file %s found.\n", xyzname);
	printf("MC particles will be initialized at random positions.\n");
	for(i=0;i<mc_n;i++){
	  for(j=0;j<DIM;j++){
	    pos[j]=rand_d(box_l[j]);
	    mc_pbox[i][j]=0;
#ifdef PBC
	    mc_pbox[i][j]=(int)floor(pos[j]/box_l[j]);
	    fold_coordinate(&pos[j],box_l[j],0);
#endif
	  }
	  (*rx)[i]=pos[0];
	  (*ry)[i]=pos[1];
	  (*rz)[i]=pos[2];
	  mc_types[i]=0;
	}
      }
    }
    if(fileflag==1 || fileflaggen==1){
      printf("MC positions read from file %s.\n", xyzname);
      if(!(fscanf(xyzfile, "%d", &nmc))){
	fprintf(stderr, "Invalid number of MC particles in xyz file.\n");
	exit(ERR_INPUT);
      }
      if(mc_n!=nmc){
	fprintf(stderr, "Number of MC particles of file %s does not match the input file.\n", xyzname);
	exit(ERR_INPUT);
      }
      /* the comment line */
      fgets(tmp_str, MAX_BUFFER, xyzfile);
      fgets(tmp_str, MAX_BUFFER, xyzfile);
      //printf("%s file says: ", xyzname);
      //puts(tmp_str);
      
      for(i=0;i<mc_n;i++){
	if(!(fscanf(xyzfile, "%d%lf%lf%lf", &type, &pos[0], &pos[1], &pos[2]))){
	  fprintf(stderr, "Wrong syntax in  %s file at line %d.\n", xyzname, i+3);
	  exit(ERR_INPUT);
	} else {
	  fgets(tmp_str, MAX_BUFFER, xyzfile);
	  /* fold coordinates */
	  for(j=0;j<DIM;j++){
	    mc_pbox[i][j]=0;
#ifdef PBC
	    mc_pbox[i][j]=(int)floor(pos[j]/box_l[j]);
	    fold_coordinate(&pos[j], box_l[j],0);
#endif
	  }
	  (*rx)[i]=pos[0];
	  (*ry)[i]=pos[1];
	  (*rz)[i]=pos[2];
	  if(type>=mc_n_types){
	    fprintf(stderr, "Type of particle %d exceeds the number of types %d.\n", type, mc_n_types);
	    exit(ERR_INPUT);
	  }
	  mc_types[i]=type;
	}
	//printf("%lf %lf %lf\n", (*rx)[i], (*ry)[i], (*rz)[i]);
      }
      fclose(xyzfile);
    } 
    break;
    
  case READ_MC_CHECK:
    printf("Positions to be read from checkpoint, and lastmc.xyz file for unfolded positions and particle types.\n");
    printf("NOT IMPLEMENTED YET.\n");
    exit(ERR_INIT);
    char chkname[11]="lastmc.xyz";
    if((xyzfile=fopen(chkname, "r"))==NULL){
      //printf("No file %s. Therefore, the MC positions will be left as they are!\n", xyzname, (int) DEFAULT_MC_MASS);
      printf("No file %s found.\n", chkname);
      exit(ERR_INPUT);
    }
    else {
      printf("File %s found.\n", chkname);
      if(!(fscanf(xyzfile, "%d", &nmc))){
	fprintf(stderr, "Invalid number of MC particles in xyz file.\n");
	exit(ERR_INPUT);
      }
      if(mc_n!=nmc){
	fprintf(stderr, "Number of MC particles of file %s does not match the input file.\n", chkname);
	exit(ERR_INPUT);
      }
      /* the comment line */
      fgets(tmp_str, MAX_BUFFER, xyzfile);
      fgets(tmp_str, MAX_BUFFER, xyzfile);
      printf("%s file says: ", chkname);
      puts(tmp_str);
      for(i=0;i<mc_n;i++){
	if(!(fscanf(xyzfile, "%d%lf%lf%lf", &type, &pos[0], &pos[1], &pos[2]))){
	  fprintf(stderr, "Wrong syntax in  %s file at line %d.\n", chkname, i+3);
	  exit(ERR_INPUT);
	} else {
	  fgets(tmp_str, MAX_BUFFER, xyzfile);
	  /* fold coordinates */
	  for(j=0;j<DIM;j++){
	    mc_pbox[i][j]=0;
#ifdef PBC
	    mc_pbox[i][j]=(int)floor(pos[j]/box_l[j]);
	    //fold_coordinate(&pos[j],box_l[j],0);
#endif
	  }
	  if(type>=mc_n_types){
	    fprintf(stderr, "Type of particle %d exceeds the number of types %d.\n", type, mc_n_types);
	    exit(ERR_INPUT);
	  }
	  mc_types[i]=type;
	}
      }
      fclose(xyzfile);
    }
    break;
    
  case DONT_READ_MC_CONF:
    printf("MC particles will be initialized at random positions.\n");
    for(i=0;i<mc_n;i++){
      for(j=0;j<DIM;j++){
	pos[j]=rand_d(box_l[j]);
	mc_pbox[i][j]=0;
#ifdef PBC
	mc_pbox[i][j]=(int)floor(pos[j]/box_l[j]);
	fold_coordinate(&pos[j],box_l[j],0);
#endif
      }
      (*rx)[i]=pos[0];
      (*ry)[i]=pos[1];
      (*rz)[i]=pos[2];
      mc_types[i]=0;
    }
    
    break;
  }
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
    
  
  
  
  for(i=0;i<nt_n;i++)
    MC_map_sugar(i, rx, ry, rz);
}

void MC_initialize_energies(int mc_n, double *rx, double *ry, double *rz){
  //MC_update_linked_lists(mc_n, rx, ry, rz);
  int nt_n=mc_n/N_PARTS_PER_NT;
  MC_initialize_verlet_lists(nt_n, rx, ry, rz);
  //MC_energy_calc(mc_n, rx, ry, rz, fsx, fsy, fsz);
   /* INITIALIZE WC LISTS */
  MC_init_wc_arrays(nt_n, rx, ry, rz);

  
}

void MC_print_positions(int mc_n, double ** rx, double **ry, double **rz){
  int i;
  for (i=0;i<mc_n;i++)
    printf("%d %lf %lf %lf\n", i, (*rx)[i], (*ry)[i], (*rz)[i]);
}

char *ew_gtl2(FILE *source) {
  static size_t st_l2=0;
  int n,i;
  static char *sch_s2=NULL;
  while ((n=getline(&sch_s2,&st_l2,source))>0) {
//    printf("%d %d %s",n,strlen(s),s);
//    assert(s[n-1]=='\n');
    if (sch_s2[n-1]=='\n')
      sch_s2[n-1]='\0';
    if (isalpha(sch_s2[0]))
      for(i=0;sch_s2[i];++i)
        if (isspace(sch_s2[i]))
          sch_s2[i--]='\0';
    if (sch_s2[0]!='#')
      break;
    }
  return n>0?sch_s2:NULL;
}

void  MC_print_params(double lx, double ly, double lz, int mc_steps, int rand_seed){
  printf("MC simulation parameters:\n");  
  printf("System size              : %lf %lf %lf\n", lx, ly, lz);
  printf("Temperature              : %lf\n", mc_target_temp);
  printf("MC total steps           : %d\n", mc_steps);
  printf("Random seed              : %d\n", rand_seed);
  printf("MC nt position displacement : %lf\n", MC_NT_XYZ);
  printf("MC ph position displacement : %lf\n", MC_PH_XYZ);
  printf("MC nt-rotation angle     : %lf\n", MC_NT_ANGLE);
  printf("\n");
}

void MC_read_params(double *lx, double *ly, double *lz, int *mc_iter, int *rand_a, int mpi_id){
  FILE *mc_params;
  char *s;
  
  MC_NT_XYZ=(double)MC_NT_XYZ_DEF;
  MC_PH_XYZ=(double)MC_PH_XYZ_DEF;
  MC_NT_ANGLE=(double)MC_NT_ANGLE_DEF;
  mc_target_temp=(double)MC_TEMP_DEF;
  *rand_a=(double)MC_RAND_SEED_DEF;
  *mc_iter=(double)MC_ITER_DEF;
  mc_params=fopen("input.mc", "r");
  if(mc_params==NULL){
    printf("No MC parameters file found. Using default values.\n");
  } else {
    if(mpi_id==0)
      printf("Reading MC parameters from file.\n");
    while (s=ew_gtl2(mc_params)) {
      if(!strcmp(s, "DIMENSIONS")) {
	assert(s=ew_gtl2(mc_params));
	assert(sscanf(s, "%lf%lf%lf", lx, ly, lz)==3);
      }
      if(!strcmp(s, "TEMPERATURE")) {
	assert(s=ew_gtl2(mc_params));
	assert(sscanf(s, "%lf", &mc_target_temp)==1);
      }
      if(!strcmp(s, "MC_NT_XYZ")) {
	assert(s=ew_gtl2(mc_params));
	assert(sscanf(s, "%lf", &MC_NT_XYZ)==1);
      }
      if(!strcmp(s, "MC_PH_XYZ")) {
	assert(s=ew_gtl2(mc_params));
	assert(sscanf(s, "%lf", &MC_PH_XYZ)==1);
      }
      
      if(!strcmp(s, "MC_STEPS")) {
	assert(s=ew_gtl2(mc_params));
	assert(sscanf(s, "%d", mc_iter)==1);
      }
      if(!strcmp(s, "RANDOM_SEED")) {
	assert(s=ew_gtl2(mc_params));
	assert(sscanf(s, "%d", rand_a)==1);
      }
      if(!strcmp(s, "MC_NT_ANGLE")) {
	assert(s=ew_gtl2(mc_params));
	assert(sscanf(s, "%lf%lf", &MC_NT_ANGLE, &MC_BB_ANGLE)==2);
	MC_NT_ANGLE_COS=cos(MC_NT_ANGLE);
	MC_NT_ANGLE_SIN=sin(MC_NT_ANGLE);
	MC_BB_ANGLE_COS=cos(MC_BB_ANGLE);
	MC_BB_ANGLE_SIN=sin(MC_BB_ANGLE);
      }
    }
  }
  fclose(mc_params);
  /***** check some parameters *****/

  /*********************************/
  if(mpi_id==0)
    MC_print_params(*lx, *ly, *lz, *mc_iter, *rand_a);
}
