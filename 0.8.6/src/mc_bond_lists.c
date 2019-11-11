#include "mc_bond_lists.h"

int MC_open_bondlist_file(){
  char filename[10]="bonds.dat";
  if((mc_bond_file = fopen(filename, "r"))==NULL) {
    printf("No bond file, no bonds loaded.\nCreate file %s if your intention is to include bonded interactions.\n", filename);
    return 0;
    //exit(ERR_INPUT);
    
  }
  return 1;
  
}

void MC_read_bondlist_from_file(int nt, int mpi_id){
  int i;
  int cbond;
  int temp;
  //if(mpi_id==0)
  //printf("Open bonds.dat with %d types of bonded interactions\n", (int) N_BONDED_INTERACTIONS);
  //first the particle index, for human readability
  if(!(fscanf(mc_bond_file, "%d", &temp))){
    fprintf(stderr, "Invalid nucleotide id at line %d of bond file.\n", nt);
    exit(ERR_INPUT);
  }
  if(temp!=nt){
    fprintf(stderr, "Invalid nucleotide id at line %d of bond file.\n", nt);
    exit(ERR_INPUT);
  }
  
  for(i=0;i<N_BONDED_INTERACTIONS;i++)
    if(!(fscanf(mc_bond_file, "%d", &mc_nbonds[nt][i]))){
      fprintf(stderr, "Invalid number of bonds for particle %d.\n", nt);
      exit(ERR_INPUT);
      
    }
  //printf("Part %d has %d and %d bonds.\n", nt,mc_nbonds[nt][0],mc_nbonds[nt][1]);//,mc_nbonds[nt][2],mc_nbonds[nt][3]);
  
  //mc_nbonds[p][i]=nbonds[i];
  fflush(stdout);
  
  /* backbone bonds */
  if(mc_nbonds[nt][0]>0){
    mc_bondlist[nt]=(int *)malloc(2*mc_nbonds[nt][0]*sizeof(int));
    for(i=0;i<mc_nbonds[nt][0];i++){
      if(!(fscanf(mc_bond_file, "%d", &cbond))){
	fprintf(stderr, "Invalid bond index for particle %d.\n", nt);
	exit(ERR_INPUT);
      }
      mc_bondlist[nt][i]=cbond;
      
      if(!(fscanf(mc_bond_file, "%d", &cbond))){
	fprintf(stderr, "Invalid bond type for particle %d.\n", nt);
	exit(ERR_INPUT);
      }
      mc_bondlist[nt][i+mc_nbonds[nt][0]]=cbond;
    }
  }
  
  /* angles */
  if(mc_nbonds[nt][1]>0){
    mc_anglelist[nt]=(int *)malloc(2*mc_nbonds[nt][1]*sizeof(int));
    mc_anglecenter[nt]=(int *)malloc(mc_nbonds[nt][1]*sizeof(int));
    for(i=0;i<mc_nbonds[nt][1];i++){
      if(!(fscanf(mc_bond_file, "%d", &cbond))){
	fprintf(stderr, "Invalid angle index for particle %d.\n", nt);
	exit(ERR_INPUT);
      }
      mc_anglelist[nt][i]=cbond;
      
      if(!(fscanf(mc_bond_file, "%d", &cbond))){
	fprintf(stderr, "Invalid angle center index for particle %d.\n", nt);
	exit(ERR_INPUT);
      }
      mc_anglecenter[nt][i]=cbond;
      
      if(!(fscanf(mc_bond_file, "%d", &cbond))){
	fprintf(stderr, "Invalid angle type for particle %d.\n", nt);
	exit(ERR_INPUT);
      }
      mc_anglelist[nt][i+mc_nbonds[nt][1]]=cbond;
    
    }
  }
#ifdef STWCBONDED
  if(mc_nbonds[nt][2]>0){
    printf("Part %d has %d stacking bonds\n", nt, mc_nbonds[nt][2]);
    mc_stbondlist[nt]=(int *)malloc(mc_nbonds[nt][2]*sizeof(int));
    for(i=0;i<mc_nbonds[nt][2];i++){
      if(!(fscanf(mc_bond_file, "%d", &cbond))){
	fprintf(stderr, "Invalid bond index for particle %d.\n", nt);
	exit(ERR_INPUT);
      }
      mc_stbondlist[nt][i]=cbond;
      
      if(!(fscanf(mc_bond_file, "%d", &cbond))){
	fprintf(stderr, "Invalid bond type for particle %d.\n", nt);
	exit(ERR_INPUT);
      }
      if(cbond!=2){
	fprintf(stderr, "Invalid bond type index for particle %d - It should be 2 for stacking.\n", nt);
	exit(ERR_INPUT);
      }
    }
  }
 
  if(mc_nbonds[nt][3]>0){
    mc_wcbondlist[nt]=(int *)malloc(mc_nbonds[nt][3]*sizeof(int));
    for(i=0;i<mc_nbonds[nt][3];i++){
      if(!(fscanf(mc_bond_file, "%d", &cbond))){
	fprintf(stderr, "Invalid bond index for particle %d.\n", nt);
	exit(ERR_INPUT);
      }
      mc_wcbondlist[nt][i]=cbond;
     if(!(fscanf(mc_bond_file, "%d", &cbond))){
	fprintf(stderr, "Invalid bond type for particle %d.\n", nt);
	exit(ERR_INPUT);
      }
      if(cbond!=3){
	fprintf(stderr, "Invalid bond type index for particle %d - It should be 3 for watson-crick.\n", nt);
	exit(ERR_INPUT);
      }
      
    }
  }
#endif
#ifdef DIHEDRALS
  /* dihedrals */
  //printf("particle %d has %d bonds, %d angles and %d dihedrals\t NBONDDEDINTERACYTIONS IS %d\n", p, mc_nbonds[p][0], mc_nbonds[p][1], mc_nbonds[p][2], (int)N_BONDED_INTERACTIONS);fflush(stdout);
  
  if(mc_nbonds[nt][2]>0){
    mc_dihedrallist[nt]=(int *)malloc(4*mc_nbonds[nt][2]*sizeof(int));
    //mc_dihedral2[p]=(int *)malloc(mc_nbonds[p][2]*sizeof(int));
    //mc_dihedral3[p]=(int *)malloc(mc_nbonds[p][2]*sizeof(int));
    for(i=0;i<mc_nbonds[nt][2];i++){
      if(!(fscanf(mc_bond_file, "%d", &cbond))){
	fprintf(stderr, "Invalid dihedral index for particle %d.\n", nt);
	exit(ERR_INPUT);
      }
      mc_dihedrallist[nt][3*i]=cbond;
      
      if(!(fscanf(mc_bond_file, "%d", &cbond))){
	fprintf(stderr, "Invalid dihedral atom 2 index for particle %d.\n", nt);
	exit(ERR_INPUT);
      }
      mc_dihedrallist[nt][3*i+1]=cbond;
      
      if(!(fscanf(mc_bond_file, "%d", &cbond))){
	fprintf(stderr, "Invalid dihedral atom 3 index for particle %d.\n", nt);
	exit(ERR_INPUT);
      }
      mc_dihedrallist[nt][3*i+2]=cbond;
      
      if(!(fscanf(mc_bond_file, "%d", &cbond))){
	fprintf(stderr, "Invalid dihedral type for particle %d.\n", nt);
	exit(ERR_INPUT);
      }
      mc_dihedrallist[nt][i+3*mc_nbonds[nt][2]]=cbond;
    
    }
  }
#endif
}

void MC_default_bonds(int nt_n){
  //int nt_n=mc_n/N_PARTS_PER_NT;
  int i;
  for(i=0;i<nt_n;i++){
    mc_nbonds[i][0]=0;
    mc_nbonds[i][1]=0;
    
    mc_bondlist[i]=NULL;
    mc_anglelist[i]=NULL;
#ifdef DIHEDRALS
    mc_dihedrallist[i]=NULL;
#endif
  }  
}

void MC_close_bondlist_file(){
    fclose(mc_bond_file);
}
