#include "mc_checkpoints.h"

void MC_initialize_save_configs(int mc_n, int mpi_id){
  int i;
  char confname[256];
  mkdir("configs",0777);
  sprintf(confname, "configs/confs.p%02d.mc", mpi_id);
  printf("Saving configurations to file %s.\n", confname);
  if((mc_configs=fopen(confname,"wb"))!=NULL){
    //the simulation begins! write the basic data
    fwrite(&mc_n, sizeof(int), 1, mc_configs);
    for(i=0;i<mc_n;i++)
      fwrite(&(mc_types[i]), sizeof(int), 1, mc_configs);
  }
  if(PDB_OUTPUT==1){
    char confname2[256];
    sprintf(confname2,"configs/confs.p%02d.pdb",mpi_id);
    printf("WARNING: WRITING TRAJECTORY IN PDB FORMAT.\n");
    if((pdb_configs=fopen(confname2,"w"))!=NULL){
      fprintf(pdb_configs,"REMARK TRAJECTORY\n");
    }
  }
}

void MC_close_configs(){
  fclose(mc_configs);
  if(PDB_OUTPUT==1)
    fclose(pdb_configs);
}


void MC_append_pdb(int mc_n, double *rx, double *ry, double *rz, double energy_t, int iter){
  int at, nt_n=mc_n/N_PARTS_PER_NT, nt, tchain=0, nintdig;
  char at_num[10], at_name[10], res_name[10], res_num[10], chain_id[10], xp[256], yp[256], zp[256], gly[10], pck[10], glp[10], frz[10];
  char pdbname[256];
  fprintf(pdb_configs, "MODEL     %d\n", iter);
  fprintf(pdb_configs, "REMARK ENERGY %lf\n", energy_t);
  for(at=0;at<mc_n;at++){
    nt=(int)(at/N_PARTS_PER_NT);
    sprintf(at_num, "%5d", at);
    if(mc_types[at]==-1) sprintf(at_name,"XVEC");
    else if(mc_types[at]==-2) sprintf(at_name,"YVEC");
    else if(mc_types[at]==-3) sprintf(at_name,"SUGR");
    else if(mc_types[at]==-4) sprintf(at_name,"PHOS");
    else sprintf(at_name,"BASE");
    
    if(mc_types[nt*N_PARTS_PER_NT]==TYP_ADENINE)  sprintf(res_name, "  A");
    else if(mc_types[nt*N_PARTS_PER_NT]==TYP_URACIL)   sprintf(res_name, "  U");
    else if(mc_types[nt*N_PARTS_PER_NT]==TYP_GUANINE)  sprintf(res_name, "  G");
    else if(mc_types[nt*N_PARTS_PER_NT]==TYP_CYTOSINE) sprintf(res_name, "  C");
    else {printf("Residue type not recognized writing pdb\n"); exit(ERR_WRITING);}
    if(nt>0 && at%N_PARTS_PER_NT==0) {if(MC_are_neighbors(nt, nt-1)==0) tchain++;}
    tchain=0;
    sprintf(chain_id, "%d", tchain);
    sprintf(res_num, "%4d", nt);
      
    sprintf(xp, "%d", (int)rx[at]);
    nintdig=strlen(xp); if(nintdig>8){printf("ERROR: particle %d out of range for x coordinate in pdb format\n", at); exit(ERR_WRITING);}if(rx[at]<0 && (int)rx[at]==0) nintdig++;
    //sprintf(xp, "%.*f ", rx[at], 6-nintdig);
    sprintf(xp, "%.*f ", 6-nintdig,rx[at]);
    sprintf(yp, "%d", (int)ry[at]);
    nintdig=strlen(yp); if(nintdig>8){printf("ERROR: particle %d out of range for y coordinate in pdb format\n", at); exit(ERR_WRITING);}if(ry[at]<0 && (int)ry[at]==0) nintdig++;
    //sprintf(yp, "%.*f ", ry[at], 6-nintdig);
    sprintf(yp, "%.*f ", 6-nintdig,ry[at]);
    sprintf(zp, "%d", (int)rz[at]);
    nintdig=strlen(zp); if(nintdig>8){printf("ERROR: particle %d out of range for z coordinate in pdb format\n", at); exit(ERR_WRITING);}if(rz[at]<0 && (int)rz[at]==0) nintdig++;
    //sprintf(zp, "%.*f ", rz[at], 6-nintdig);
    sprintf(zp, "%.*f ", 6-nintdig,rz[at]);
    /* sprintf(xp, "%8.3f", rx[at]); */
    /* sprintf(yp, "%8.3f", ry[at]); */
    /* sprintf(zp, "%8.3f", rz[at]); */
    
    fprintf(pdb_configs, "ATOM  %s %s %s %s%s    %s%s%s", at_num, at_name, res_name, chain_id, res_num, xp, yp, zp);
    //printf( "ATOM  %s     %d  %d   %d", res_name,  at, nt, mc_types[nt*N_PARTS_PER_NT]);
    if(at%N_PARTS_PER_NT==0){
      if(mc_glyc[nt]==GLYC_A) sprintf(gly,"A");
      if(mc_glyc[nt]==GLYC_H) sprintf(gly,"H");
      if(mc_glyc[nt]==GLYC_S) sprintf(gly,"S");
      if(mc_puck[nt]==PUCK_3) sprintf(pck,"3");
      if(mc_puck[nt]==PUCK_2) sprintf(pck,"2");
      if(glp_is_flippable[nt]==GLP_FIXED) sprintf(glp,"N");
      if(glp_is_flippable[nt]==GLP_GLYC) sprintf(glp,"G");
      if(glp_is_flippable[nt]==GLP_PUCK) sprintf(glp,"P");
      if(glp_is_flippable[nt]==GLP_BOTH) sprintf(glp,"A");
      if(fr_is_mobile[nt]==FR_MOB_FROZ) sprintf(frz,"N");
      if(fr_is_mobile[nt]==FR_MOB_BASE) sprintf(frz,"B");
      if(fr_is_mobile[nt]==FR_MOB_PHOS) sprintf(frz,"P");
      if(fr_is_mobile[nt]==FR_MOB_FULL) sprintf(frz,"A");
      fprintf(pdb_configs, "  %s%s%s%s", gly, pck, glp, frz);
    }
    fprintf(pdb_configs, "\n");
  }
  fprintf(pdb_configs,"ENDMDL\n");
}


void MC_write_pdb(char * name, int mc_n, double *rx, double *ry, double *rz, double energy_t, int mpi_id){
  int at, nt_n=mc_n/N_PARTS_PER_NT, nt, tchain=0, nintdig;
  FILE *pdbfile;
  char at_num[10], at_name[10], res_name[10], res_num[10], chain_id[10], xp[256], yp[256], zp[256], gly[10], pck[10], glp[10], frz[10];
  char pdbname[256];
  sprintf(pdbname, "%s.p%02d.pdb", name, mpi_id);
  if((pdbfile=fopen(pdbname, "w"))==NULL){
    printf("Unable to create pdb file.\n");
    exit(ERR_WRITING);
  }
  else{
    fprintf(pdbfile, "REMARK ENERGY %lf\n", energy_t);
    for(at=0;at<mc_n;at++){
      nt=(int)(at/N_PARTS_PER_NT);
      sprintf(at_num, "%5d", at);
      if(mc_types[at]==-1) sprintf(at_name,"XVEC");
      else if(mc_types[at]==-2) sprintf(at_name,"YVEC");
      else if(mc_types[at]==-3) sprintf(at_name,"SUGR");
      else if(mc_types[at]==-4) sprintf(at_name,"PHOS");
      else sprintf(at_name,"BASE");
      
      if(mc_types[nt*N_PARTS_PER_NT]==TYP_ADENINE)  sprintf(res_name, "  A");
      else if(mc_types[nt*N_PARTS_PER_NT]==TYP_URACIL)   sprintf(res_name, "  U");
      else if(mc_types[nt*N_PARTS_PER_NT]==TYP_GUANINE)  sprintf(res_name, "  G");
      else if(mc_types[nt*N_PARTS_PER_NT]==TYP_CYTOSINE) sprintf(res_name, "  C");
      else {printf("Residue type not recognized writing pdb\n"); exit(ERR_WRITING);}
      if(nt>0 && at%N_PARTS_PER_NT==0) {if(MC_are_neighbors(nt, nt-1)==0) tchain++;}
      //tchain=0;
      sprintf(chain_id, "%d", tchain);
      sprintf(res_num, "%4d", nt);
      
      sprintf(xp, "%d", (int)rx[at]);
      nintdig=strlen(xp); if(nintdig>8){printf("ERROR: particle %d out of range for x coordinate in pdb format\n", at); exit(ERR_WRITING);}if(rx[at]<0 && (int)rx[at]==0) nintdig++;
      //sprintf(xp, "%.*f ", rx[at], 6-nintdig);
      sprintf(xp, "%.*f ", 6-nintdig,rx[at]);
      sprintf(yp, "%d", (int)ry[at]);
      nintdig=strlen(yp); if(nintdig>8){printf("ERROR: particle %d out of range for y coordinate in pdb format\n", at); exit(ERR_WRITING);}if(ry[at]<0 && (int)ry[at]==0) nintdig++;
      //sprintf(yp, "%.*f ", ry[at], 6-nintdig);
      sprintf(yp, "%.*f ", 6-nintdig,ry[at]);
      sprintf(zp, "%d", (int)rz[at]);
      nintdig=strlen(zp); if(nintdig>8){printf("ERROR: particle %d out of range for z coordinate in pdb format\n", at); exit(ERR_WRITING);}if(rz[at]<0 && (int)rz[at]==0) nintdig++;
      //sprintf(zp, "%.*f ", rz[at], 6-nintdig);
      sprintf(zp, "%.*f ", 6-nintdig,rz[at]);
      /* sprintf(xp, "%8.3f", rx[at]); */
      /* sprintf(yp, "%8.3f", ry[at]); */
      /* sprintf(zp, "%8.3f", rz[at]); */
      
      fprintf(pdbfile, "ATOM  %s %s %s %s%s    %s%s%s", at_num, at_name, res_name, chain_id, res_num, xp, yp, zp);
      //printf( "ATOM  %s     %d  %d   %d", res_name,  at, nt, mc_types[nt*N_PARTS_PER_NT]);
      if(at%N_PARTS_PER_NT==0){
      	if(mc_glyc[nt]==GLYC_A) sprintf(gly,"A");
      	if(mc_glyc[nt]==GLYC_H) sprintf(gly,"H");
      	if(mc_glyc[nt]==GLYC_S) sprintf(gly,"S");
      	if(mc_puck[nt]==PUCK_3) sprintf(pck,"3");
      	if(mc_puck[nt]==PUCK_2) sprintf(pck,"2");
      	if(glp_is_flippable[nt]==GLP_FIXED) sprintf(glp,"N");
      	if(glp_is_flippable[nt]==GLP_GLYC) sprintf(glp,"G");
      	if(glp_is_flippable[nt]==GLP_PUCK) sprintf(glp,"P");
      	if(glp_is_flippable[nt]==GLP_BOTH) sprintf(glp,"A");
      	if(fr_is_mobile[nt]==FR_MOB_FROZ) sprintf(frz,"N");
      	if(fr_is_mobile[nt]==FR_MOB_BASE) sprintf(frz,"B");
      	if(fr_is_mobile[nt]==FR_MOB_PHOS) sprintf(frz,"P");
      	if(fr_is_mobile[nt]==FR_MOB_FULL) sprintf(frz,"A");
      	fprintf(pdbfile, "  %s%s%s%s", gly, pck, glp, frz);
      }
      fprintf(pdbfile, "\n");
    }
  }
  fclose(pdbfile);
}

void MC_save_configuration(int mc_n, double *rx, double *ry, double *rz, double energy_t, int iter){
  if(PDB_OUTPUT==1)
    MC_append_pdb(mc_n, rx, ry, rz, energy_t,iter);
  int i, tempi, nt;
  double tempf;
  int nt_n=mc_n/N_PARTS_PER_NT;
  
  //we start saving the energy and temperature
  tempf=(double)(energy_t);
  fwrite(&tempf,sizeof(double),1,mc_configs);
  fwrite(&mc_target_temp,sizeof(double),1,mc_configs);
  //and then the configuration
  for(i=0;i<mc_n;i++){
    //nt=i/N_PARTS_PER_NT;
    tempf=(double)(get_unf_coo_x(rx,i));
    //rx[i]+mc_pbox[i][0]*box_l[0]);
    fwrite(&tempf,sizeof(double),1,mc_configs);
    tempf=(double)(get_unf_coo_y(ry,i));
    fwrite(&tempf,sizeof(double),1,mc_configs);
    tempf=(double)(get_unf_coo_z(rz,i));
    fwrite(&tempf,sizeof(double),1,mc_configs);
  }
  for(nt=0;nt<nt_n;nt++){
    tempi=(int)(mc_glyc[nt]);
    fwrite(&tempi,sizeof(int),1,mc_configs);
    tempi=(int)(mc_puck[nt]);
    fwrite(&tempi,sizeof(int),1,mc_configs);
  }
}

int MC_read_checkpoint(int *mc_n, double **rx, double **ry, double **rz, int *rand_a, int mpi_id, char *chkname, int read_flag, double * add_data){
  //read binary configuration file, as saved in a checkpoint - read types, positions, glp, frz, and bonds
  FILE *chkfile;
  int nt_n;
  int i, d, tempi, *chains, tempiter=0;
  double tempd;
  double tx, ty, tz;
  size_t frout;
  long tidum2, tiy, tiv[NTAB];
  
  if((chkfile=fopen(chkname, "rb"))==NULL){
    printf("Unable to read binary file %s\n", chkname);
  }
  else{
    frout=fread(mc_n, sizeof(int), 1, chkfile);
    nt_n=*mc_n/N_PARTS_PER_NT;
    chains=(int *)malloc(nt_n*sizeof(int));
    printf("File %s open successfully with %d nucleotides.\n", chkname, nt_n);
    MC_initialize_global(*mc_n, *rand_a,mpi_id); //here is initialized the random seed
    MC_initialize_arrays(*mc_n, rx, ry, rz); //here the mpi_id is used only for printing purposes and initializing the random seed
    for(i=0;i<*mc_n;i++){
      frout=fread(&tempi, sizeof(int),1,chkfile);
      mc_types[i]=tempi;
    }
    frout=fread(&tempd, sizeof(double),1,chkfile); //energy
    frout=fread(&tempd, sizeof(double),1,chkfile); //temperature
    if(read_flag==1)
      mc_target_temp=tempd;
    
    for(i=0;i<*mc_n;i++){
      frout=fread(&tempd,sizeof(double),1,chkfile);(*rx)[i]=tempd;//tx=tempd;
      frout=fread(&tempd,sizeof(double),1,chkfile);(*ry)[i]=tempd;//ty=tempd;
      frout=fread(&tempd,sizeof(double),1,chkfile);(*rz)[i]=tempd;//tz=tempd;
      //printf("%d  %lf %lf %lf\n", i, tx, ty, tz);
    }
    for(i=0;i<nt_n;i++){
      frout=fread(&tempi,sizeof(int),1,chkfile);mc_glyc[i]=tempi;
      frout=fread(&tempi,sizeof(int),1,chkfile);mc_puck[i]=tempi;
      frout=fread(&tempi,sizeof(int),1,chkfile);glp_is_flippable[i]=tempi;
      frout=fread(&tempi,sizeof(int),1,chkfile);fr_is_mobile[i]=tempi;
      frout=fread(&tempi,sizeof(int),1,chkfile);chains[i]=tempi;
      //printf("%d  %d %d   %d    %d %d\n", i, mc_glyc[i], mc_puck[i], chains[i], glp_is_flippable[i], fr_is_mobile[i]);
    }
    frout=fread(&tempd,sizeof(double),1,chkfile);//idum
    //fread(idum,sizeof(double),1,chkpnt);
    frout=fread(&tidum2, sizeof(long),1,chkfile);
    frout=fread(&tiy, sizeof(long),1,chkfile);
    frout=fread(tiv, sizeof(long),NTAB,chkfile);
    
    if(read_flag==1){
      *idum=tempd;
      idum2=tidum2;
      iy=tiy;
      for(i=0;i<NTAB;i++)
	iv[i]=tiv[i];
    } else{
      //just for safety
      idum2=123456789;
      iy=0;
    }
    
    frout=fread(&tempiter,sizeof(int),1,chkfile);
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
  int nt=0;
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
    /* printf("bonds %d : %d  ", i, mc_nbonds[i][0]); */
    /* for(nt=0;nt<mc_nbonds[i][0];nt++) */
    /*   printf("%d ", mc_bondlist[i][nt]); */
    /* printf("\n");*/
  }
  
#ifdef TANN
//in this case, add_data is the annealing parameters
  char word[3];
  int tint;
  if(fread(word, sizeof(char),3,chkfile)==3 && (word[0]=='A' && word[1]=='N' && word[2]=='N')){
    frout=fread(&tint, sizeof(int),1,chkfile);
    //add_data must already be allocated!
    frout=fread(add_data, sizeof(double),ANN_NPARAMS,chkfile);  
    if(tint!=ANN_NPARAMS){
      printf("Number of parameters of annealing (%d) does not match the number of parameters of checkpoint file (%d).\n", tint, (int)ANN_NPARAMS);
      exit(ERR_INPUT);
    }
  }
#endif
  return tempiter;
}

void MC_save_checkpoint(int mc_n, double *rx, double *ry, double *rz, int iter, double energy_t,int mpi_id, double *add_data){
  int i, tempi, nt;
  double tempf;
  int nt_n=mc_n/N_PARTS_PER_NT;
  char confname[256];
  FILE *chkpnt;
  char word[3]="END";
  int tchain=0;
  if(iter<0) sprintf(confname, "configs/chk.last.p%02d.mc", mpi_id);
  else sprintf(confname, "configs/chk.%010d.p%02d.mc", iter, mpi_id);
  //printf("Saving configurations to file %s.\n", confname);
  if((chkpnt=fopen(confname,"wb"))!=NULL){
    fwrite(&mc_n, sizeof(int), 1, chkpnt);
    for(i=0;i<mc_n;i++)
      fwrite(&(mc_types[i]), sizeof(int), 1, chkpnt);
    tempf=(double)(energy_t);
    fwrite(&tempf,sizeof(double),1,chkpnt);
    fwrite(&mc_target_temp,sizeof(double),1,chkpnt); 
    //and then the configuration
    for(i=0;i<mc_n;i++){
      tempf=(double)(get_unf_coo_x(rx,i));
      fwrite(&tempf,sizeof(double),1,chkpnt);
      tempf=(double)(get_unf_coo_y(ry,i));
      fwrite(&tempf,sizeof(double),1,chkpnt);
      tempf=(double)(get_unf_coo_z(rz,i));
      fwrite(&tempf,sizeof(double),1,chkpnt);
    }
    for(nt=0;nt<nt_n;nt++){
      if(nt>0) {if(MC_are_neighbors(nt, nt-1)==0) tchain++;}
      tempi=(int)(mc_glyc[nt]);
      fwrite(&tempi,sizeof(int),1,chkpnt);
      tempi=(int)(mc_puck[nt]);
      fwrite(&tempi,sizeof(int),1,chkpnt);
      tempi=(int)(glp_is_flippable[nt]);
      fwrite(&tempi,sizeof(int),1,chkpnt);
      tempi=(int)(fr_is_mobile[nt]);
      fwrite(&tempi,sizeof(int),1,chkpnt);
      
      fwrite(&tchain,sizeof(int),1,chkpnt);
    }
    //finally, the random seed
    fwrite(idum,sizeof(double),1,chkpnt);
    fwrite(&idum2, sizeof(long),1,chkpnt);
    fwrite(&iy, sizeof(long),1,chkpnt);
    fwrite(iv, sizeof(long),NTAB,chkpnt);
    //and the time step
    int tempiter=iter;
    if(iter<0) tempiter=0;
    fwrite(&tempiter,sizeof(int),1,chkpnt);
   

#ifdef TANN
    //in this case, add_data is the annealing parameters
    char word[3]="ANN";
    int tint=(int)ANN_NPARAMS;
    fwrite(word, sizeof(char),3,chkpnt);
    fwrite(&tint, sizeof(int),1,chkpnt);
    fwrite(add_data, sizeof(double),ANN_NPARAMS,chkpnt);
#endif
    fclose(chkpnt);
  }
  else{
    printf("Unable to write checkpoint!\n");
    exit(ERR_WRITING);
  }
}


void MC_save_xyz_configuration(int mc_n, double *rx, double *ry, double *rz, double final_iter, double energy_t, int mpi_id){
  FILE * mcconf;
  char filename[256];
  sprintf(filename, "configs/lastmc.p%02d.xyz", mpi_id);
  int i;

  if(!(mcconf=fopen(filename , "w"))){
    fprintf(stderr, "Unable to write final MC configuration file.\n");
    //exit(ERR_INPUT);
  } else {
    //double unfolded_pos[DIM];
    fprintf(mcconf, "%d\n", mc_n);
    fprintf(mcconf, "Energy %lf Temperature %lf Step %f\n", energy_t, mc_target_temp, final_iter);
    for(i=0;i<mc_n;i++){
      //unfold_position(mc_particles[i].pos,mc_particles[i].box,unfolded_pos);
      fprintf(mcconf, "%d\t", mc_types[i]);
      fprintf(mcconf, "%f\t", get_unf_coo_x(rx,i));
      fprintf(mcconf, "%f\t", get_unf_coo_y(ry,i));
      fprintf(mcconf, "%f\t", get_unf_coo_z(rz,i));
      fprintf(mcconf, "\n");
    }
    fclose(mcconf);
  }

    //we also save gly_pck configuration
  FILE * glpconf;
  char glpname[256];
  sprintf(glpname, "configs/gly_pck.last.p%02d.dat", mpi_id);
  if(!(glpconf=fopen(glpname, "w"))){
    fprintf(stderr, "Unable to write current glp configuration file.\n");
    exit(ERR_INPUT);
  }else{
    char cglyc;
    int puck, glpflag;
    for(i=0;i<mc_n/N_PARTS_PER_NT;i++){
      if(mc_puck[i]==PUCK_2) puck=2;
      else puck=3;
      if(mc_glyc[i]==GLYC_A) cglyc='A';
      else if(mc_glyc[i]==GLYC_H) cglyc='H';
      else cglyc='S';
      
      fprintf(glpconf, "%d %c %d %d\n", i, cglyc, puck, glp_is_flippable[i]);
      
    }
    fclose(glpconf);
  }

}

void MC_save_current_configuration(int mc_n, double *rx, double *ry, double *rz, double iter, int index, double energy_t, int mpi_id){
  FILE * mcconf;
  char filename[256];
  int i;
  sprintf(filename, "configs/cmc.%010d.p%02d.xyz",index, mpi_id);
  //  char filename[11]="lastmc.xyz";
  
  if(!(mcconf=fopen(filename , "w"))){
    fprintf(stderr, "Unable to write current MC configuration file.\n");
    exit(ERR_INPUT);
  } else {
    //double unfolded_pos[DIM];
    fprintf(mcconf, "%d\n", mc_n);
    fprintf(mcconf, "Energy %lf  Temperature %lf Step %f\n", energy_t, mc_target_temp, iter);
    for(i=0;i<mc_n;i++){
      //unfold_position(mc_particles[i].pos,mc_particles[i].box,unfolded_pos);
      fprintf(mcconf, "%d\t", mc_types[i]);
      fprintf(mcconf, "%f\t", get_unf_coo_x(rx,i));
      //rx[i]+mc_pbox[i][0]*box_l[0]);
      fprintf(mcconf, "%f\t", get_unf_coo_y(ry,i));
      fprintf(mcconf, "%f\t", get_unf_coo_z(rz,i));
      fprintf(mcconf, "\n");
    }
    fclose(mcconf);
  }
  
  //we also save gly_pck configuration
  FILE * glpconf;
  char glpname[256];
  sprintf(glpname, "configs/gly_pck.%010d.p%02d.dat", index, mpi_id);
  if(!(glpconf=fopen(glpname, "w"))){
    fprintf(stderr, "Unable to write current glp configuration file.\n");
    exit(ERR_INPUT);
  }else{
    char cglyc;
    int puck, glpflag;
    for(i=0;i<mc_n/N_PARTS_PER_NT;i++){
      if(mc_puck[i]==PUCK_2) puck=2;
      else puck=3;
      if(mc_glyc[i]==GLYC_A) cglyc='A';
      else if(mc_glyc[i]==GLYC_H) cglyc='H';
      else cglyc='S';
      
      fprintf(glpconf, "%d %c %d %d\n", i, cglyc, puck, glp_is_flippable[i]);
      
    }
    fclose(glpconf);
  }
}


void MC_min_energ_xyz_configuration(int mc_n, double *rx, double *ry, double *rz, double energ, double temperature, int init){
  FILE * mcconf;
  char filename[256];
  sprintf(filename, "configs/min_energ.p%02d.xyz", init);
  
  if(!(mcconf=fopen(filename , "w"))){
    fprintf(stderr, "Unable to write minimum energy MC configuration file.\n");
    //exit(ERR_INPUT);
  } else {
    int i;
    //double unfolded_pos[DIM];
    fprintf(mcconf, "%d\n", mc_n);
    fprintf(mcconf, "%lf  %lf\n", energ, temperature);
    for(i=0;i<mc_n;i++){
      //unfold_position(mc_particles[i].pos,mc_particles[i].box,unfolded_pos);
      fprintf(mcconf, "%d\t", mc_types[i]);
      fprintf(mcconf, "%f\t",  get_unf_coo_x(rx,i));
      //rx[i]+mc_pbox[i][0]*box_l[0]);
      fprintf(mcconf, "%f\t",  get_unf_coo_y(ry,i));
      fprintf(mcconf, "%f\t",  get_unf_coo_z(rz,i));
      fprintf(mcconf, "\n");
    }
    fclose(mcconf);
  }
}
