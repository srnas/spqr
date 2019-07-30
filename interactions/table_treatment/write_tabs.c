#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mc.h"
#include "mc_energies.h"

void MC_read_write_energy_tables(){
  size_t binnum=sizeof(double);
  mc_n_types=N_BASES;
  FILE *outfile;
  char tablename[256];
  char outname[256];
  int i,j, typ_ind, ener, ntot, stf;
  //int ntypsq=mc_n_types*mc_n_types;
  int ntypsq=N_BASES*N_BASES;
  char ctemp1, ctemp2, ctemp3, ctemp4, ctemp5;
  sprintf(outname, "intrac.btb");
  if((outfile=fopen(outname, "wb"))==NULL){
    printf("ERROR: Could not open %s file for writing tables.\n", outname);
    exit(ERR_WRITING);
  }
  
  /* NON BONDED POTENTIAL WELLS */
  sprintf(tablename, "tab_energs/table_wells.tab");
  FILE *table;
  if((table=fopen(tablename,"r"))==NULL){
    printf("file not found!\n");
    exit(1);
  }
  else{
    fscanf(table,"%c%c", &ctemp1, &ctemp2);
    if(ctemp1!='b' || ctemp2 != 'b'){
      printf("Wrong syntax at table_wells.tab file (bb)!\n%c %c\n", ctemp1, ctemp2);
      exit(1);}
    fscanf(table, "%lf", &(BB_PREF));
    fwrite(&BB_PREF, binnum, 1, outfile); 
    fscanf(table, "%lf", &(BB_PREF_A));
    fwrite(&BB_PREF_A, binnum, 1, outfile);
    fscanf(table,"%c%c%c", &ctemp1, &ctemp2, &ctemp3);
    if(ctemp2!='g' || ctemp3 != 'l'){
      printf("Wrong syntax at table_wells.tab file (gl)!\n%c %c\n", ctemp1, ctemp2);
      exit(1);
      
    }
    for(j=0;j<N_PUCK_STATES;j++)
      for(i=0;i<N_GLYC_STATES;i++){
	fscanf(table, "%lf", &(glp_well_R[i][j]));
	fwrite(&(glp_well_R[i][j]), binnum, 1, outfile);
      }
    for(j=0;j<N_PUCK_STATES;j++)
      for(i=0;i<N_GLYC_STATES;i++){
	fscanf(table, "%lf", &(glp_well_Y[i][j]));
	fwrite(&(glp_well_Y[i][j]), binnum, 1, outfile);
      }    
    
    /* NON BONDED POTENTIAL WELLS */
    fscanf(table,"%c%c%c%c", &ctemp1, &ctemp2, &ctemp3, &ctemp4);
    if(ctemp2!='s' || ctemp3 != 't' || ctemp4!='b'){
      printf("Wrong syntax at table_wells.tab file (st , 1)!\n%c %c %c %c\n", ctemp1, ctemp2, ctemp3, ctemp4);
      exit(1);}
    for(stf=0;stf<N_STFACES;stf++){
      for(i=0;i<N_BASES;i++)
	for(j=0;j<N_BASES;j++){
	  fscanf(table, "%lf", &(b_st_well[N_BASES*i+j][stf]));
	  fwrite(&(b_st_well[N_BASES*i+j][stf]), binnum, 1, outfile);
	}
    }
    fscanf(table,"%c%c%c",&ctemp3, &ctemp1, &ctemp2);
    if(ctemp1!='s' || ctemp2 != 't'){
      printf("Wrong syntax at table_wells.tab file (st , 2)!\n%c %c %c", ctemp1, ctemp2, ctemp3);
      exit(1);}
    for(i=0;i<N_BASES;i++){
      for(j=0;j<N_BASES;j++){
	fscanf(table, "%lf", &(nb_st_well[N_BASES*i+j]));
	fwrite(&(nb_st_well[N_BASES*i+j]), binnum, 1, outfile);
      }
    }
    
    for(i=0;i<N_BASES;i++)
      for(j=0;j<N_BASES;j++){
	if(nb_st_well[N_BASES*i+j]!=nb_st_well[N_BASES*j+i]){
	  printf("Non-bonded matrix (stacking) is not symmetric!\n");
	  exit(1);
	}
      }
    fscanf(table,"%c%c%c", &ctemp3,&ctemp1, &ctemp2);
    if(ctemp1!='w' || ctemp2 != 'c'){
      printf("Wrong syntax at table_wells.tab file (wc)!\n %c %c %c", ctemp1, ctemp2, ctemp3);
      exit(1);}
    for(i=0;i<N_BASES;i++)
      for(j=0;j<N_BASES;j++){
	fscanf(table,"%lf", &(nb_wc_well[N_BASES*i+j][0*WC_FACES+0]));
	fwrite(&(nb_wc_well[N_BASES*i+j][0*WC_FACES+0]), binnum, 1, outfile);
	fscanf(table,"%lf", &(nb_wc_well[N_BASES*i+j][1*WC_FACES+1]));
	fwrite(&(nb_wc_well[N_BASES*i+j][1*WC_FACES+1]), binnum, 1, outfile);
	fscanf(table,"%lf", &(nb_wc_well[N_BASES*i+j][2*WC_FACES+2]));
	fwrite(&(nb_wc_well[N_BASES*i+j][2*WC_FACES+2]), binnum, 1, outfile);
	
	fscanf(table,"%lf", &(nb_wc_well[N_BASES*i+j][1*WC_FACES+2]));
	fwrite(&(nb_wc_well[N_BASES*i+j][1*WC_FACES+2]), binnum, 1, outfile);
	nb_wc_well[N_BASES*j+i][2*WC_FACES+1]=nb_wc_well[N_BASES*i+j][1*WC_FACES+2];
	fscanf(table,"%lf", &(nb_wc_well[N_BASES*i+j][1*WC_FACES+0]));
	fwrite(&(nb_wc_well[N_BASES*i+j][1*WC_FACES+0]), binnum, 1, outfile);
	nb_wc_well[N_BASES*j+i][0*WC_FACES+1]=nb_wc_well[N_BASES*i+j][1*WC_FACES+0];
	fscanf(table,"%lf", &(nb_wc_well[N_BASES*i+j][2*WC_FACES+0]));
	fwrite(&(nb_wc_well[N_BASES*i+j][2*WC_FACES+0]), binnum, 1, outfile);
	nb_wc_well[N_BASES*j+i][0*WC_FACES+2]=nb_wc_well[N_BASES*i+j][2*WC_FACES+0];
      }
    
    //BASE-PAIRING IN PARALLEL CONFORMATION
    fscanf(table,"%c%c%c", &ctemp3,&ctemp1, &ctemp2);
    if(ctemp1!='w' || ctemp2 != 'P' ){
      printf("Wrong syntax at table_wells.tab file (wc)!\n %c %c %c", ctemp1, ctemp2, ctemp3);
      exit(1);}
    for(i=0;i<N_BASES;i++)
      for(j=0;j<N_BASES;j++){
	fscanf(table,"%lf", &(nb_wc_well_F[N_BASES*i+j][0*WC_FACES+0]));
	fwrite(&(nb_wc_well_F[N_BASES*i+j][0*WC_FACES+0]), binnum, 1, outfile);
	fscanf(table,"%lf", &(nb_wc_well_F[N_BASES*i+j][1*WC_FACES+1]));
	fwrite(&(nb_wc_well_F[N_BASES*i+j][1*WC_FACES+1]), binnum, 1, outfile);
	fscanf(table,"%lf", &(nb_wc_well_F[N_BASES*i+j][2*WC_FACES+2]));
	fwrite(&(nb_wc_well_F[N_BASES*i+j][2*WC_FACES+2]), binnum, 1, outfile);
	
	fscanf(table,"%lf", &(nb_wc_well_F[N_BASES*i+j][1*WC_FACES+2]));
	fwrite(&(nb_wc_well_F[N_BASES*i+j][1*WC_FACES+2]), binnum, 1, outfile);
	nb_wc_well_F[N_BASES*j+i][2*WC_FACES+1]=nb_wc_well_F[N_BASES*i+j][1*WC_FACES+2];
	fscanf(table,"%lf", &(nb_wc_well_F[N_BASES*i+j][1*WC_FACES+0]));
	fwrite(&(nb_wc_well_F[N_BASES*i+j][1*WC_FACES+0]), binnum, 1, outfile);
	nb_wc_well_F[N_BASES*j+i][0*WC_FACES+1]=nb_wc_well_F[N_BASES*i+j][1*WC_FACES+0];
	fscanf(table,"%lf", &(nb_wc_well_F[N_BASES*i+j][2*WC_FACES+0]));
	fwrite(&(nb_wc_well_F[N_BASES*i+j][2*WC_FACES+0]), binnum, 1, outfile);
	nb_wc_well_F[N_BASES*j+i][0*WC_FACES+2]=nb_wc_well_F[N_BASES*i+j][2*WC_FACES+0];
      }
    
    /************** BASE PHOSPHATE ***************/
    fscanf(table,"%c%c%c", &ctemp3,&ctemp1, &ctemp2);
    //printf("reading wc %c  %c\n", ctemp1, ctemp2);
    if(ctemp1!='b' || ctemp2 != 'p'){
      printf("Wrong syntax at table_wells.tab file (wc)!\n %c %c %c", ctemp1, ctemp2, ctemp3);
      exit(1);}
    for(i=0;i<N_BASES;i++){
      for(j=0;j<WC_FACES;j++){
	fscanf(table,"%lf", &(nb_bp_well[i][j]));
	fwrite(&(nb_bp_well[i][j]), binnum, 1, outfile);
	printf("%lf ", nb_bp_well[i][j]);
      }
      printf("\n");
    }
    fscanf(table, "%lf", &(nb_bp_spec_well[2][0]));
    fwrite(&(nb_bp_spec_well[2][0]), binnum, 1, outfile);
    fscanf(table, "%lf", &(nb_bp_spec_well[2][1]));
    fwrite(&(nb_bp_spec_well[2][1]), binnum, 1, outfile);
    fscanf(table, "%lf", &(nb_bp_spec_well[2][2]));
    fwrite(&(nb_bp_spec_well[2][2]), binnum, 1, outfile);
    fclose(table);
  }
  
  /************************ BASE-PAIR SECOND DIHEDRALS **************************/
  
  sprintf(tablename, "tab_energs/table_wc_secdih.tab");
  if((table=fopen(tablename,"r"))==NULL){
    printf("file not found!\n");
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
	fwrite(&(wc_secdih_min[N_BASES*i+j][0*WC_FACES+0]), binnum, 1, outfile);
	fscanf(table,"%lf", &(wc_secdih_min[N_BASES*i+j][1*WC_FACES+1]));
	fwrite(&(wc_secdih_min[N_BASES*i+j][1*WC_FACES+1]), binnum, 1, outfile);
	fscanf(table,"%lf", &(wc_secdih_min[N_BASES*i+j][2*WC_FACES+2]));
	fwrite(&(wc_secdih_min[N_BASES*i+j][2*WC_FACES+2]), binnum, 1, outfile);
	
	fscanf(table,"%lf", &(wc_secdih_min[N_BASES*i+j][1*WC_FACES+2]));
	fwrite(&(wc_secdih_min[N_BASES*i+j][1*WC_FACES+2]), binnum, 1, outfile);
	wc_secdih_min[N_BASES*j+i][2*WC_FACES+1]=wc_secdih_min[N_BASES*i+j][1*WC_FACES+2];
	fscanf(table,"%lf", &(wc_secdih_min[N_BASES*i+j][1*WC_FACES+0]));
	fwrite(&(wc_secdih_min[N_BASES*i+j][1*WC_FACES+0]), binnum, 1, outfile);
	wc_secdih_min[N_BASES*j+i][0*WC_FACES+1]=wc_secdih_min[N_BASES*i+j][1*WC_FACES+0];
	fscanf(table,"%lf", &(wc_secdih_min[N_BASES*i+j][2*WC_FACES+0]));
	fwrite(&(wc_secdih_min[N_BASES*i+j][2*WC_FACES+0]), binnum, 1, outfile);
	wc_secdih_min[N_BASES*j+i][0*WC_FACES+2]=wc_secdih_min[N_BASES*i+j][2*WC_FACES+0];
      }
    
    fscanf(table,"%c%c%c%c%c", &ctemp5,&ctemp1, &ctemp2, &ctemp3, &ctemp4);
    //printf("reading wc %c  %c\n", ctemp1, ctemp2);
    if(ctemp1!='A' || ctemp2 != 'm' || ctemp3 != 'a' || ctemp4!='x'){
      printf("Wrong syntax at table_wells.tab file (wc)!\n %c %c %c", ctemp1, ctemp2, ctemp3);
      exit(1);}
    for(i=0;i<N_BASES;i++)
      for(j=0;j<N_BASES;j++){
	fscanf(table,"%lf", &(wc_secdih_max[N_BASES*i+j][0*WC_FACES+0]));
	fwrite(&(wc_secdih_max[N_BASES*i+j][0*WC_FACES+0] ), binnum, 1, outfile);
	fscanf(table,"%lf", &(wc_secdih_max[N_BASES*i+j][1*WC_FACES+1]));
	fwrite(&(wc_secdih_max[N_BASES*i+j][1*WC_FACES+1] ), binnum, 1, outfile);
	fscanf(table,"%lf", &(wc_secdih_max[N_BASES*i+j][2*WC_FACES+2]));
	fwrite(&(wc_secdih_max[N_BASES*i+j][2*WC_FACES+2] ), binnum, 1, outfile);
	
	fscanf(table,"%lf", &(wc_secdih_max[N_BASES*i+j][1*WC_FACES+2]));
	fwrite(&(wc_secdih_max[N_BASES*i+j][1*WC_FACES+2] ), binnum, 1, outfile);
	wc_secdih_max[N_BASES*j+i][2*WC_FACES+1]=wc_secdih_max[N_BASES*i+j][1*WC_FACES+2];
	fscanf(table,"%lf", &(wc_secdih_max[N_BASES*i+j][1*WC_FACES+0]));
	fwrite(&(wc_secdih_max[N_BASES*i+j][1*WC_FACES+0] ), binnum, 1, outfile);
	wc_secdih_max[N_BASES*j+i][0*WC_FACES+1]=wc_secdih_max[N_BASES*i+j][1*WC_FACES+0];
	fscanf(table,"%lf", &(wc_secdih_max[N_BASES*i+j][2*WC_FACES+0]));
	fwrite(&(wc_secdih_max[N_BASES*i+j][2*WC_FACES+0] ), binnum, 1, outfile);
	wc_secdih_max[N_BASES*j+i][0*WC_FACES+2]=wc_secdih_max[N_BASES*i+j][2*WC_FACES+0];
      }
    
    fscanf(table,"%c%c%c%c%c", &ctemp5,&ctemp1, &ctemp2, &ctemp3, &ctemp4);
    if(ctemp1!='P' || ctemp2 != 'm' || ctemp3 != 'i' || ctemp4!='n'){
      printf("Wrong syntax at table_wells.tab file (wc)!\n %c %c %c", ctemp1, ctemp2, ctemp3);
      exit(1);}
    for(i=0;i<N_BASES;i++)
      for(j=0;j<N_BASES;j++){
	fscanf(table,"%lf", &(wc_secdih_min_F[N_BASES*i+j][0*WC_FACES+0]));
	fwrite(&(wc_secdih_min_F[N_BASES*i+j][0*WC_FACES+0] ), binnum, 1, outfile);
	fscanf(table,"%lf", &(wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+1]));
	fwrite(&(wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+1] ), binnum, 1, outfile);
	fscanf(table,"%lf", &(wc_secdih_min_F[N_BASES*i+j][2*WC_FACES+2]));
	fwrite(&(wc_secdih_min_F[N_BASES*i+j][2*WC_FACES+2] ), binnum, 1, outfile);
	
	fscanf(table,"%lf", &(wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+2]));
	fwrite(&(wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+2] ), binnum, 1, outfile);
	wc_secdih_min_F[N_BASES*j+i][2*WC_FACES+1]=wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+2];
	fscanf(table,"%lf", &(wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+0]));
	fwrite(&(wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+0] ), binnum, 1, outfile);
	wc_secdih_min_F[N_BASES*j+i][0*WC_FACES+1]=wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+0];
	fscanf(table,"%lf", &(wc_secdih_min_F[N_BASES*i+j][2*WC_FACES+0]));
	fwrite(&(wc_secdih_min_F[N_BASES*i+j][2*WC_FACES+0] ), binnum, 1, outfile);
	wc_secdih_min_F[N_BASES*j+i][0*WC_FACES+2]=wc_secdih_min_F[N_BASES*i+j][2*WC_FACES+0];
      }
    
    fscanf(table,"%c%c%c%c%c", &ctemp5,&ctemp1, &ctemp2, &ctemp3, &ctemp4);
    if(ctemp1!='P' || ctemp2 != 'm' || ctemp3 != 'a' || ctemp4!='x'){
      printf("Wrong syntax at table_wells.tab file (wc)!\n %c %c %c", ctemp1, ctemp2, ctemp3);
      exit(1);}
    for(i=0;i<N_BASES;i++)
      for(j=0;j<N_BASES;j++){
	fscanf(table,"%lf", &(wc_secdih_max_F[N_BASES*i+j][0*WC_FACES+0]));
	fwrite(&(wc_secdih_max_F[N_BASES*i+j][0*WC_FACES+0] ), binnum, 1, outfile);
	fscanf(table,"%lf", &(wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+1]));
	fwrite(&(wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+1] ), binnum, 1, outfile);
	fscanf(table,"%lf", &(wc_secdih_max_F[N_BASES*i+j][2*WC_FACES+2]));
	fwrite(&(wc_secdih_max_F[N_BASES*i+j][2*WC_FACES+2] ), binnum, 1, outfile);
	
	fscanf(table,"%lf", &(wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+2]));
	fwrite(&(wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+2] ), binnum, 1, outfile);
	wc_secdih_max_F[N_BASES*j+i][2*WC_FACES+1]=wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+2];
	fscanf(table,"%lf", &(wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+0]));
	fwrite(&(wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+0] ), binnum, 1, outfile);
	wc_secdih_max_F[N_BASES*j+i][0*WC_FACES+1]=wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+0];
	fscanf(table,"%lf", &(wc_secdih_max_F[N_BASES*i+j][2*WC_FACES+0]));
	fwrite(&(wc_secdih_max_F[N_BASES*i+j][2*WC_FACES+0] ), binnum, 1, outfile);
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
  if((table=fopen(tablename, "r"))==NULL){
    printf("No backbone interaction between sugars.\n");
    table_ssB1_N_33=0;
  }
  else{
    printf("[ANGLE P3 P3] ");
    fscanf(table,"%d", &(table_ssB1_N_33));
    fwrite(&(table_ssB1_N_33 ), sizeof(int), 1, outfile);
    fscanf(table,"%lf", &(table_ssB1_params_33[0]));
    fwrite(&(table_ssB1_params_33[0] ), binnum, 1, outfile);
    fscanf(table,"%lf", &(table_ssB1_params_33[1]));
    fwrite(&(table_ssB1_params_33[1] ), binnum, 1, outfile);
    table_ssB1_33=(double *)malloc(sizeof(double)*table_ssB1_N_33);
    for(ener=0;ener<table_ssB1_N_33;ener++){
      fscanf(table, "%lf", &(table_ssB1_33[ener]));
      fwrite(&(table_ssB1_33[ener] ), binnum, 1, outfile);
    }
    fclose(table);
  }
  
  //TYPE 32
  sprintf(tablename, "tab_energs/table_ssB1_p32.tab");
  //FILE *table;
  if((table=fopen(tablename, "r"))==NULL){
    printf("No backbone interaction between sugars.\n");
    table_ssB1_N_32=0;
    exit(ERR_INPUT);
  }
  else{
    printf("[ANGLE P3 P2] ");
    fscanf(table,"%d", &(table_ssB1_N_32));
    fwrite(&(table_ssB1_N_32 ), sizeof(int), 1, outfile);
    fscanf(table,"%lf", &(table_ssB1_params_32[0]));
    fwrite(&(table_ssB1_params_32[0] ), binnum, 1, outfile);
    fscanf(table,"%lf", &(table_ssB1_params_32[1]));
    fwrite(&(table_ssB1_params_32[1] ), binnum, 1, outfile);
    table_ssB1_32=(double *)malloc(sizeof(double)*table_ssB1_N_32);
    for(ener=0;ener<table_ssB1_N_32;ener++){
      fscanf(table, "%lf", &(table_ssB1_32[ener]));
      fwrite(&(table_ssB1_32[ener] ), binnum, 1, outfile);
    }
    fclose(table);
  }
  
  //TYPE 23
  sprintf(tablename, "tab_energs/table_ssB1_p23.tab");
  if((table=fopen(tablename, "r"))==NULL){
    printf("No backbone interaction between sugars.\n");
    table_ssB1_N_23=0;
    exit(ERR_INPUT);
  }
  else{
    printf("[ANGLE P2 P3] ");
    fscanf(table,"%d", &(table_ssB1_N_23));
    fwrite(&(table_ssB1_N_23 ), sizeof(int), 1, outfile);
    fscanf(table,"%lf", &(table_ssB1_params_23[0]));
    fwrite(&(table_ssB1_params_23[0] ), binnum, 1, outfile);
    fscanf(table,"%lf", &(table_ssB1_params_23[1]));
    fwrite(&(table_ssB1_params_23[1] ), binnum, 1, outfile);
    
    table_ssB1_23=(double *)malloc(sizeof(double)*table_ssB1_N_23);
    for(ener=0;ener<table_ssB1_N_23;ener++){
      fscanf(table, "%lf", &(table_ssB1_23[ener]));
      fwrite(&(table_ssB1_23[ener] ), binnum, 1, outfile);
    }
    fclose(table);
  }
  //TYPE 22
  sprintf(tablename, "tab_energs/table_ssB1_p22.tab");
  if((table=fopen(tablename, "r"))==NULL){
    printf("No backbone interaction between sugars.\n");
    table_ssB1_N_22=0;
    exit(ERR_INPUT);
  }
  else{
    printf("[ANGLE P2 P2] ");
    fscanf(table,"%d", &(table_ssB1_N_22));
    fwrite(&(table_ssB1_N_22 ), sizeof(int), 1, outfile);
    fscanf(table,"%lf", &(table_ssB1_params_22[0]));
    fwrite(&(table_ssB1_params_22[0] ), binnum, 1, outfile);
    fscanf(table,"%lf", &(table_ssB1_params_22[1]));
    fwrite(&(table_ssB1_params_22[1] ), binnum, 1, outfile);
    table_ssB1_22=(double *)malloc(sizeof(double)*table_ssB1_N_22);
    for(ener=0;ener<table_ssB1_N_22;ener++){
      fscanf(table, "%lf", &(table_ssB1_22[ener]));
      fwrite(&(table_ssB1_22[ener] ), binnum, 1, outfile);
    }
    fclose(table);
  }
  
  /* BASE - PHOSPHATE , INTRA NT */
  //glyc ANTI
  //puck3
  printf("\nReading INTRA-BP interactions...\n");
  table_bpI_A3=(double **)malloc(sizeof(double *)*mc_n_types);
  for(i=0;i<N_BASES;i++){
    typ_ind=i;
    sprintf(tablename, "tab_energs/table_bpI_%d_gA3.tab", i);
    if((table=fopen(tablename, "r"))==NULL){
      printf("No backbone intra BP interaction for %d, glycosidic conformation ANTI.\n", i);
      table_bpI_N_A3[typ_ind][0]=-1;
      table_bpI_N_A3[typ_ind][1]=-1;
      table_bpI_N_A3[typ_ind][2]=-1;
      exit(ERR_INPUT);
    }
    else{
      printf("[%d ANTI P3] ", i);
      fscanf(table,"%d%d%d", &(table_bpI_N_A3[typ_ind][0]), &(table_bpI_N_A3[typ_ind][1]), &(table_bpI_N_A3[typ_ind][2])); //NX, NY, NZ
      fwrite(&(table_bpI_N_A3[typ_ind][0] ), sizeof(int), 1, outfile);
      fwrite(&(table_bpI_N_A3[typ_ind][1] ), sizeof(int), 1, outfile);
      fwrite(&(table_bpI_N_A3[typ_ind][2] ), sizeof(int), 1, outfile);
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_A3[typ_ind][0][0]), &(table_bpI_params_A3[typ_ind][1][0]), &(table_bpI_params_A3[typ_ind][2][0]));//DX, DY, DZ
      fwrite(&(table_bpI_params_A3[typ_ind][0][0] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_A3[typ_ind][1][0] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_A3[typ_ind][2][0] ), binnum, 1, outfile);
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_A3[typ_ind][0][1]), &(table_bpI_params_A3[typ_ind][1][1]), &(table_bpI_params_A3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fwrite(&(table_bpI_params_A3[typ_ind][0][1] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_A3[typ_ind][1][1] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_A3[typ_ind][2][1] ), binnum, 1, outfile);
      ntot=table_bpI_N_A3[typ_ind][0]*table_bpI_N_A3[typ_ind][1]*table_bpI_N_A3[typ_ind][2];
      table_bpI_A3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	fscanf(table, "%lf", &(table_bpI_A3[typ_ind][ener]));
	fwrite(&(table_bpI_A3[typ_ind][ener] ), binnum, 1, outfile);
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
      printf("No backbone intra BP interaction for %d, glycosidic conformation ANTI.\n", i);
      table_bpI_N_A2[typ_ind][0]=-1;
      table_bpI_N_A2[typ_ind][1]=-1;
      table_bpI_N_A2[typ_ind][2]=-1;
    }
    else{
      printf("[%d ANTI P2] ", i);
      fscanf(table,"%d%d%d", &(table_bpI_N_A2[typ_ind][0]), &(table_bpI_N_A2[typ_ind][1]), &(table_bpI_N_A2[typ_ind][2])); //NX, NY, NZ
      fwrite(&(table_bpI_N_A2[typ_ind][0] ), sizeof(int), 1, outfile);
      fwrite(&(table_bpI_N_A2[typ_ind][1] ), sizeof(int), 1, outfile);
      fwrite(&(table_bpI_N_A2[typ_ind][2] ), sizeof(int), 1, outfile);
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_A2[typ_ind][0][0]), &(table_bpI_params_A2[typ_ind][1][0]), &(table_bpI_params_A2[typ_ind][2][0]));//DX, DY, DZ
      fwrite(&(table_bpI_params_A2[typ_ind][0][0] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_A2[typ_ind][1][0] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_A2[typ_ind][2][0] ), binnum, 1, outfile);
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_A2[typ_ind][0][1]), &(table_bpI_params_A2[typ_ind][1][1]), &(table_bpI_params_A2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fwrite(&(table_bpI_params_A2[typ_ind][0][1] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_A2[typ_ind][1][1] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_A2[typ_ind][2][1] ), binnum, 1, outfile);
      ntot=table_bpI_N_A2[typ_ind][0]*table_bpI_N_A2[typ_ind][1]*table_bpI_N_A2[typ_ind][2];
      table_bpI_A2[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	fscanf(table, "%lf", &(table_bpI_A2[typ_ind][ener]));
	fwrite(&(table_bpI_A2[typ_ind][ener] ), binnum, 1, outfile);
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
      printf("No backbone intra SP interaction for %d, glycosidic conformation HIGH ANTI.\n", i);
      table_bpI_N_H3[typ_ind][0]=-1;
      table_bpI_N_H3[typ_ind][1]=-1;
      table_bpI_N_H3[typ_ind][2]=-1;
      exit(ERR_INPUT);
    }
    else{
      printf("[%d HIGH ANTI P3] ", i);
      fscanf(table,"%d%d%d", &(table_bpI_N_H3[typ_ind][0]), &(table_bpI_N_H3[typ_ind][1]), &(table_bpI_N_H3[typ_ind][2])); //NX, NY, NZ
      fwrite(&(table_bpI_N_H3[typ_ind][0] ), sizeof(int), 1, outfile);
      fwrite(&(table_bpI_N_H3[typ_ind][1] ), sizeof(int), 1, outfile);
      fwrite(&(table_bpI_N_H3[typ_ind][2] ), sizeof(int), 1, outfile);
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_H3[typ_ind][0][0]), &(table_bpI_params_H3[typ_ind][1][0]), &(table_bpI_params_H3[typ_ind][2][0]));//DX, DY, DZ
      fwrite(&(table_bpI_params_H3[typ_ind][0][0] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_H3[typ_ind][1][0] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_H3[typ_ind][2][0] ), binnum, 1, outfile);
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_H3[typ_ind][0][1]), &(table_bpI_params_H3[typ_ind][1][1]), &(table_bpI_params_H3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fwrite(&(table_bpI_params_H3[typ_ind][0][1] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_H3[typ_ind][1][1] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_H3[typ_ind][2][1] ), binnum, 1, outfile);
      
      ntot=table_bpI_N_H3[typ_ind][0]*table_bpI_N_H3[typ_ind][1]*table_bpI_N_H3[typ_ind][2];
      table_bpI_H3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	fscanf(table, "%lf", &(table_bpI_H3[typ_ind][ener]));
	fwrite(&(table_bpI_H3[typ_ind][ener] ), binnum, 1, outfile);
      }
      fclose(table);
    }
  }
  table_bpI_H2=(double **)malloc(sizeof(double *)*mc_n_types);
  for(i=0;i<N_BASES;i++){
    typ_ind=i;
    sprintf(tablename, "tab_energs/table_bpI_%d_gH2.tab", i);
    if((table=fopen(tablename, "r"))==NULL){
      printf("No backbone intra SP interaction for %d, glycosidic conformation HIGH ANTI.\n", i);
      table_bpI_N_H2[typ_ind][0]=-1;
      table_bpI_N_H2[typ_ind][1]=-1;
      table_bpI_N_H2[typ_ind][2]=-1;
    }
    else{
      printf("[%d HIGH ANTI P2] ", i);
      fscanf(table,"%d%d%d", &(table_bpI_N_H2[typ_ind][0]), &(table_bpI_N_H2[typ_ind][1]), &(table_bpI_N_H2[typ_ind][2])); //NX, NY, NZ
      fwrite(&(table_bpI_N_H2[typ_ind][0] ), sizeof(int), 1, outfile);
      fwrite(&(table_bpI_N_H2[typ_ind][1] ), sizeof(int), 1, outfile);
      fwrite(&(table_bpI_N_H2[typ_ind][2] ), sizeof(int), 1, outfile);
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_H2[typ_ind][0][0]), &(table_bpI_params_H2[typ_ind][1][0]), &(table_bpI_params_H2[typ_ind][2][0]));//DX, DY, DZ
      fwrite(&(table_bpI_params_H2[typ_ind][0][0] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_H2[typ_ind][1][0] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_H2[typ_ind][2][0] ), binnum, 1, outfile);
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_H2[typ_ind][0][1]), &(table_bpI_params_H2[typ_ind][1][1]), &(table_bpI_params_H2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fwrite(&(table_bpI_params_H2[typ_ind][0][1] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_H2[typ_ind][1][1] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_H2[typ_ind][2][1] ), binnum, 1, outfile);
      ntot=table_bpI_N_H2[typ_ind][0]*table_bpI_N_H2[typ_ind][1]*table_bpI_N_H2[typ_ind][2];
      table_bpI_H2[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	fscanf(table, "%lf", &(table_bpI_H2[typ_ind][ener]));
	fwrite(&(table_bpI_H2[typ_ind][ener] ), binnum, 1, outfile);
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
      printf("No backbone intra SP interaction for %d, glycosidic conformation SYN.\n", i);
      table_bpI_N_S3[typ_ind][0]=-1;
      table_bpI_N_S3[typ_ind][1]=-1;
      table_bpI_N_S3[typ_ind][2]=-1;
    }
    else{
      printf("[%d SYN A3] ", i);
      fscanf(table,"%d%d%d", &(table_bpI_N_S3[typ_ind][0]), &(table_bpI_N_S3[typ_ind][1]), &(table_bpI_N_S3[typ_ind][2])); //NX, NY, NZ
      fwrite(&(table_bpI_N_S3[typ_ind][0] ), sizeof(int), 1, outfile);
      fwrite(&(table_bpI_N_S3[typ_ind][1] ), sizeof(int), 1, outfile);
      fwrite(&(table_bpI_N_S3[typ_ind][2] ), sizeof(int), 1, outfile);
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_S3[typ_ind][0][0]), &(table_bpI_params_S3[typ_ind][1][0]), &(table_bpI_params_S3[typ_ind][2][0]));//DX, DY, DZ
      fwrite(&(table_bpI_params_S3[typ_ind][0][0] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_S3[typ_ind][1][0] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_S3[typ_ind][2][0] ), binnum, 1, outfile);
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_S3[typ_ind][0][1]), &(table_bpI_params_S3[typ_ind][1][1]), &(table_bpI_params_S3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fwrite(&(table_bpI_params_S3[typ_ind][0][1] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_S3[typ_ind][1][1] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_S3[typ_ind][2][1] ), binnum, 1, outfile);
      ntot=table_bpI_N_S3[typ_ind][0]*table_bpI_N_S3[typ_ind][1]*table_bpI_N_S3[typ_ind][2];
      table_bpI_S3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	fscanf(table, "%lf", &(table_bpI_S3[typ_ind][ener]));
	fwrite(&(table_bpI_S3[typ_ind][ener] ), binnum, 1, outfile);
      }
      fclose(table);
    }
  }
  table_bpI_S2=(double **)malloc(sizeof(double *)*mc_n_types);
  for(i=0;i<N_BASES;i++){
    typ_ind=i;
    sprintf(tablename, "tab_energs/table_bpI_%d_gS2.tab", i);
    if((table=fopen(tablename, "r"))==NULL){
      printf("No backbone intra SP interaction for %d, glycosidic conformation SYN.\n", i);
      table_bpI_N_S2[typ_ind][0]=-1;
      table_bpI_N_S2[typ_ind][1]=-1;
      table_bpI_N_S2[typ_ind][2]=-1;
    }
    else{
      printf("[%d SYN P2] ", i);
      fscanf(table,"%d%d%d", &(table_bpI_N_S2[typ_ind][0]), &(table_bpI_N_S2[typ_ind][1]), &(table_bpI_N_S2[typ_ind][2])); //NX, NY, NZ
      fwrite(&(table_bpI_N_S2[typ_ind][0] ), sizeof(int), 1, outfile);
      fwrite(&(table_bpI_N_S2[typ_ind][1] ), sizeof(int), 1, outfile);
      fwrite(&(table_bpI_N_S2[typ_ind][2] ), sizeof(int), 1, outfile);
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_S2[typ_ind][0][0]), &(table_bpI_params_S2[typ_ind][1][0]), &(table_bpI_params_S2[typ_ind][2][0]));//DX, DY, DZ
      fwrite(&(table_bpI_params_S2[typ_ind][0][0] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_S2[typ_ind][1][0] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_S2[typ_ind][2][0] ), binnum, 1, outfile);
      fscanf(table,"%lf%lf%lf", &(table_bpI_params_S2[typ_ind][0][1]), &(table_bpI_params_S2[typ_ind][1][1]), &(table_bpI_params_S2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fwrite(&(table_bpI_params_S2[typ_ind][0][1] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_S2[typ_ind][1][1] ), binnum, 1, outfile);
      fwrite(&(table_bpI_params_S2[typ_ind][2][1] ), binnum, 1, outfile);
      ntot=table_bpI_N_S2[typ_ind][0]*table_bpI_N_S2[typ_ind][1]*table_bpI_N_S2[typ_ind][2];
      table_bpI_S2[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	fscanf(table, "%lf", &(table_bpI_S2[typ_ind][ener]));
	fwrite(&(table_bpI_S2[typ_ind][ener] ), binnum, 1, outfile);
      }
      fclose(table);
    }
  }
  
  /* SUGAR - PHOSPHATE , INTER NT */
  //GLYC ANTI
  printf("\nReading INTER-BP interactions...\n");
  table_bpB_A3=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      sprintf(tablename, "tab_energs/table_bpB_%d%d_gA3.tab", i, j);
      if((table=fopen(tablename, "r"))==NULL){
	printf("No backbone inter SP interaction for %d and %d.\n", i,j);
	table_bpB_N_A3[typ_ind][0]=-1;
	table_bpB_N_A3[typ_ind][1]=-1;
	table_bpB_N_A3[typ_ind][2]=-1;
      }
      else{
	printf("[%d %d ANTI P3] ", i,j);
	fscanf(table,"%d%d%d", &(table_bpB_N_A3[typ_ind][0]), &(table_bpB_N_A3[typ_ind][1]), &(table_bpB_N_A3[typ_ind][2])); //NX, NY, NZ
	fwrite(&(table_bpB_N_A3[typ_ind][0] ), sizeof(int), 1, outfile);
	fwrite(&(table_bpB_N_A3[typ_ind][1] ), sizeof(int), 1, outfile);
	fwrite(&(table_bpB_N_A3[typ_ind][2] ), sizeof(int), 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_A3[typ_ind][0][0]), &(table_bpB_params_A3[typ_ind][1][0]), &(table_bpB_params_A3[typ_ind][2][0]));//DX, DY, DZ
	fwrite(&(table_bpB_params_A3[typ_ind][0][0] ), binnum, 1, outfile);
	fwrite(&(table_bpB_params_A3[typ_ind][1][0] ), binnum, 1, outfile);
	fwrite(&(table_bpB_params_A3[typ_ind][2][0] ), binnum, 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_A3[typ_ind][0][1]), &(table_bpB_params_A3[typ_ind][1][1]), &(table_bpB_params_A3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	fwrite(&(table_bpB_params_A3[typ_ind][0][1] ), binnum, 1, outfile);
	fwrite(&(table_bpB_params_A3[typ_ind][1][1] ), binnum, 1, outfile);
	fwrite(&(table_bpB_params_A3[typ_ind][2][1] ), binnum, 1, outfile);
	ntot=table_bpB_N_A3[typ_ind][0]*table_bpB_N_A3[typ_ind][1]*table_bpB_N_A3[typ_ind][2];
	table_bpB_A3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_bpB_A3[typ_ind][ener]));
	  fwrite(&(table_bpB_A3[typ_ind][ener] ), binnum, 1, outfile);
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
      if((table=fopen(tablename, "r"))==NULL){
	printf("No backbone inter SP interaction for %d and %d.\n", i,j);
	table_bpB_N_A2[typ_ind][0]=-1;
	table_bpB_N_A2[typ_ind][1]=-1;
	table_bpB_N_A2[typ_ind][2]=-1;
      }
      else{
	printf("[%d %d ANTI P2] ", i,j);
	fscanf(table,"%d%d%d", &(table_bpB_N_A2[typ_ind][0]), &(table_bpB_N_A2[typ_ind][1]), &(table_bpB_N_A2[typ_ind][2])); //NX, NY, NZ
	fwrite(&(table_bpB_N_A2[typ_ind][0]) , sizeof(int), 1, outfile);
	fwrite(&(table_bpB_N_A2[typ_ind][1]) , sizeof(int), 1, outfile);
	fwrite(&(table_bpB_N_A2[typ_ind][2]) , sizeof(int), 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_A2[typ_ind][0][0]), &(table_bpB_params_A2[typ_ind][1][0]), &(table_bpB_params_A2[typ_ind][2][0]));//DX, DY, DZ
	fwrite(&(table_bpB_params_A2[typ_ind][0][0] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_A2[typ_ind][1][0] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_A2[typ_ind][2][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_A2[typ_ind][0][1]), &(table_bpB_params_A2[typ_ind][1][1]), &(table_bpB_params_A2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	fwrite(&(table_bpB_params_A2[typ_ind][0][1] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_A2[typ_ind][1][1] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_A2[typ_ind][2][1] ) , binnum, 1, outfile);
	ntot=table_bpB_N_A2[typ_ind][0]*table_bpB_N_A2[typ_ind][1]*table_bpB_N_A2[typ_ind][2];
	table_bpB_A2[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_bpB_A2[typ_ind][ener]));
	  fwrite(&(table_bpB_A2[typ_ind][ener] ) , binnum, 1, outfile);
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
	printf("No backbone inter SP interaction for %d and %d.\n", i,j);
	table_bpB_N_H3[typ_ind][0]=-1;
	table_bpB_N_H3[typ_ind][1]=-1;
	table_bpB_N_H3[typ_ind][2]=-1;
      }
      else{
	printf("[%d %d HIGH ANTI P3] ", i,j);
	fscanf(table,"%d%d%d", &(table_bpB_N_H3[typ_ind][0]), &(table_bpB_N_H3[typ_ind][1]), &(table_bpB_N_H3[typ_ind][2])); //NX, NY, NZ
	fwrite(&(table_bpB_N_H3[typ_ind][0] ) , sizeof(int), 1, outfile);
	fwrite(&(table_bpB_N_H3[typ_ind][1] ) , sizeof(int), 1, outfile);
	fwrite(&(table_bpB_N_H3[typ_ind][2] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_H3[typ_ind][0][0]), &(table_bpB_params_H3[typ_ind][1][0]), &(table_bpB_params_H3[typ_ind][2][0]));//DX, DY, DZ
	fwrite(&(table_bpB_params_H3[typ_ind][0][0] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_H3[typ_ind][1][0] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_H3[typ_ind][2][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_H3[typ_ind][0][1]), &(table_bpB_params_H3[typ_ind][1][1]), &(table_bpB_params_H3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	fwrite(&(table_bpB_params_H3[typ_ind][0][1] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_H3[typ_ind][1][1] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_H3[typ_ind][2][1] ) , binnum, 1, outfile);
	ntot=table_bpB_N_H3[typ_ind][0]*table_bpB_N_H3[typ_ind][1]*table_bpB_N_H3[typ_ind][2];
	table_bpB_H3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_bpB_H3[typ_ind][ener]));
	  fwrite(&(table_bpB_H3[typ_ind][ener] ) , binnum, 1, outfile);
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
	printf("No backbone inter SP interaction for %d and %d.\n", i,j);
	table_bpB_N_H2[typ_ind][0]=-1;
	table_bpB_N_H2[typ_ind][1]=-1;
	table_bpB_N_H2[typ_ind][2]=-1;
      }
      else{
	printf("[%d %d HIGH ANTI P2] ", i,j);
	fscanf(table,"%d%d%d", &(table_bpB_N_H2[typ_ind][0]), &(table_bpB_N_H2[typ_ind][1]), &(table_bpB_N_H2[typ_ind][2])); //NX, NY, NZ
	fwrite(&(table_bpB_N_H2[typ_ind][0] ) , sizeof(int), 1, outfile);
	fwrite(&(table_bpB_N_H2[typ_ind][1] ) , sizeof(int), 1, outfile);
	fwrite(&(table_bpB_N_H2[typ_ind][2] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_H2[typ_ind][0][0]), &(table_bpB_params_H2[typ_ind][1][0]), &(table_bpB_params_H2[typ_ind][2][0]));//DX, DY, DZ
	fwrite(&(table_bpB_params_H2[typ_ind][0][0] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_H2[typ_ind][1][0] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_H2[typ_ind][2][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_H2[typ_ind][0][1]), &(table_bpB_params_H2[typ_ind][1][1]), &(table_bpB_params_H2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	fwrite(&(table_bpB_params_H2[typ_ind][0][1] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_H2[typ_ind][1][1] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_H2[typ_ind][2][1] ) , binnum, 1, outfile);
	ntot=table_bpB_N_H2[typ_ind][0]*table_bpB_N_H2[typ_ind][1]*table_bpB_N_H2[typ_ind][2];
	table_bpB_H2[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_bpB_H2[typ_ind][ener]));
	  fwrite(&(table_bpB_H2[typ_ind][ener] ) , binnum, 1, outfile);
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
	printf("No backbone inter SP interaction for %d and %d.\n", i,j);
	table_bpB_N_S3[typ_ind][0]=-1;
	table_bpB_N_S3[typ_ind][1]=-1;
	table_bpB_N_S3[typ_ind][2]=-1;
      }
      else{
	printf("[%d %d SYN P3] ", i,j);
	fscanf(table,"%d%d%d", &(table_bpB_N_S3[typ_ind][0]), &(table_bpB_N_S3[typ_ind][1]), &(table_bpB_N_S3[typ_ind][2])); //NX, NY, NZ
	fwrite(&(table_bpB_N_S3[typ_ind][0] ) , sizeof(int), 1, outfile);
	fwrite(&(table_bpB_N_S3[typ_ind][1] ) , sizeof(int), 1, outfile);
	fwrite(&(table_bpB_N_S3[typ_ind][2] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_S3[typ_ind][0][0]), &(table_bpB_params_S3[typ_ind][1][0]), &(table_bpB_params_S3[typ_ind][2][0]));//DX, DY, DZ
	fwrite(&(table_bpB_params_S3[typ_ind][0][0] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_S3[typ_ind][1][0] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_S3[typ_ind][2][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_S3[typ_ind][0][1]), &(table_bpB_params_S3[typ_ind][1][1]), &(table_bpB_params_S3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	fwrite(&(table_bpB_params_S3[typ_ind][0][1] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_S3[typ_ind][1][1] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_S3[typ_ind][2][1] ) , binnum, 1, outfile);
	ntot=table_bpB_N_S3[typ_ind][0]*table_bpB_N_S3[typ_ind][1]*table_bpB_N_S3[typ_ind][2];
	table_bpB_S3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_bpB_S3[typ_ind][ener]));
	  fwrite(&(table_bpB_S3[typ_ind][ener] ) , binnum, 1, outfile);
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
	printf("No backbone inter SP interaction for %d and %d.\n", i,j);
	table_bpB_N_S2[typ_ind][0]=-1;
	table_bpB_N_S2[typ_ind][1]=-1;
	table_bpB_N_S2[typ_ind][2]=-1;
      }
      else{
	printf("[%d %d SYN P2] ", i,j);
	fscanf(table,"%d%d%d", &(table_bpB_N_S2[typ_ind][0]), &(table_bpB_N_S2[typ_ind][1]), &(table_bpB_N_S2[typ_ind][2])); //NX, NY, NZ
	fwrite(&(table_bpB_N_S2[typ_ind][0] ) , sizeof(int), 1, outfile);
	fwrite(&(table_bpB_N_S2[typ_ind][1] ) , sizeof(int), 1, outfile);
	fwrite(&(table_bpB_N_S2[typ_ind][2] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_S2[typ_ind][0][0]), &(table_bpB_params_S2[typ_ind][1][0]), &(table_bpB_params_S2[typ_ind][2][0]));//DX, DY, DZ
	fwrite(&(table_bpB_params_S2[typ_ind][0][0] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_S2[typ_ind][1][0] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_S2[typ_ind][2][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_bpB_params_S2[typ_ind][0][1]), &(table_bpB_params_S2[typ_ind][1][1]), &(table_bpB_params_S2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	fwrite(&(table_bpB_params_S2[typ_ind][0][1] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_S2[typ_ind][1][1] ) , binnum, 1, outfile);
	fwrite(&(table_bpB_params_S2[typ_ind][2][1] ) , binnum, 1, outfile);
	ntot=table_bpB_N_S2[typ_ind][0]*table_bpB_N_S2[typ_ind][1]*table_bpB_N_S2[typ_ind][2];
	table_bpB_S2[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_bpB_S2[typ_ind][ener]));
	  fwrite(&(table_bpB_S2[typ_ind][ener] ) , binnum, 1, outfile);
	}
	fclose(table);
      }
    }
  }
  
  printf("\nReading STACKING interactions...\n");
  /* STACKING - BONDED */
  //s35
  table_nnB_0_s35=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      sprintf(tablename, "tab_energs/table_nnB_%d%d_0_s35.tab", i,j);
      if((table=fopen(tablename, "r"))==NULL){
	printf("No stacking interaction between %d and %d.\n", i,j);
	table_nnB_N_0_s35[typ_ind][0]=-1;
	table_nnB_N_0_s35[typ_ind][1]=-1;
	table_nnB_N_0_s35[typ_ind][2]=-1;
      }
      else{
	printf("[%d %d s35] ", i,j);
	fscanf(table,"%d%d%d", &(table_nnB_N_0_s35[typ_ind][0]), &(table_nnB_N_0_s35[typ_ind][1]), &(table_nnB_N_0_s35[typ_ind][2])); //NX, NY, NZ
	fwrite(&(table_nnB_N_0_s35[typ_ind][0] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnB_N_0_s35[typ_ind][1] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnB_N_0_s35[typ_ind][2] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s35[typ_ind][0][0]), &(table_nnB_params_0_s35[typ_ind][1][0]), &(table_nnB_params_0_s35[typ_ind][2][0]));//DX, DY, DZ
	fwrite(&(table_nnB_params_0_s35[typ_ind][0][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnB_params_0_s35[typ_ind][1][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnB_params_0_s35[typ_ind][2][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s35[typ_ind][0][1]), &(table_nnB_params_0_s35[typ_ind][1][1]), &(table_nnB_params_0_s35[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	fwrite(&(table_nnB_params_0_s35[typ_ind][0][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnB_params_0_s35[typ_ind][1][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnB_params_0_s35[typ_ind][2][1] ) , binnum, 1, outfile);
	ntot=table_nnB_N_0_s35[typ_ind][0]*table_nnB_N_0_s35[typ_ind][1]*table_nnB_N_0_s35[typ_ind][2];
	table_nnB_0_s35[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnB_0_s35[typ_ind][ener]));
	  fwrite(&(table_nnB_0_s35[typ_ind][ener] ) , binnum, 1, outfile);
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
	printf("No stacking interaction between %d and %d.\n", i,j);
	table_nnB_N_0_s53[typ_ind][0]=-1;
	table_nnB_N_0_s53[typ_ind][1]=-1;
	table_nnB_N_0_s53[typ_ind][2]=-1;
      }
      else{
	printf("[%d %d s53] ", i,j);
	fscanf(table,"%d%d%d", &(table_nnB_N_0_s53[typ_ind][0]), &(table_nnB_N_0_s53[typ_ind][1]), &(table_nnB_N_0_s53[typ_ind][2])); //NX, NY, NZ
	fwrite(&(table_nnB_N_0_s53[typ_ind][0] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnB_N_0_s53[typ_ind][1] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnB_N_0_s53[typ_ind][2] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s53[typ_ind][0][0]), &(table_nnB_params_0_s53[typ_ind][1][0]), &(table_nnB_params_0_s53[typ_ind][2][0]));//DX, DY, DZ
	fwrite(&(table_nnB_params_0_s53[typ_ind][0][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnB_params_0_s53[typ_ind][1][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnB_params_0_s53[typ_ind][2][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s53[typ_ind][0][1]), &(table_nnB_params_0_s53[typ_ind][1][1]), &(table_nnB_params_0_s53[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	fwrite(&(table_nnB_params_0_s53[typ_ind][0][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnB_params_0_s53[typ_ind][1][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnB_params_0_s53[typ_ind][2][1] ) , binnum, 1, outfile);
	ntot=table_nnB_N_0_s53[typ_ind][0]*table_nnB_N_0_s53[typ_ind][1]*table_nnB_N_0_s53[typ_ind][2];
	table_nnB_0_s53[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnB_0_s53[typ_ind][ener]));
	  fwrite(&(table_nnB_0_s53[typ_ind][ener] ) , binnum, 1, outfile);
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
	printf("No stacking interaction between %d and %d.\n", i,j);
	table_nnB_N_0_s33[typ_ind][0]=-1;
	table_nnB_N_0_s33[typ_ind][1]=-1;
	table_nnB_N_0_s33[typ_ind][2]=-1;
      }
      else{
	printf("[%d %d s33] ", i,j);
	fscanf(table,"%d%d%d", &(table_nnB_N_0_s33[typ_ind][0]), &(table_nnB_N_0_s33[typ_ind][1]), &(table_nnB_N_0_s33[typ_ind][2])); //NX, NY, NZ
	fwrite(&(table_nnB_N_0_s33[typ_ind][0] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnB_N_0_s33[typ_ind][1] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnB_N_0_s33[typ_ind][2] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s33[typ_ind][0][0]), &(table_nnB_params_0_s33[typ_ind][1][0]), &(table_nnB_params_0_s33[typ_ind][2][0]));//DX, DY, DZ
	fwrite(&(table_nnB_params_0_s33[typ_ind][0][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnB_params_0_s33[typ_ind][1][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnB_params_0_s33[typ_ind][2][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s33[typ_ind][0][1]), &(table_nnB_params_0_s33[typ_ind][1][1]), &(table_nnB_params_0_s33[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	fwrite(&(table_nnB_params_0_s33[typ_ind][0][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnB_params_0_s33[typ_ind][1][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnB_params_0_s33[typ_ind][2][1] ) , binnum, 1, outfile);
	ntot=table_nnB_N_0_s33[typ_ind][0]*table_nnB_N_0_s33[typ_ind][1]*table_nnB_N_0_s33[typ_ind][2];
	table_nnB_0_s33[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnB_0_s33[typ_ind][ener]));
	  fwrite(&(table_nnB_0_s33[typ_ind][ener] ) , binnum, 1, outfile);
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
	printf("No stacking interaction between %d and %d.\n", i,j);
	table_nnB_N_0_s55[typ_ind][0]=-1;
	table_nnB_N_0_s55[typ_ind][1]=-1;
	table_nnB_N_0_s55[typ_ind][2]=-1;
      }
      else{
	printf("[%d %d s55] ", i,j);
	fscanf(table,"%d%d%d", &(table_nnB_N_0_s55[typ_ind][0]), &(table_nnB_N_0_s55[typ_ind][1]), &(table_nnB_N_0_s55[typ_ind][2])); //NX, NY, NZ
	fwrite(&(table_nnB_N_0_s55[typ_ind][0] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnB_N_0_s55[typ_ind][1] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnB_N_0_s55[typ_ind][2] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s55[typ_ind][0][0]), &(table_nnB_params_0_s55[typ_ind][1][0]), &(table_nnB_params_0_s55[typ_ind][2][0]));//DX, DY, DZ
	fwrite(&(table_nnB_params_0_s55[typ_ind][0][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnB_params_0_s55[typ_ind][1][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnB_params_0_s55[typ_ind][2][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s55[typ_ind][0][1]), &(table_nnB_params_0_s55[typ_ind][1][1]), &(table_nnB_params_0_s55[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	fwrite(&(table_nnB_params_0_s55[typ_ind][0][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnB_params_0_s55[typ_ind][1][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnB_params_0_s55[typ_ind][2][1] ) , binnum, 1, outfile);
	ntot=table_nnB_N_0_s55[typ_ind][0]*table_nnB_N_0_s55[typ_ind][1]*table_nnB_N_0_s55[typ_ind][2];
	table_nnB_0_s55[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnB_0_s55[typ_ind][ener]));
	  fwrite(&(table_nnB_0_s55[typ_ind][ener] ) , binnum, 1, outfile);
	}
	fclose(table);
      }
    }
  }
  //////////////////////////////////
  printf("\nReading STACKING-DIHEDRAL interactions...\n");
  table_nnB_1_s33=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //typ_ind2=N_BASES*j+i;
      sprintf(tablename, "tab_energs/table_nnB_%d%d_1_s33.tab", i,j);
      //FILE *table;
      if((table=fopen(tablename, "r"))==NULL){
	printf("[NO %d %d 33] ", i,j);
	table_nnB_N_1_s33[typ_ind]=-1;
	//table_nnB_N_1[typ_ind2]=-1;
	//exit(ERR_INPUT);
      }
      else{
	printf("[%d %d 33] ", i,j);
	fscanf(table,"%d", &(table_nnB_N_1_s33[typ_ind])); //N_ETA
	fwrite(&(table_nnB_N_1_s33[typ_ind] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf", &(table_nnB_params_1_s33[typ_ind][0]));//D_ETA
	fwrite(&(table_nnB_params_1_s33[typ_ind][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf", &(table_nnB_params_1_s33[typ_ind][1])); //ETAMIN
	fwrite(&(table_nnB_params_1_s33[typ_ind][1] ) , binnum, 1, outfile);
	ntot=table_nnB_N_1_s33[typ_ind];
	table_nnB_1_s33[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnB_1_s33[typ_ind][ener]));
	  fwrite(&(table_nnB_1_s33[typ_ind][ener] ) , binnum, 1, outfile);
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
	printf("[NO %d %d 35] ", i,j);
	table_nnB_N_1_s35[typ_ind]=-1;
	//table_nnB_N_1[typ_ind2]=-1;
	//exit(ERR_INPUT);
      }
      else{
	printf("[%d %d 35] ", i,j);
	fscanf(table,"%d", &(table_nnB_N_1_s35[typ_ind])); //N_ETA
	fwrite(&(table_nnB_N_1_s35[typ_ind] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf", &(table_nnB_params_1_s35[typ_ind][0]));//D_ETA
	fwrite(&(table_nnB_params_1_s35[typ_ind][0] ) , binnum, 1,  outfile);
	fscanf(table,"%lf", &(table_nnB_params_1_s35[typ_ind][1])); //ETAMIN
	fwrite(&(table_nnB_params_1_s35[typ_ind][1] ) , binnum, 1, outfile);
	ntot=table_nnB_N_1_s35[typ_ind];
	table_nnB_1_s35[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnB_1_s35[typ_ind][ener]));
	  fwrite(&(table_nnB_1_s35[typ_ind][ener] ) , binnum, 1, outfile);
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
	printf("[NO %d %d 53] ", i,j);
	table_nnB_N_1_s53[typ_ind]=-1;
	//table_nnB_N_1[typ_ind2]=-1;
	//exit(ERR_INPUT);
      }
      else{
	printf("[%d %d 53] ", i,j);
	fscanf(table,"%d", &(table_nnB_N_1_s53[typ_ind])); //N_ETA
	fwrite(&(table_nnB_N_1_s53[typ_ind] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf", &(table_nnB_params_1_s53[typ_ind][0]));//D_ETA
	fwrite(&(table_nnB_params_1_s53[typ_ind][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf", &(table_nnB_params_1_s53[typ_ind][1])); //ETAMIN
	fwrite(&(table_nnB_params_1_s53[typ_ind][1] ) , binnum, 1, outfile);
	ntot=table_nnB_N_1_s53[typ_ind];
	table_nnB_1_s53[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnB_1_s53[typ_ind][ener]));
	  fwrite(&(table_nnB_1_s53[typ_ind][ener] ) , binnum, 1, outfile);
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
	printf("[NO %d %d 55] ", i,j);
	table_nnB_N_1_s55[typ_ind]=-1;
	//table_nnB_N_1[typ_ind2]=-1;
	//exit(ERR_INPUT);
      }
      else{
	printf("[%d %d 55] ", i,j);
	fscanf(table,"%d", &(table_nnB_N_1_s55[typ_ind])); //N_ETA
	fwrite(&(table_nnB_N_1_s55[typ_ind] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf", &(table_nnB_params_1_s55[typ_ind][0]));//D_ETA
	fwrite(&(table_nnB_params_1_s55[typ_ind][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf", &(table_nnB_params_1_s55[typ_ind][1])); //ETAMIN
	fwrite(&(table_nnB_params_1_s55[typ_ind][1] ) , binnum, 1, outfile);
	ntot=table_nnB_N_1_s55[typ_ind];
	table_nnB_1_s55[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnB_1_s55[typ_ind][ener]));
	  fwrite(&(table_nnB_1_s55[typ_ind][ener] ) , binnum, 1, outfile);
	}
	fclose(table);
      }
    }
  }
  
  printf("\nReading NON-BONDED STACKING interactions...\n");
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
	printf("No generic non-bonded stacking interaction.\n");
	table_nnN_N_0s3[typ_ind][0]=-1;
	table_nnN_N_0s3[typ_ind][1]=-1;
	table_nnN_N_0s3[typ_ind][2]=-1;
      }
      else{
	printf("[%d %d s3] ", i,j);
	fscanf(table,"%d%d%d", &(table_nnN_N_0s3[typ_ind][0]), &(table_nnN_N_0s3[typ_ind][1]), &(table_nnN_N_0s3[typ_ind][2])); //NX, NY, NZ
	fwrite(&(table_nnN_N_0s3[typ_ind][0] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnN_N_0s3[typ_ind][1] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnN_N_0s3[typ_ind][2] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_0s3[typ_ind][0][0]), &(table_nnN_params_0s3[typ_ind][1][0]), &(table_nnN_params_0s3[typ_ind][2][0]));//DX, DY, DZ
	fwrite(&(table_nnN_params_0s3[typ_ind][0][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_0s3[typ_ind][1][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_0s3[typ_ind][2][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_0s3[typ_ind][0][1]), &(table_nnN_params_0s3[typ_ind][1][1]), &(table_nnN_params_0s3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	fwrite(&(table_nnN_params_0s3[typ_ind][0][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_0s3[typ_ind][1][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_0s3[typ_ind][2][1] ) , binnum, 1, outfile);
	ntot=table_nnN_N_0s3[typ_ind][0]*table_nnN_N_0s3[typ_ind][1]*table_nnN_N_0s3[typ_ind][2];
	table_nnN_0s3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnN_0s3[typ_ind][ener]));
	  fwrite(&(table_nnN_0s3[typ_ind][ener] ) , binnum, 1, outfile);
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
	printf("No generic non-bonded stacking_inv interaction.\n");
	table_nnN_N_0s5[typ_ind][0]=-1;
	table_nnN_N_0s5[typ_ind][1]=-1;
	table_nnN_N_0s5[typ_ind][2]=-1;
      }
      else{
	printf("[%d %d s5] ", i,j);
	fscanf(table,"%d%d%d", &(table_nnN_N_0s5[typ_ind][0]), &(table_nnN_N_0s5[typ_ind][1]), &(table_nnN_N_0s5[typ_ind][2])); //NX, NY, NZ
	fwrite(&(table_nnN_N_0s5[typ_ind][0] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnN_N_0s5[typ_ind][1] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnN_N_0s5[typ_ind][2] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_0s5[typ_ind][0][0]), &(table_nnN_params_0s5[typ_ind][1][0]), &(table_nnN_params_0s5[typ_ind][2][0]));//DX, DY, DZ
	fwrite(&(table_nnN_params_0s5[typ_ind][0][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_0s5[typ_ind][1][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_0s5[typ_ind][2][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_0s5[typ_ind][0][1]), &(table_nnN_params_0s5[typ_ind][1][1]), &(table_nnN_params_0s5[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	fwrite(&(table_nnN_params_0s5[typ_ind][0][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_0s5[typ_ind][1][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_0s5[typ_ind][2][1] ) , binnum, 1, outfile);
	ntot=table_nnN_N_0s5[typ_ind][0]*table_nnN_N_0s5[typ_ind][1]*table_nnN_N_0s5[typ_ind][2];
	table_nnN_0s5[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnN_0s5[typ_ind][ener]));
	  fwrite(&(table_nnN_0s5[typ_ind][ener] ) , binnum, 1, outfile);
	}
	fclose(table);
      }
    }
  }
  
  /* WATSON-CRICK */
  printf("\nReading BASE-PAIR interactions...\n");
  table_nnN_2=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //typ_ind2=N_BASES*j+i;
      sprintf(tablename, "tab_energs/table_nnN_%d%d_2.tab", i,j);
      //FILE *table;
      if((table=fopen(tablename, "r"))==NULL){
	printf("[NO %d %d] ", i,j);
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
	printf("[%d %d] ", i,j);
	fscanf(table,"%d%d%d", &(table_nnN_N_2[typ_ind][0]), &(table_nnN_N_2[typ_ind][1]), &(table_nnN_N_2[typ_ind][2])); //NX, NY, NZ
	fwrite(&(table_nnN_N_2[typ_ind][0] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnN_N_2[typ_ind][1] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnN_N_2[typ_ind][2] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_2[typ_ind][0][0]), &(table_nnN_params_2[typ_ind][1][0]), &(table_nnN_params_2[typ_ind][2][0]));//DX, DY, DZ
	fwrite(&(table_nnN_params_2[typ_ind][0][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_2[typ_ind][1][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_2[typ_ind][2][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_2[typ_ind][0][1]), &(table_nnN_params_2[typ_ind][1][1]), &(table_nnN_params_2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	fwrite(&(table_nnN_params_2[typ_ind][0][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_2[typ_ind][1][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_2[typ_ind][2][1] ) , binnum, 1, outfile);
	ntot=table_nnN_N_2[typ_ind][0]*table_nnN_N_2[typ_ind][1]*table_nnN_N_2[typ_ind][2];
	table_nnN_N_2[typ_ind][3]=ntot;
	table_nnN_2[typ_ind]=(double *)malloc(sizeof(double)*ntot*WC_FACES);
	for(ener=0;ener<ntot*WC_FACES;ener++){
	  fscanf(table, "%lf", &(table_nnN_2[typ_ind][ener]));
	  fwrite(&(table_nnN_2[typ_ind][ener] ) , binnum, 1, outfile);
	  //table_nnN_2[typ_ind2][ener]=table_nnN_2[typ_ind][ener];
	}
	fclose(table);
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
	printf("[NO inv %d %d] ", i,j);
	table_nnN_N_2_inv[typ_ind][0]=-1;
	table_nnN_N_2_inv[typ_ind][1]=-1;
	table_nnN_N_2_inv[typ_ind][2]=-1;
	table_nnN_N_2_inv[typ_ind][3]=-1;
	//table_nnN_N_2_inv[typ_ind2][0]=-1;
	//table_nnN_N_2_inv[typ_ind2][1]=-1;
	//table_nnN_N_2_inv[typ_ind2][2]=-1;
      }
      else{
	printf("[inv %d %d] ", i,j);
	fscanf(table,"%d%d%d", &(table_nnN_N_2_inv[typ_ind][0]), &(table_nnN_N_2_inv[typ_ind][1]), &(table_nnN_N_2_inv[typ_ind][2])); //NX, NY, NZ
	fwrite(&(table_nnN_N_2_inv[typ_ind][0] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnN_N_2_inv[typ_ind][1] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnN_N_2_inv[typ_ind][2] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_2_inv[typ_ind][0][0]), &(table_nnN_params_2_inv[typ_ind][1][0]), &(table_nnN_params_2_inv[typ_ind][2][0]));//DX, DY, DZ
	fwrite(&(table_nnN_params_2_inv[typ_ind][0][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_2_inv[typ_ind][1][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_2_inv[typ_ind][2][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_2_inv[typ_ind][0][1]), &(table_nnN_params_2_inv[typ_ind][1][1]), &(table_nnN_params_2_inv[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	fwrite(&(table_nnN_params_2_inv[typ_ind][0][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_2_inv[typ_ind][1][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_2_inv[typ_ind][2][1] ) , binnum, 1, outfile);
	ntot=table_nnN_N_2_inv[typ_ind][0]*table_nnN_N_2_inv[typ_ind][1]*table_nnN_N_2_inv[typ_ind][2];
	table_nnN_N_2_inv[typ_ind][3]=ntot;
	table_nnN_2_inv[typ_ind]=(double *)malloc(sizeof(double)*ntot*WC_FACES);
	for(ener=0;ener<ntot*WC_FACES;ener++){
	  fscanf(table, "%lf", &(table_nnN_2_inv[typ_ind][ener]));
	  //table_nnN_2_inv[typ_ind2][ener]=table_nnN_2_inv[typ_ind][ener];
	  fwrite(&(table_nnN_2_inv[typ_ind][ener] ) , binnum, 1, outfile);
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
	printf("[NO DIH %d %d] ", i,j);
	table_nnN_N_3[typ_ind]=-1;
	//table_nnN_N_3[typ_ind2]=-1;
	exit(ERR_INPUT);
      }
      else{
	printf("[DIH %d %d] ", i,j);
	fscanf(table,"%d", &(table_nnN_N_3[typ_ind])); //N_THETA
	fwrite(&(table_nnN_N_3[typ_ind] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf", &(table_nnN_params_3[typ_ind][0]));//D_THETA
	fwrite(&(table_nnN_params_3[typ_ind][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf", &(table_nnN_params_3[typ_ind][1])); //THETAMIN
	fwrite(&(table_nnN_params_3[typ_ind][1] ) , binnum, 1, outfile);
	ntot=table_nnN_N_3[typ_ind];
	table_nnN_3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	//table_nnN_3[typ_ind2]=(double *)malloc(sizeof(double)*ntot);
	
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnN_3[typ_ind][ener]));
	  fwrite(&(table_nnN_3[typ_ind][ener] ) , binnum, 1, outfile);
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
	printf("[NO F %d %d] ", i,j);
	table_nnN_N_2_F[typ_ind][0]=-1;
	table_nnN_N_2_F[typ_ind][1]=-1;
	table_nnN_N_2_F[typ_ind][2]=-1;
	table_nnN_N_2_F[typ_ind][3]=-1;
      }
      else{
	printf("[F %d %d] ", i,j);
	fscanf(table,"%d%d%d", &(table_nnN_N_2_F[typ_ind][0]), &(table_nnN_N_2_F[typ_ind][1]), &(table_nnN_N_2_F[typ_ind][2])); //NX, NY, NZ
	fwrite(&(table_nnN_N_2_F[typ_ind][0] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnN_N_2_F[typ_ind][1] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnN_N_2_F[typ_ind][2] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_2_F[typ_ind][0][0]), &(table_nnN_params_2_F[typ_ind][1][0]), &(table_nnN_params_2_F[typ_ind][2][0]));//DX, DY, DZ
	fwrite(&(table_nnN_params_2_F[typ_ind][0][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_2_F[typ_ind][1][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_2_F[typ_ind][2][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_2_F[typ_ind][0][1]), &(table_nnN_params_2_F[typ_ind][1][1]), &(table_nnN_params_2_F[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	fwrite(&(table_nnN_params_2_F[typ_ind][0][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_2_F[typ_ind][1][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_2_F[typ_ind][2][1] ) , binnum, 1, outfile);
	ntot=table_nnN_N_2_F[typ_ind][0]*table_nnN_N_2_F[typ_ind][1]*table_nnN_N_2_F[typ_ind][2];
	table_nnN_N_2_F[typ_ind][3]=ntot;
	table_nnN_2_F[typ_ind]=(double *)malloc(sizeof(double)*ntot*WC_FACES);
	for(ener=0;ener<ntot*WC_FACES;ener++){
	  fscanf(table, "%lf", &(table_nnN_2_F[typ_ind][ener]));
	  fwrite(&(table_nnN_2_F[typ_ind][ener] ) , binnum, 1, outfile);
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
	printf("[NO F inv %d %d] ", i,j);
	table_nnN_N_2_inv_F[typ_ind][0]=-1;
	table_nnN_N_2_inv_F[typ_ind][1]=-1;
	table_nnN_N_2_inv_F[typ_ind][2]=-1;
	table_nnN_N_2_inv_F[typ_ind][3]=-1;
      }
      else{
	printf("[F inv %d %d] ", i,j);
	fscanf(table,"%d%d%d", &(table_nnN_N_2_inv_F[typ_ind][0]), &(table_nnN_N_2_inv_F[typ_ind][1]), &(table_nnN_N_2_inv_F[typ_ind][2])); //NX, NY, NZ
	fwrite(&(table_nnN_N_2_inv_F[typ_ind][0] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnN_N_2_inv_F[typ_ind][1] ) , sizeof(int), 1, outfile);
	fwrite(&(table_nnN_N_2_inv_F[typ_ind][2] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_2_inv_F[typ_ind][0][0]), &(table_nnN_params_2_inv_F[typ_ind][1][0]), &(table_nnN_params_2_inv_F[typ_ind][2][0]));//DX, DY, DZ
	fwrite(&(table_nnN_params_2_inv_F[typ_ind][0][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_2_inv_F[typ_ind][1][0] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_2_inv_F[typ_ind][2][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf%lf%lf", &(table_nnN_params_2_inv_F[typ_ind][0][1]), &(table_nnN_params_2_inv_F[typ_ind][1][1]), &(table_nnN_params_2_inv_F[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	fwrite(&(table_nnN_params_2_inv_F[typ_ind][0][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_2_inv_F[typ_ind][1][1] ) , binnum, 1, outfile);
	fwrite(&(table_nnN_params_2_inv_F[typ_ind][2][1] ) , binnum, 1, outfile);
	ntot=table_nnN_N_2_inv_F[typ_ind][0]*table_nnN_N_2_inv_F[typ_ind][1]*table_nnN_N_2_inv_F[typ_ind][2];
	table_nnN_N_2_inv_F[typ_ind][3]=ntot;
	table_nnN_2_inv_F[typ_ind]=(double *)malloc(sizeof(double)*ntot*WC_FACES);
	for(ener=0;ener<ntot*WC_FACES;ener++){
	  fscanf(table, "%lf", &(table_nnN_2_inv_F[typ_ind][ener]));
	  fwrite(&(table_nnN_2_inv_F[typ_ind][ener] ) , binnum, 1, outfile);
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
	printf("[NO F DIH %d %d] ", i,j);
	table_nnN_N_3_F[typ_ind]=-1;
	exit(ERR_INPUT);
      }
      else{
	printf("[F DIH %d %d] ", i,j);
	fscanf(table,"%d", &(table_nnN_N_3_F[typ_ind])); //N_THETA
	fwrite(&(table_nnN_N_3_F[typ_ind] ) , sizeof(int), 1, outfile);
	fscanf(table,"%lf", &(table_nnN_params_3_F[typ_ind][0]));//D_THETA
	fwrite(&(table_nnN_params_3_F[typ_ind][0] ) , binnum, 1, outfile);
	fscanf(table,"%lf", &(table_nnN_params_3_F[typ_ind][1])); //THETAMIN
	fwrite(&(table_nnN_params_3_F[typ_ind][1] ) , binnum, 1, outfile);
	ntot=table_nnN_N_3_F[typ_ind];
	table_nnN_3_F[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  fscanf(table, "%lf", &(table_nnN_3_F[typ_ind][ener]));
	  fwrite(&(table_nnN_3_F[typ_ind][ener] ) , binnum, 1, outfile);
	}
	fclose(table);
      }
    }
  }
  /***********************/
  
  printf("\nReading BASE-PHOSPHATE interactions...\n");
  // BASE-PHOSPHATE
  table_npN_0=(double **)malloc(sizeof(double *)*N_BASES);
  for(i=0;i<N_BASES;i++){
    sprintf(tablename, "tab_energs/table_npN_%d_0.tab", i);
    if((table=fopen(tablename, "r"))==NULL){
      printf("[NO BPH %d] ", i);
      table_npN_N_0[i][0]=-1;
      table_npN_N_0[i][1]=-1;
      table_npN_N_0[i][2]=-1;
    }
    else{
      printf("[BPH %d] ", i);
      fscanf(table,"%d%d%d", &(table_npN_N_0[i][0]), &(table_npN_N_0[i][1]), &(table_npN_N_0[i][2])); //NX, NY, NZ
      fwrite(&(table_npN_N_0[i][0] ) , sizeof(int), 1, outfile);
      fwrite(&(table_npN_N_0[i][1] ) , sizeof(int), 1, outfile);
      fwrite(&(table_npN_N_0[i][2] ) , sizeof(int), 1, outfile);
      fscanf(table,"%lf%lf%lf", &(table_npN_params_0[i][0][0]), &(table_npN_params_0[i][1][0]), &(table_npN_params_0[i][2][0]));//DX, DY, DZ
      fwrite(&(table_npN_params_0[i][0][0] ) , binnum, 1, outfile);
      fwrite(&(table_npN_params_0[i][1][0] ) , binnum, 1, outfile);
      fwrite(&(table_npN_params_0[i][2][0] ) , binnum, 1, outfile);
      fscanf(table,"%lf%lf%lf", &(table_npN_params_0[i][0][1]), &(table_npN_params_0[i][1][1]), &(table_npN_params_0[i][2][1])); //XMIN, YMIN, ZMIN
      fwrite(&(table_npN_params_0[i][0][1] ) , binnum, 1, outfile);
      fwrite(&(table_npN_params_0[i][1][1] ) , binnum, 1, outfile);
      fwrite(&(table_npN_params_0[i][2][1] ) , binnum, 1, outfile);
      ntot=table_npN_N_0[i][0]*table_npN_N_0[i][1]*table_npN_N_0[i][2];
      table_npN_0[i]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	fscanf(table, "%lf", &(table_npN_0[i][ener]));
	fwrite(&(table_npN_0[i][ener] ) , binnum, 1, outfile);
      }
      fclose(table);
    }
  }
  printf("\n");
  fclose(outfile);
}





void MC_read_bin_energy_tables(){
  mc_n_types=N_BASES;
  size_t binnum=sizeof(float);
  FILE *outfile;
  char tablename[256];
  char outname[256];
  int i,j, typ_ind, ener, ntot, stf;
  int ntypsq=mc_n_types*mc_n_types;
  char ctemp1, ctemp2, ctemp3, ctemp4, ctemp5;
  sprintf(outname, "intrac.btb");
  if((outfile=fopen(outname, "rb"))==NULL){
    printf("ERROR: Could not open %s file for reading tables.\n", outname);
    exit(ERR_INPUT);
  }
  
  /* NON BONDED POTENTIAL WELLS */
  //sprintf(tablename, "tab_energs/table_wells.tab");
  //FILE *table;
  //if((table=fopen(tablename,"r"))==NULL){
  //  printf("file not found!\n");
  //  exit(1);
  // }
  //else{
  //fscanf(table,"%c%c", &ctemp1, &ctemp2);
  //if(ctemp1!='b' || ctemp2 != 'b'){
  //  printf("Wrong syntax at table_wells.tab file (bb)!\n%c %c\n", ctemp1, ctemp2);
  //  exit(1);}
  //fscanf(table, "%lf", &(BB_PREF));
  fread(&BB_PREF, binnum, 1, outfile); 
  //fscanf(table, "%lf", &(BB_PREF_A));
  fread(&BB_PREF_A, binnum, 1, outfile);
  //fscanf(table,"%c%c%c", &ctemp1, &ctemp2, &ctemp3);
  //if(ctemp2!='g' || ctemp3 != 'l'){
  //   printf("Wrong syntax at table_wells.tab file (gl)!\n%c %c\n", ctemp1, ctemp2);
  //   exit(1);
  // }
  for(j=0;j<N_PUCK_STATES;j++)
    for(i=0;i<N_GLYC_STATES;i++){
      //fscanf(table, "%lf", &(glp_well_R[i][j]));
      fread(&(glp_well_R[i][j]), binnum, 1, outfile);
    }
  for(j=0;j<N_PUCK_STATES;j++)
    for(i=0;i<N_GLYC_STATES;i++){
      //fscanf(table, "%lf", &(glp_well_Y[i][j]));
      fread(&(glp_well_Y[i][j]), binnum, 1, outfile);
    }
  
  /* NON BONDED POTENTIAL WELLS */
  //fscanf(table,"%c%c%c%c", &ctemp1, &ctemp2, &ctemp3, &ctemp4);
  //if(ctemp2!='s' || ctemp3 != 't' || ctemp4!='b'){
  //  printf("Wrong syntax at table_wells.tab file (st , 1)!\n%c %c %c %c\n", ctemp1, ctemp2, ctemp3, ctemp4);
  //  exit(1);}
  for(stf=0;stf<N_STFACES;stf++){
    for(i=0;i<N_BASES;i++)
      for(j=0;j<N_BASES;j++){
	//fscanf(table, "%lf", &(b_st_well[N_BASES*i+j][stf]));
	fread(&(b_st_well[N_BASES*i+j][stf]), binnum, 1, outfile);
      }
  }
  //fscanf(table,"%c%c%c",&ctemp3, &ctemp1, &ctemp2);
  //if(ctemp1!='s' || ctemp2 != 't'){
  //  printf("Wrong syntax at table_wells.tab file (st , 2)!\n%c %c %c", ctemp1, ctemp2, ctemp3);
  //  exit(1);}
  for(i=0;i<N_BASES;i++)
    for(j=0;j<N_BASES;j++){
      //fscanf(table, "%lf", &(nb_st_well[N_BASES*i+j]));
      fread(&(nb_st_well[N_BASES*i+j]), binnum, 1, outfile);
    }
  for(i=0;i<N_BASES;i++)
    for(j=0;j<N_BASES;j++){
      if(nb_st_well[N_BASES*i+j]!=nb_st_well[N_BASES*j+i]){
	printf("Non-bonded matrix (stacking) is not symmetric!\n");
	exit(1);
      }
    }
  //fscanf(table,"%c%c%c", &ctemp3,&ctemp1, &ctemp2);
  //if(ctemp1!='w' || ctemp2 != 'c'){
  //printf("Wrong syntax at table_wells.tab file (wc)!\n %c %c %c", ctemp1, ctemp2, ctemp3);
  //exit(1);}
  for(i=0;i<N_BASES;i++)
    for(j=0;j<N_BASES;j++){
      //fscanf(table,"%lf", &(nb_wc_well[N_BASES*i+j][0*WC_FACES+0]));
      fread(&(nb_wc_well[N_BASES*i+j][0*WC_FACES+0]), binnum, 1, outfile);
      //fscanf(table,"%lf", &(nb_wc_well[N_BASES*i+j][1*WC_FACES+1]));
      fread(&(nb_wc_well[N_BASES*i+j][1*WC_FACES+1]), binnum, 1, outfile);
      //fscanf(table,"%lf", &(nb_wc_well[N_BASES*i+j][2*WC_FACES+2]));
      fread(&(nb_wc_well[N_BASES*i+j][2*WC_FACES+2]), binnum, 1, outfile);
      
      //fscanf(table,"%lf", &(nb_wc_well[N_BASES*i+j][1*WC_FACES+2]));
      fread(&(nb_wc_well[N_BASES*i+j][1*WC_FACES+2]), binnum, 1, outfile);
      nb_wc_well[N_BASES*j+i][2*WC_FACES+1]=nb_wc_well[N_BASES*i+j][1*WC_FACES+2];
      //fscanf(table,"%lf", &(nb_wc_well[N_BASES*i+j][1*WC_FACES+0]));
      fread(&(nb_wc_well[N_BASES*i+j][1*WC_FACES+0]), binnum, 1, outfile);
      nb_wc_well[N_BASES*j+i][0*WC_FACES+1]=nb_wc_well[N_BASES*i+j][1*WC_FACES+0];
      //fscanf(table,"%lf", &(nb_wc_well[N_BASES*i+j][2*WC_FACES+0]));
      fread(&(nb_wc_well[N_BASES*i+j][2*WC_FACES+0]), binnum, 1, outfile);
      nb_wc_well[N_BASES*j+i][0*WC_FACES+2]=nb_wc_well[N_BASES*i+j][2*WC_FACES+0];
    }
    
  //BASE-PAIRING IN PARALLEL CONFORMATION
  //fscanf(table,"%c%c%c", &ctemp3,&ctemp1, &ctemp2);
  //if(ctemp1!='w' || ctemp2 != 'P' ){
  //  printf("Wrong syntax at table_wells.tab file (wc)!\n %c %c %c", ctemp1, ctemp2, ctemp3);
  //  exit(1);}
  for(i=0;i<N_BASES;i++)
    for(j=0;j<N_BASES;j++){
      //fscanf(table,"%lf", &(nb_wc_well_F[N_BASES*i+j][0*WC_FACES+0]));
      fread(&(nb_wc_well_F[N_BASES*i+j][0*WC_FACES+0]), binnum, 1, outfile);
      //fscanf(table,"%lf", &(nb_wc_well_F[N_BASES*i+j][1*WC_FACES+1]));
      fread(&(nb_wc_well_F[N_BASES*i+j][1*WC_FACES+1]), binnum, 1, outfile);
      //fscanf(table,"%lf", &(nb_wc_well_F[N_BASES*i+j][2*WC_FACES+2]));
      fread(&(nb_wc_well_F[N_BASES*i+j][2*WC_FACES+2]), binnum, 1, outfile);
      
      //fscanf(table,"%lf", &(nb_wc_well_F[N_BASES*i+j][1*WC_FACES+2]));
      fread(&(nb_wc_well_F[N_BASES*i+j][1*WC_FACES+2]), binnum, 1, outfile);
      nb_wc_well_F[N_BASES*j+i][2*WC_FACES+1]=nb_wc_well_F[N_BASES*i+j][1*WC_FACES+2];
      //fscanf(table,"%lf", &(nb_wc_well_F[N_BASES*i+j][1*WC_FACES+0]));
      fread(&(nb_wc_well_F[N_BASES*i+j][1*WC_FACES+0]), binnum, 1, outfile);
      nb_wc_well_F[N_BASES*j+i][0*WC_FACES+1]=nb_wc_well_F[N_BASES*i+j][1*WC_FACES+0];
      //fscanf(table,"%lf", &(nb_wc_well_F[N_BASES*i+j][2*WC_FACES+0]));
      fread(&(nb_wc_well_F[N_BASES*i+j][2*WC_FACES+0]), binnum, 1, outfile);
      nb_wc_well_F[N_BASES*j+i][0*WC_FACES+2]=nb_wc_well_F[N_BASES*i+j][2*WC_FACES+0];
    }
  
  /************** BASE PHOSPHATE ***************/
  //fscanf(table,"%c%c%c", &ctemp3,&ctemp1, &ctemp2);
  //printf("reading wc %c  %c\n", ctemp1, ctemp2);
  //if(ctemp1!='b' || ctemp2 != 'p'){
  //printf("Wrong syntax at table_wells.tab file (wc)!\n %c %c %c", ctemp1, ctemp2, ctemp3);
  //exit(1);}
  for(i=0;i<N_BASES;i++)
    for(j=0;j<WC_FACES;j++){
      //fscanf(table,"%lf", &(nb_bp_well[i][j]));
      fread(&(nb_bp_well[i][j]), binnum, 1, outfile);
    }
  //fscanf(table, "%lf", &(nb_bp_spec_well[2][0]));
  fread(&(nb_bp_spec_well[2][0]), binnum, 1, outfile);
  //fscanf(table, "%lf", &(nb_bp_spec_well[2][1]));
  fread(&(nb_bp_spec_well[2][1]), binnum, 1, outfile);
  //fscanf(table, "%lf", &(nb_bp_spec_well[2][2]));
  fread(&(nb_bp_spec_well[2][2]), binnum, 1, outfile);
  //fclose(table);
  //}
  
  /************************ BASE-PAIR SECOND DIHEDRALS **************************/
  
  //sprintf(tablename, "tab_energs/table_wc_secdih.tab");
  //if((table=fopen(tablename,"r"))==NULL){
  //  printf("file not found!\n");
  //printf("No potential wells for non bonded interactions loaded.\n");
  // exit(1);
  //}
  //else{
  //fscanf(table,"%c%c%c%c", &ctemp1,&ctemp2, &ctemp3, &ctemp4);
  //printf("reading wc %c  %c\n", ctemp1, ctemp2);
  //if(ctemp1!='A' || ctemp2 != 'm' || ctemp3 != 'i' || ctemp4!='n'){
  //printf("Wrong syntax at table_wc_secdih.tab file (wc)!\n %c %c %c %c", ctemp1, ctemp2, ctemp3, ctemp4);
  //  exit(1);}
  for(i=0;i<N_BASES;i++)
    for(j=0;j<N_BASES;j++){
      //fscanf(table,"%lf", &(wc_secdih_min[N_BASES*i+j][0*WC_FACES+0]));
      fread(&(wc_secdih_min[N_BASES*i+j][0*WC_FACES+0]), binnum, 1, outfile);
      //fscanf(table,"%lf", &(wc_secdih_min[N_BASES*i+j][1*WC_FACES+1]));
      fread(&(wc_secdih_min[N_BASES*i+j][1*WC_FACES+1]), binnum, 1, outfile);
      //fscanf(table,"%lf", &(wc_secdih_min[N_BASES*i+j][2*WC_FACES+2]));
      fread(&(wc_secdih_min[N_BASES*i+j][2*WC_FACES+2]), binnum, 1, outfile);
      
      //fscanf(table,"%lf", &(wc_secdih_min[N_BASES*i+j][1*WC_FACES+2]));
      fread(&(wc_secdih_min[N_BASES*i+j][1*WC_FACES+2]), binnum, 1, outfile);
      wc_secdih_min[N_BASES*j+i][2*WC_FACES+1]=wc_secdih_min[N_BASES*i+j][1*WC_FACES+2];
      //fscanf(table,"%lf", &(wc_secdih_min[N_BASES*i+j][1*WC_FACES+0]));
      fread(&(wc_secdih_min[N_BASES*i+j][1*WC_FACES+0]), binnum, 1, outfile);
      wc_secdih_min[N_BASES*j+i][0*WC_FACES+1]=wc_secdih_min[N_BASES*i+j][1*WC_FACES+0];
      //fscanf(table,"%lf", &(wc_secdih_min[N_BASES*i+j][2*WC_FACES+0]));
      fread(&(wc_secdih_min[N_BASES*i+j][2*WC_FACES+0]), binnum, 1, outfile);
      wc_secdih_min[N_BASES*j+i][0*WC_FACES+2]=wc_secdih_min[N_BASES*i+j][2*WC_FACES+0];
    }
  //fscanf(table,"%c%c%c%c%c", &ctemp5,&ctemp1, &ctemp2, &ctemp3, &ctemp4);
  //printf("reading wc %c  %c\n", ctemp1, ctemp2);
  //if(ctemp1!='A' || ctemp2 != 'm' || ctemp3 != 'a' || ctemp4!='x'){
  //  printf("Wrong syntax at table_wells.tab file (wc)!\n %c %c %c", ctemp1, ctemp2, ctemp3);
  //  exit(1);}
  for(i=0;i<N_BASES;i++)
    for(j=0;j<N_BASES;j++){
      //fscanf(table,"%lf", &(wc_secdih_max[N_BASES*i+j][0*WC_FACES+0]));
      fread(&(wc_secdih_max[N_BASES*i+j][0*WC_FACES+0] ), binnum, 1, outfile);
      //fscanf(table,"%lf", &(wc_secdih_max[N_BASES*i+j][1*WC_FACES+1]));
      fread(&(wc_secdih_max[N_BASES*i+j][1*WC_FACES+1] ), binnum, 1, outfile);
      //fscanf(table,"%lf", &(wc_secdih_max[N_BASES*i+j][2*WC_FACES+2]));
      fread(&(wc_secdih_max[N_BASES*i+j][2*WC_FACES+2] ), binnum, 1, outfile);
      
      //fscanf(table,"%lf", &(wc_secdih_max[N_BASES*i+j][1*WC_FACES+2]));
      fread(&(wc_secdih_max[N_BASES*i+j][1*WC_FACES+2] ), binnum, 1, outfile);
      wc_secdih_max[N_BASES*j+i][2*WC_FACES+1]=wc_secdih_max[N_BASES*i+j][1*WC_FACES+2];
      //fscanf(table,"%lf", &(wc_secdih_max[N_BASES*i+j][1*WC_FACES+0]));
      fread(&(wc_secdih_max[N_BASES*i+j][1*WC_FACES+0] ), binnum, 1, outfile);
      wc_secdih_max[N_BASES*j+i][0*WC_FACES+1]=wc_secdih_max[N_BASES*i+j][1*WC_FACES+0];
      //fscanf(table,"%lf", &(wc_secdih_max[N_BASES*i+j][2*WC_FACES+0]));
      fread(&(wc_secdih_max[N_BASES*i+j][2*WC_FACES+0] ), binnum, 1, outfile);
      wc_secdih_max[N_BASES*j+i][0*WC_FACES+2]=wc_secdih_max[N_BASES*i+j][2*WC_FACES+0];
    }
  //    fscanf(table,"%c%c%c%c%c", &ctemp5,&ctemp1, &ctemp2, &ctemp3, &ctemp4);
  //if(ctemp1!='P' || ctemp2 != 'm' || ctemp3 != 'i' || ctemp4!='n'){
  //printf("Wrong syntax at table_wells.tab file (wc)!\n %c %c %c", ctemp1, ctemp2, ctemp3);
  //exit(1);}
  for(i=0;i<N_BASES;i++)
    for(j=0;j<N_BASES;j++){
      //fscanf(table,"%lf", &(wc_secdih_min_F[N_BASES*i+j][0*WC_FACES+0]));
      fread(&(wc_secdih_min_F[N_BASES*i+j][0*WC_FACES+0] ), binnum, 1, outfile);
      //fscanf(table,"%lf", &(wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+1]));
      fread(&(wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+1] ), binnum, 1, outfile);
      //fscanf(table,"%lf", &(wc_secdih_min_F[N_BASES*i+j][2*WC_FACES+2]));
      fread(&(wc_secdih_min_F[N_BASES*i+j][2*WC_FACES+2] ), binnum, 1, outfile);
      
      //fscanf(table,"%lf", &(wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+2]));
      fread(&(wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+2] ), binnum, 1, outfile);
      wc_secdih_min_F[N_BASES*j+i][2*WC_FACES+1]=wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+2];
      //fscanf(table,"%lf", &(wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+0]));
      fread(&(wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+0] ), binnum, 1, outfile);
      wc_secdih_min_F[N_BASES*j+i][0*WC_FACES+1]=wc_secdih_min_F[N_BASES*i+j][1*WC_FACES+0];
      //fscanf(table,"%lf", &(wc_secdih_min_F[N_BASES*i+j][2*WC_FACES+0]));
      fread(&(wc_secdih_min_F[N_BASES*i+j][2*WC_FACES+0] ), binnum, 1, outfile);
      wc_secdih_min_F[N_BASES*j+i][0*WC_FACES+2]=wc_secdih_min_F[N_BASES*i+j][2*WC_FACES+0];
    }
  //fscanf(table,"%c%c%c%c%c", &ctemp5,&ctemp1, &ctemp2, &ctemp3, &ctemp4);
  //if(ctemp1!='P' || ctemp2 != 'm' || ctemp3 != 'a' || ctemp4!='x'){
  //printf("Wrong syntax at table_wells.tab file (wc)!\n %c %c %c", ctemp1, ctemp2, ctemp3);
  //exit(1);}
  for(i=0;i<N_BASES;i++)
    for(j=0;j<N_BASES;j++){
      //fscanf(table,"%lf", &(wc_secdih_max_F[N_BASES*i+j][0*WC_FACES+0]));
      fread(&(wc_secdih_max_F[N_BASES*i+j][0*WC_FACES+0] ), binnum, 1, outfile);
      //fscanf(table,"%lf", &(wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+1]));
      fread(&(wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+1] ), binnum, 1, outfile);
      //fscanf(table,"%lf", &(wc_secdih_max_F[N_BASES*i+j][2*WC_FACES+2]));
      fread(&(wc_secdih_max_F[N_BASES*i+j][2*WC_FACES+2] ), binnum, 1, outfile);
      
      //fscanf(table,"%lf", &(wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+2]));
      fread(&(wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+2] ), binnum, 1, outfile);
      wc_secdih_max_F[N_BASES*j+i][2*WC_FACES+1]=wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+2];
      //fscanf(table,"%lf", &(wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+0]));
      fread(&(wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+0] ), binnum, 1, outfile);
      wc_secdih_max_F[N_BASES*j+i][0*WC_FACES+1]=wc_secdih_max_F[N_BASES*i+j][1*WC_FACES+0];
      //fscanf(table,"%lf", &(wc_secdih_max_F[N_BASES*i+j][2*WC_FACES+0]));
      fread(&(wc_secdih_max_F[N_BASES*i+j][2*WC_FACES+0] ), binnum, 1, outfile);
      wc_secdih_max_F[N_BASES*j+i][0*WC_FACES+2]=wc_secdih_max_F[N_BASES*i+j][2*WC_FACES+0];
    }
  //    fclose(table);
  //}
  
  /*** BACKBONE ***/
  /* SUGAR - PHOSPHATE - SUGAR */
  //we have to read the nine GLYCOSIDIC states : AA (0) , AH (1) , AS (2) , HA (3) , HH (4) , HS (5) , SA (6) , SH (7) , SS (8)
  //BUT WE USE THE PUCKERS
  //TYPE 0
  //sprintf(tablename, "tab_energs/table_ssB1_p33.tab");
  //if((table=fopen(tablename, "r"))==NULL){
  //printf("No backbone interaction between sugars.\n");
  //table_ssB1_N_33=0;
  //}
  //else{
  //printf("[ANGLE P3 P3] ");
  //fscanf(table,"%d", &(table_ssB1_N_33));
  fread(&(table_ssB1_N_33 ), sizeof(int), 1, outfile);
  //fscanf(table,"%lf", &(table_ssB1_params_33[0]));
  fread(&(table_ssB1_params_33[0] ), binnum, 1, outfile);
  //fscanf(table,"%lf", &(table_ssB1_params_33[1]));
  fread(&(table_ssB1_params_33[1] ), binnum, 1, outfile);
  table_ssB1_33=(double *)malloc(sizeof(double)*table_ssB1_N_33);
  for(ener=0;ener<table_ssB1_N_33;ener++){
    //fscanf(table, "%lf", &(table_ssB1_33[ener]));
    fread(&(table_ssB1_33[ener] ), binnum, 1, outfile);
  }
  //fclose(table);
  //}
  
  //TYPE 32
  //sprintf(tablename, "tab_energs/table_ssB1_p32.tab");
  //FILE *table;
  //if((table=fopen(tablename, "r"))==NULL){
  //printf("No backbone interaction between sugars.\n");
  //table_ssB1_N_32=0;
  //}
  //else{
  //printf("[ANGLE P3 P2] ");
  //fscanf(table,"%d", &(table_ssB1_N_32));
  fread(&(table_ssB1_N_32 ), sizeof(int), 1, outfile);
  //fscanf(table,"%lf", &(table_ssB1_params_32[0]));
  fread(&(table_ssB1_params_32[0] ), binnum, 1, outfile);
  //fscanf(table,"%lf", &(table_ssB1_params_32[1]));
  fread(&(table_ssB1_params_32[1] ), binnum, 1, outfile);
  table_ssB1_32=(double *)malloc(sizeof(double)*table_ssB1_N_32);
  for(ener=0;ener<table_ssB1_N_32;ener++){
    //fscanf(table, "%lf", &(table_ssB1_32[ener]));
    fread(&(table_ssB1_32[ener] ), binnum, 1, outfile);
  }
  //fclose(table);
  //}
  
  //TYPE 23
  //sprintf(tablename, "tab_energs/table_ssB1_p23.tab");
  //if((table=fopen(tablename, "r"))==NULL){
  //printf("No backbone interaction between sugars.\n");
  //table_ssB1_N_23=0;
  //}
  //else{
  //printf("[ANGLE P2 P3] ");
  //fscanf(table,"%d", &(table_ssB1_N_23));
  fread(&(table_ssB1_N_23 ), sizeof(int), 1, outfile);
  //fscanf(table,"%lf", &(table_ssB1_params_23[0]));
  fread(&(table_ssB1_params_23[0] ), binnum, 1, outfile);
  //fscanf(table,"%lf", &(table_ssB1_params_23[1]));
  fread(&(table_ssB1_params_23[1] ), binnum, 1, outfile);
  
  table_ssB1_23=(double *)malloc(sizeof(double)*table_ssB1_N_23);
  for(ener=0;ener<table_ssB1_N_23;ener++){
    //fscanf(table, "%lf", &(table_ssB1_23[ener]));
    fread(&(table_ssB1_23[ener] ), binnum, 1, outfile);
  }
  //fclose(table);
  //}
  //TYPE 22
  //sprintf(tablename, "tab_energs/table_ssB1_p22.tab");
  //if((table=fopen(tablename, "r"))==NULL){
  // printf("No backbone interaction between sugars.\n");
  // table_ssB1_N_22=0;
  //}
  //else{
  //printf("[ANGLE P2 P2] ");
  //fscanf(table,"%d", &(table_ssB1_N_22));
  fread(&(table_ssB1_N_22 ), sizeof(int), 1, outfile);
  //fscanf(table,"%lf", &(table_ssB1_params_22[0]));
  fread(&(table_ssB1_params_22[0] ), binnum, 1, outfile);
  //fscanf(table,"%lf", &(table_ssB1_params_22[1]));
  fread(&(table_ssB1_params_22[1] ), binnum, 1, outfile);
  table_ssB1_22=(double *)malloc(sizeof(double)*table_ssB1_N_22);
  for(ener=0;ener<table_ssB1_N_22;ener++){
    //fscanf(table, "%lf", &(table_ssB1_22[ener]));
    fread(&(table_ssB1_22[ener] ), binnum, 1, outfile);
  }
  //fclose(table);
  //}
  
  /* BASE - PHOSPHATE , INTRA NT */
  //glyc ANTI
  //puck3
  //printf("\nReading INTRA-BP interactions...\n");
  table_bpI_A3=(double **)malloc(sizeof(double *)*mc_n_types);
  for(i=0;i<N_BASES;i++){
    typ_ind=i;
    //sprintf(tablename, "tab_energs/table_bpI_%d_gA3.tab", i);
    //if((table=fopen(tablename, "r"))==NULL){
    // printf("No backbone intra BP interaction for %d, glycosidic conformation ANTI.\n", i);
    // table_bpI_N_A3[typ_ind][0]=-1;
    // table_bpI_N_A3[typ_ind][1]=-1;
    // table_bpI_N_A3[typ_ind][2]=-1;
    // exit(ERR_INPUT);
    //}
    //else{
    //printf("[%d ANTI P3] ", i);
    //fscanf(table,"%d%d%d", &(table_bpI_N_A3[typ_ind][0]), &(table_bpI_N_A3[typ_ind][1]), &(table_bpI_N_A3[typ_ind][2])); //NX, NY, NZ
    fread(&(table_bpI_N_A3[typ_ind][0] ), sizeof(int), 1, outfile);
    fread(&(table_bpI_N_A3[typ_ind][1] ), sizeof(int), 1, outfile);
    fread(&(table_bpI_N_A3[typ_ind][2] ), sizeof(int), 1, outfile);
    //fscanf(table,"%lf%lf%lf", &(table_bpI_params_A3[typ_ind][0][0]), &(table_bpI_params_A3[typ_ind][1][0]), &(table_bpI_params_A3[typ_ind][2][0]));//DX, DY, DZ
    fread(&(table_bpI_params_A3[typ_ind][0][0] ), binnum, 1, outfile);
    fread(&(table_bpI_params_A3[typ_ind][1][0] ), binnum, 1, outfile);
    fread(&(table_bpI_params_A3[typ_ind][2][0] ), binnum, 1, outfile);
    //fscanf(table,"%lf%lf%lf", &(table_bpI_params_A3[typ_ind][0][1]), &(table_bpI_params_A3[typ_ind][1][1]), &(table_bpI_params_A3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
    fread(&(table_bpI_params_A3[typ_ind][0][1] ), binnum, 1, outfile);
    fread(&(table_bpI_params_A3[typ_ind][1][1] ), binnum, 1, outfile);
    fread(&(table_bpI_params_A3[typ_ind][2][1] ), binnum, 1, outfile);
    ntot=table_bpI_N_A3[typ_ind][0]*table_bpI_N_A3[typ_ind][1]*table_bpI_N_A3[typ_ind][2];
    table_bpI_A3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
    for(ener=0;ener<ntot;ener++){
      //fscanf(table, "%lf", &(table_bpI_A3[typ_ind][ener]));
      fread(&(table_bpI_A3[typ_ind][ener] ), binnum, 1, outfile);
    }
    //fclose(table);
    //}
  }
      
  //puck2
  table_bpI_A2=(double **)malloc(sizeof(double *)*mc_n_types);
  for(i=0;i<N_BASES;i++){
    typ_ind=i;
    //sprintf(tablename, "tab_energs/table_bpI_%d_gA2.tab", i);
    //if((table=fopen(tablename, "r"))==NULL){
    //printf("No backbone intra BP interaction for %d, glycosidic conformation ANTI.\n", i);
    //table_bpI_N_A2[typ_ind][0]=-1;
    //table_bpI_N_A2[typ_ind][1]=-1;
    //table_bpI_N_A2[typ_ind][2]=-1;
    //}
    //else{
    //printf("[%d ANTI P2] ", i);
    //fscanf(table,"%d%d%d", &(table_bpI_N_A2[typ_ind][0]), &(table_bpI_N_A2[typ_ind][1]), &(table_bpI_N_A2[typ_ind][2])); //NX, NY, NZ
    fread(&(table_bpI_N_A2[typ_ind][0] ), sizeof(int), 1, outfile);
    fread(&(table_bpI_N_A2[typ_ind][1] ), sizeof(int), 1, outfile);
    fread(&(table_bpI_N_A2[typ_ind][2] ), sizeof(int), 1, outfile);
    //fscanf(table,"%lf%lf%lf", &(table_bpI_params_A2[typ_ind][0][0]), &(table_bpI_params_A2[typ_ind][1][0]), &(table_bpI_params_A2[typ_ind][2][0]));//DX, DY, DZ
    fread(&(table_bpI_params_A2[typ_ind][0][0] ), binnum, 1, outfile);
    fread(&(table_bpI_params_A2[typ_ind][1][0] ), binnum, 1, outfile);
    fread(&(table_bpI_params_A2[typ_ind][2][0] ), binnum, 1, outfile);
    //fscanf(table,"%lf%lf%lf", &(table_bpI_params_A2[typ_ind][0][1]), &(table_bpI_params_A2[typ_ind][1][1]), &(table_bpI_params_A2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
    fread(&(table_bpI_params_A2[typ_ind][0][1] ), binnum, 1, outfile);
    fread(&(table_bpI_params_A2[typ_ind][1][1] ), binnum, 1, outfile);
    fread(&(table_bpI_params_A2[typ_ind][2][1] ), binnum, 1, outfile);
    ntot=table_bpI_N_A2[typ_ind][0]*table_bpI_N_A2[typ_ind][1]*table_bpI_N_A2[typ_ind][2];
    table_bpI_A2[typ_ind]=(double *)malloc(sizeof(double)*ntot);
    for(ener=0;ener<ntot;ener++){
      //fscanf(table, "%lf", &(table_bpI_A2[typ_ind][ener]));
      fread(&(table_bpI_A2[typ_ind][ener] ), binnum, 1, outfile);
    }
    //fclose(table);
    //}
  }
  
  //glyc HIGH ANTI
  table_bpI_H3=(double **)malloc(sizeof(double *)*mc_n_types);
  for(i=0;i<N_BASES;i++){
    typ_ind=i;
    //sprintf(tablename, "tab_energs/table_bpI_%d_gH3.tab", i);
    //if((table=fopen(tablename, "r"))==NULL){
    //printf("No backbone intra SP interaction for %d, glycosidic conformation HIGH ANTI.\n", i);
    //table_bpI_N_H3[typ_ind][0]=-1;
    //table_bpI_N_H3[typ_ind][1]=-1;
    //table_bpI_N_H3[typ_ind][2]=-1;
    //}
    //else{
    //printf("[%d HIGH ANTI P3] ", i);
    //fscanf(table,"%d%d%d", &(table_bpI_N_H3[typ_ind][0]), &(table_bpI_N_H3[typ_ind][1]), &(table_bpI_N_H3[typ_ind][2])); //NX, NY, NZ
    fread(&(table_bpI_N_H3[typ_ind][0] ), sizeof(int), 1, outfile);
    fread(&(table_bpI_N_H3[typ_ind][1] ), sizeof(int), 1, outfile);
    fread(&(table_bpI_N_H3[typ_ind][2] ), sizeof(int), 1, outfile);
    //fscanf(table,"%lf%lf%lf", &(table_bpI_params_H3[typ_ind][0][0]), &(table_bpI_params_H3[typ_ind][1][0]), &(table_bpI_params_H3[typ_ind][2][0]));//DX, DY, DZ
    fread(&(table_bpI_params_H3[typ_ind][0][0] ), binnum, 1, outfile);
    fread(&(table_bpI_params_H3[typ_ind][1][0] ), binnum, 1, outfile);
    fread(&(table_bpI_params_H3[typ_ind][2][0] ), binnum, 1, outfile);
    //fscanf(table,"%lf%lf%lf", &(table_bpI_params_H3[typ_ind][0][1]), &(table_bpI_params_H3[typ_ind][1][1]), &(table_bpI_params_H3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
    fread(&(table_bpI_params_H3[typ_ind][0][1] ), binnum, 1, outfile);
    fread(&(table_bpI_params_H3[typ_ind][1][1] ), binnum, 1, outfile);
    fread(&(table_bpI_params_H3[typ_ind][2][1] ), binnum, 1, outfile);
    
    ntot=table_bpI_N_H3[typ_ind][0]*table_bpI_N_H3[typ_ind][1]*table_bpI_N_H3[typ_ind][2];
    table_bpI_H3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
    for(ener=0;ener<ntot;ener++){
      //fscanf(table, "%lf", &(table_bpI_H3[typ_ind][ener]));
      fread(&(table_bpI_H3[typ_ind][ener] ), binnum, 1, outfile);
    }
    //fclose(table);
    //}
  }
  table_bpI_H2=(double **)malloc(sizeof(double *)*mc_n_types);
  for(i=0;i<N_BASES;i++){
    typ_ind=i;
    //sprintf(tablename, "tab_energs/table_bpI_%d_gH2.tab", i);
    //if((table=fopen(tablename, "r"))==NULL){
    //printf("No backbone intra SP interaction for %d, glycosidic conformation HIGH ANTI.\n", i);
    //table_bpI_N_H2[typ_ind][0]=-1;
    //table_bpI_N_H2[typ_ind][1]=-1;
    //table_bpI_N_H2[typ_ind][2]=-1;
    //}
    //else{
    //printf("[%d HIGH ANTI P2] ", i);
    //fscanf(table,"%d%d%d", &(table_bpI_N_H2[typ_ind][0]), &(table_bpI_N_H2[typ_ind][1]), &(table_bpI_N_H2[typ_ind][2])); //NX, NY, NZ
    fread(&(table_bpI_N_H2[typ_ind][0] ), sizeof(int), 1, outfile);
    fread(&(table_bpI_N_H2[typ_ind][1] ), sizeof(int), 1, outfile);
    fread(&(table_bpI_N_H2[typ_ind][2] ), sizeof(int), 1, outfile);
    //fscanf(table,"%lf%lf%lf", &(table_bpI_params_H2[typ_ind][0][0]), &(table_bpI_params_H2[typ_ind][1][0]), &(table_bpI_params_H2[typ_ind][2][0]));//DX, DY, DZ
    fread(&(table_bpI_params_H2[typ_ind][0][0] ), binnum, 1, outfile);
    fread(&(table_bpI_params_H2[typ_ind][1][0] ), binnum, 1, outfile);
    fread(&(table_bpI_params_H2[typ_ind][2][0] ), binnum, 1, outfile);
    //fscanf(table,"%lf%lf%lf", &(table_bpI_params_H2[typ_ind][0][1]), &(table_bpI_params_H2[typ_ind][1][1]), &(table_bpI_params_H2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
    fread(&(table_bpI_params_H2[typ_ind][0][1] ), binnum, 1, outfile);
    fread(&(table_bpI_params_H2[typ_ind][1][1] ), binnum, 1, outfile);
    fread(&(table_bpI_params_H2[typ_ind][2][1] ), binnum, 1, outfile);
    ntot=table_bpI_N_H2[typ_ind][0]*table_bpI_N_H2[typ_ind][1]*table_bpI_N_H2[typ_ind][2];
    table_bpI_H2[typ_ind]=(double *)malloc(sizeof(double)*ntot);
    for(ener=0;ener<ntot;ener++){
      //fscanf(table, "%lf", &(table_bpI_H2[typ_ind][ener]));
      fread(&(table_bpI_H2[typ_ind][ener] ), binnum, 1, outfile);
    }
    //fclose(table);
    //}
  }
  //glyc SYN
  table_bpI_S3=(double **)malloc(sizeof(double *)*mc_n_types);
  for(i=0;i<N_BASES;i++){
    typ_ind=i;
    //sprintf(tablename, "tab_energs/table_bpI_%d_gS3.tab", i);
    //if((table=fopen(tablename, "r"))==NULL){
    //printf("No backbone intra SP interaction for %d, glycosidic conformation SYN.\n", i);
    //table_bpI_N_S3[typ_ind][0]=-1;
    //table_bpI_N_S3[typ_ind][1]=-1;
    //table_bpI_N_S3[typ_ind][2]=-1;
    //}
    //else{
    //printf("[%d SYN A3] ", i);
    //fscanf(table,"%d%d%d", &(table_bpI_N_S3[typ_ind][0]), &(table_bpI_N_S3[typ_ind][1]), &(table_bpI_N_S3[typ_ind][2])); //NX, NY, NZ
    fread(&(table_bpI_N_S3[typ_ind][0] ), sizeof(int), 1, outfile);
    fread(&(table_bpI_N_S3[typ_ind][1] ), sizeof(int), 1, outfile);
    fread(&(table_bpI_N_S3[typ_ind][2] ), sizeof(int), 1, outfile);
    //fscanf(table,"%lf%lf%lf", &(table_bpI_params_S3[typ_ind][0][0]), &(table_bpI_params_S3[typ_ind][1][0]), &(table_bpI_params_S3[typ_ind][2][0]));//DX, DY, DZ
    fread(&(table_bpI_params_S3[typ_ind][0][0] ), binnum, 1, outfile);
    fread(&(table_bpI_params_S3[typ_ind][1][0] ), binnum, 1, outfile);
    fread(&(table_bpI_params_S3[typ_ind][2][0] ), binnum, 1, outfile);
    //fscanf(table,"%lf%lf%lf", &(table_bpI_params_S3[typ_ind][0][1]), &(table_bpI_params_S3[typ_ind][1][1]), &(table_bpI_params_S3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
    fread(&(table_bpI_params_S3[typ_ind][0][1] ), binnum, 1, outfile);
    fread(&(table_bpI_params_S3[typ_ind][1][1] ), binnum, 1, outfile);
    fread(&(table_bpI_params_S3[typ_ind][2][1] ), binnum, 1, outfile);
    ntot=table_bpI_N_S3[typ_ind][0]*table_bpI_N_S3[typ_ind][1]*table_bpI_N_S3[typ_ind][2];
    table_bpI_S3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
    for(ener=0;ener<ntot;ener++){
      //fscanf(table, "%lf", &(table_bpI_S3[typ_ind][ener]));
      fread(&(table_bpI_S3[typ_ind][ener] ), binnum, 1, outfile);
    }
    //fclose(table);
    //}
  }
  table_bpI_S2=(double **)malloc(sizeof(double *)*mc_n_types);
  for(i=0;i<N_BASES;i++){
    typ_ind=i;
    //sprintf(tablename, "tab_energs/table_bpI_%d_gS2.tab", i);
    //if((table=fopen(tablename, "r"))==NULL){
    //printf("No backbone intra SP interaction for %d, glycosidic conformation SYN.\n", i);
    //table_bpI_N_S2[typ_ind][0]=-1;
    //table_bpI_N_S2[typ_ind][1]=-1;
    //table_bpI_N_S2[typ_ind][2]=-1;
    //}
    //else{
    //printf("[%d SYN P2] ", i);
    //fscanf(table,"%d%d%d", &(table_bpI_N_S2[typ_ind][0]), &(table_bpI_N_S2[typ_ind][1]), &(table_bpI_N_S2[typ_ind][2])); //NX, NY, NZ
    fread(&(table_bpI_N_S2[typ_ind][0] ), sizeof(int), 1, outfile);
    fread(&(table_bpI_N_S2[typ_ind][1] ), sizeof(int), 1, outfile);
    fread(&(table_bpI_N_S2[typ_ind][2] ), sizeof(int), 1, outfile);
    //fscanf(table,"%lf%lf%lf", &(table_bpI_params_S2[typ_ind][0][0]), &(table_bpI_params_S2[typ_ind][1][0]), &(table_bpI_params_S2[typ_ind][2][0]));//DX, DY, DZ
    fread(&(table_bpI_params_S2[typ_ind][0][0] ), binnum, 1, outfile);
    fread(&(table_bpI_params_S2[typ_ind][1][0] ), binnum, 1, outfile);
    fread(&(table_bpI_params_S2[typ_ind][2][0] ), binnum, 1, outfile);
    //fscanf(table,"%lf%lf%lf", &(table_bpI_params_S2[typ_ind][0][1]), &(table_bpI_params_S2[typ_ind][1][1]), &(table_bpI_params_S2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
    fread(&(table_bpI_params_S2[typ_ind][0][1] ), binnum, 1, outfile);
    fread(&(table_bpI_params_S2[typ_ind][1][1] ), binnum, 1, outfile);
    fread(&(table_bpI_params_S2[typ_ind][2][1] ), binnum, 1, outfile);
    ntot=table_bpI_N_S2[typ_ind][0]*table_bpI_N_S2[typ_ind][1]*table_bpI_N_S2[typ_ind][2];
    table_bpI_S2[typ_ind]=(double *)malloc(sizeof(double)*ntot);
    for(ener=0;ener<ntot;ener++){
      //fscanf(table, "%lf", &(table_bpI_S2[typ_ind][ener]));
      fread(&(table_bpI_S2[typ_ind][ener] ), binnum, 1, outfile);
    }
    //fclose(table);
    //}
  }
  
  /* SUGAR - PHOSPHATE , INTER NT */
  //GLYC ANTI
  //printf("\nReading INTER-BP interactions...\n");
  table_bpB_A3=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //sprintf(tablename, "tab_energs/table_bpB_%d%d_gA3.tab", i, j);
      //if((table=fopen(tablename, "r"))==NULL){
      //printf("No backbone inter SP interaction for %d and %d.\n", i,j);
      //table_bpB_N_A3[typ_ind][0]=-1;
      //table_bpB_N_A3[typ_ind][1]=-1;
      //table_bpB_N_A3[typ_ind][2]=-1;
      //}
      //else{
      //printf("[%d %d ANTI P3] ", i,j);
      //fscanf(table,"%d%d%d", &(table_bpB_N_A3[typ_ind][0]), &(table_bpB_N_A3[typ_ind][1]), &(table_bpB_N_A3[typ_ind][2])); //NX, NY, NZ
      fread(&(table_bpB_N_A3[typ_ind][0] ), sizeof(int), 1, outfile);
      fread(&(table_bpB_N_A3[typ_ind][1] ), sizeof(int), 1, outfile);
      fread(&(table_bpB_N_A3[typ_ind][2] ), sizeof(int), 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_bpB_params_A3[typ_ind][0][0]), &(table_bpB_params_A3[typ_ind][1][0]), &(table_bpB_params_A3[typ_ind][2][0]));//DX, DY, DZ
      fread(&(table_bpB_params_A3[typ_ind][0][0] ), binnum, 1, outfile);
      fread(&(table_bpB_params_A3[typ_ind][1][0] ), binnum, 1, outfile);
      fread(&(table_bpB_params_A3[typ_ind][2][0] ), binnum, 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_bpB_params_A3[typ_ind][0][1]), &(table_bpB_params_A3[typ_ind][1][1]), &(table_bpB_params_A3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fread(&(table_bpB_params_A3[typ_ind][0][1] ), binnum, 1, outfile);
      fread(&(table_bpB_params_A3[typ_ind][1][1] ), binnum, 1, outfile);
      fread(&(table_bpB_params_A3[typ_ind][2][1] ), binnum, 1, outfile);
      ntot=table_bpB_N_A3[typ_ind][0]*table_bpB_N_A3[typ_ind][1]*table_bpB_N_A3[typ_ind][2];
      table_bpB_A3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	//fscanf(table, "%lf", &(table_bpB_A3[typ_ind][ener]));
	fread(&(table_bpB_A3[typ_ind][ener] ), binnum, 1, outfile);
      }
      //fclose(table);
      //}
    }
  }
  table_bpB_A2=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //sprintf(tablename, "tab_energs/table_bpB_%d%d_gA2.tab", i, j);
      //if((table=fopen(tablename, "r"))==NULL){
      //printf("No backbone inter SP interaction for %d and %d.\n", i,j);
      //table_bpB_N_A2[typ_ind][0]=-1;
      //table_bpB_N_A2[typ_ind][1]=-1;
      //table_bpB_N_A2[typ_ind][2]=-1;
      //}
      //else{
      //printf("[%d %d ANTI P2] ", i,j);
      //fscanf(table,"%d%d%d", &(table_bpB_N_A2[typ_ind][0]), &(table_bpB_N_A2[typ_ind][1]), &(table_bpB_N_A2[typ_ind][2])); //NX, NY, NZ
      fread(&(table_bpB_N_A2[typ_ind][0]) , sizeof(int), 1, outfile);
      fread(&(table_bpB_N_A2[typ_ind][1]) , sizeof(int), 1, outfile);
      fread(&(table_bpB_N_A2[typ_ind][2]) , sizeof(int), 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_bpB_params_A2[typ_ind][0][0]), &(table_bpB_params_A2[typ_ind][1][0]), &(table_bpB_params_A2[typ_ind][2][0]));//DX, DY, DZ
      fread(&(table_bpB_params_A2[typ_ind][0][0] ) , binnum, 1, outfile);
      fread(&(table_bpB_params_A2[typ_ind][1][0] ) , binnum, 1, outfile);
      fread(&(table_bpB_params_A2[typ_ind][2][0] ) , binnum, 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_bpB_params_A2[typ_ind][0][1]), &(table_bpB_params_A2[typ_ind][1][1]), &(table_bpB_params_A2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fread(&(table_bpB_params_A2[typ_ind][0][1] ) , binnum, 1, outfile);
      fread(&(table_bpB_params_A2[typ_ind][1][1] ) , binnum, 1, outfile);
      fread(&(table_bpB_params_A2[typ_ind][2][1] ) , binnum, 1, outfile);
      ntot=table_bpB_N_A2[typ_ind][0]*table_bpB_N_A2[typ_ind][1]*table_bpB_N_A2[typ_ind][2];
      table_bpB_A2[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	//fscanf(table, "%lf", &(table_bpB_A2[typ_ind][ener]));
	fread(&(table_bpB_A2[typ_ind][ener] ) , binnum, 1, outfile);
      }
      //fclose(table);
      //}
    }
  }
      
  //GLYC HIGH ANTI
  table_bpB_H3=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //sprintf(tablename, "tab_energs/table_bpB_%d%d_gH3.tab", i, j);
      //if((table=fopen(tablename, "r"))==NULL){
      //printf("No backbone inter SP interaction for %d and %d.\n", i,j);
      //table_bpB_N_H3[typ_ind][0]=-1;
      //table_bpB_N_H3[typ_ind][1]=-1;
      //table_bpB_N_H3[typ_ind][2]=-1;
      //}
      //else{
	//printf("[%d %d HIGH ANTI P3] ", i,j);
	//fscanf(table,"%d%d%d", &(table_bpB_N_H3[typ_ind][0]), &(table_bpB_N_H3[typ_ind][1]), &(table_bpB_N_H3[typ_ind][2])); //NX, NY, NZ
	fread(&(table_bpB_N_H3[typ_ind][0] ) , sizeof(int), 1, outfile);
	fread(&(table_bpB_N_H3[typ_ind][1] ) , sizeof(int), 1, outfile);
	fread(&(table_bpB_N_H3[typ_ind][2] ) , sizeof(int), 1, outfile);
	//fscanf(table,"%lf%lf%lf", &(table_bpB_params_H3[typ_ind][0][0]), &(table_bpB_params_H3[typ_ind][1][0]), &(table_bpB_params_H3[typ_ind][2][0]));//DX, DY, DZ
	fread(&(table_bpB_params_H3[typ_ind][0][0] ) , binnum, 1, outfile);
	fread(&(table_bpB_params_H3[typ_ind][1][0] ) , binnum, 1, outfile);
	fread(&(table_bpB_params_H3[typ_ind][2][0] ) , binnum, 1, outfile);
	//fscanf(table,"%lf%lf%lf", &(table_bpB_params_H3[typ_ind][0][1]), &(table_bpB_params_H3[typ_ind][1][1]), &(table_bpB_params_H3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
	fread(&(table_bpB_params_H3[typ_ind][0][1] ) , binnum, 1, outfile);
	fread(&(table_bpB_params_H3[typ_ind][1][1] ) , binnum, 1, outfile);
	fread(&(table_bpB_params_H3[typ_ind][2][1] ) , binnum, 1, outfile);
	ntot=table_bpB_N_H3[typ_ind][0]*table_bpB_N_H3[typ_ind][1]*table_bpB_N_H3[typ_ind][2];
	table_bpB_H3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
	for(ener=0;ener<ntot;ener++){
	  //fscanf(table, "%lf", &(table_bpB_H3[typ_ind][ener]));
	  fread(&(table_bpB_H3[typ_ind][ener] ) , binnum, 1, outfile);
	}
	//    fclose(table);
	//}
    }
  }
  table_bpB_H2=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //sprintf(tablename, "tab_energs/table_bpB_%d%d_gH2.tab", i, j);
      //if((table=fopen(tablename, "r"))==NULL){
      //printf("No backbone inter SP interaction for %d and %d.\n", i,j);
      //table_bpB_N_H2[typ_ind][0]=-1;
      //table_bpB_N_H2[typ_ind][1]=-1;
      //table_bpB_N_H2[typ_ind][2]=-1;
      //}
      //else{
      //printf("[%d %d HIGH ANTI P2] ", i,j);
      //fscanf(table,"%d%d%d", &(table_bpB_N_H2[typ_ind][0]), &(table_bpB_N_H2[typ_ind][1]), &(table_bpB_N_H2[typ_ind][2])); //NX, NY, NZ
      fread(&(table_bpB_N_H2[typ_ind][0] ) , sizeof(int), 1, outfile);
      fread(&(table_bpB_N_H2[typ_ind][1] ) , sizeof(int), 1, outfile);
      fread(&(table_bpB_N_H2[typ_ind][2] ) , sizeof(int), 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_bpB_params_H2[typ_ind][0][0]), &(table_bpB_params_H2[typ_ind][1][0]), &(table_bpB_params_H2[typ_ind][2][0]));//DX, DY, DZ
      fread(&(table_bpB_params_H2[typ_ind][0][0] ) , binnum, 1, outfile);
      fread(&(table_bpB_params_H2[typ_ind][1][0] ) , binnum, 1, outfile);
      fread(&(table_bpB_params_H2[typ_ind][2][0] ) , binnum, 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_bpB_params_H2[typ_ind][0][1]), &(table_bpB_params_H2[typ_ind][1][1]), &(table_bpB_params_H2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fread(&(table_bpB_params_H2[typ_ind][0][1] ) , binnum, 1, outfile);
      fread(&(table_bpB_params_H2[typ_ind][1][1] ) , binnum, 1, outfile);
      fread(&(table_bpB_params_H2[typ_ind][2][1] ) , binnum, 1, outfile);
      ntot=table_bpB_N_H2[typ_ind][0]*table_bpB_N_H2[typ_ind][1]*table_bpB_N_H2[typ_ind][2];
      table_bpB_H2[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	//fscanf(table, "%lf", &(table_bpB_H2[typ_ind][ener]));
	fread(&(table_bpB_H2[typ_ind][ener] ) , binnum, 1, outfile);
      }
      //fclose(table);
      //}
    }
  }
  
  //GLYC SYN
  table_bpB_S3=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //sprintf(tablename, "tab_energs/table_bpB_%d%d_gS3.tab", i, j);
      //if((table=fopen(tablename, "r"))==NULL){
      //printf("No backbone inter SP interaction for %d and %d.\n", i,j);
      //table_bpB_N_S3[typ_ind][0]=-1;
      //table_bpB_N_S3[typ_ind][1]=-1;
      //table_bpB_N_S3[typ_ind][2]=-1;
      //}
      //else{
      //printf("[%d %d SYN P3] ", i,j);
      //fscanf(table,"%d%d%d", &(table_bpB_N_S3[typ_ind][0]), &(table_bpB_N_S3[typ_ind][1]), &(table_bpB_N_S3[typ_ind][2])); //NX, NY, NZ
      fread(&(table_bpB_N_S3[typ_ind][0] ) , sizeof(int), 1, outfile);
      fread(&(table_bpB_N_S3[typ_ind][1] ) , sizeof(int), 1, outfile);
      fread(&(table_bpB_N_S3[typ_ind][2] ) , sizeof(int), 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_bpB_params_S3[typ_ind][0][0]), &(table_bpB_params_S3[typ_ind][1][0]), &(table_bpB_params_S3[typ_ind][2][0]));//DX, DY, DZ
      fread(&(table_bpB_params_S3[typ_ind][0][0] ) , binnum, 1, outfile);
      fread(&(table_bpB_params_S3[typ_ind][1][0] ) , binnum, 1, outfile);
      fread(&(table_bpB_params_S3[typ_ind][2][0] ) , binnum, 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_bpB_params_S3[typ_ind][0][1]), &(table_bpB_params_S3[typ_ind][1][1]), &(table_bpB_params_S3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fread(&(table_bpB_params_S3[typ_ind][0][1] ) , binnum, 1, outfile);
      fread(&(table_bpB_params_S3[typ_ind][1][1] ) , binnum, 1, outfile);
      fread(&(table_bpB_params_S3[typ_ind][2][1] ) , binnum, 1, outfile);
      ntot=table_bpB_N_S3[typ_ind][0]*table_bpB_N_S3[typ_ind][1]*table_bpB_N_S3[typ_ind][2];
      table_bpB_S3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	//fscanf(table, "%lf", &(table_bpB_S3[typ_ind][ener]));
	fread(&(table_bpB_S3[typ_ind][ener] ) , binnum, 1, outfile);
      }
      //fclose(table);
      //}
    }
  }
  table_bpB_S2=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //sprintf(tablename, "tab_energs/table_bpB_%d%d_gS2.tab", i, j);
      //if((table=fopen(tablename, "r"))==NULL){
      //printf("No backbone inter SP interaction for %d and %d.\n", i,j);
      //table_bpB_N_S2[typ_ind][0]=-1;
      //table_bpB_N_S2[typ_ind][1]=-1;
      //table_bpB_N_S2[typ_ind][2]=-1;
      //}
      //else{
      //printf("[%d %d SYN P2] ", i,j);
      //fscanf(table,"%d%d%d", &(table_bpB_N_S2[typ_ind][0]), &(table_bpB_N_S2[typ_ind][1]), &(table_bpB_N_S2[typ_ind][2])); //NX, NY, NZ
      fread(&(table_bpB_N_S2[typ_ind][0] ) , sizeof(int), 1, outfile);
      fread(&(table_bpB_N_S2[typ_ind][1] ) , sizeof(int), 1, outfile);
      fread(&(table_bpB_N_S2[typ_ind][2] ) , sizeof(int), 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_bpB_params_S2[typ_ind][0][0]), &(table_bpB_params_S2[typ_ind][1][0]), &(table_bpB_params_S2[typ_ind][2][0]));//DX, DY, DZ
      fread(&(table_bpB_params_S2[typ_ind][0][0] ) , binnum, 1, outfile);
      fread(&(table_bpB_params_S2[typ_ind][1][0] ) , binnum, 1, outfile);
      fread(&(table_bpB_params_S2[typ_ind][2][0] ) , binnum, 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_bpB_params_S2[typ_ind][0][1]), &(table_bpB_params_S2[typ_ind][1][1]), &(table_bpB_params_S2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fread(&(table_bpB_params_S2[typ_ind][0][1] ) , binnum, 1, outfile);
      fread(&(table_bpB_params_S2[typ_ind][1][1] ) , binnum, 1, outfile);
      fread(&(table_bpB_params_S2[typ_ind][2][1] ) , binnum, 1, outfile);
      ntot=table_bpB_N_S2[typ_ind][0]*table_bpB_N_S2[typ_ind][1]*table_bpB_N_S2[typ_ind][2];
      table_bpB_S2[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	//fscanf(table, "%lf", &(table_bpB_S2[typ_ind][ener]));
	fread(&(table_bpB_S2[typ_ind][ener] ) , binnum, 1, outfile);
      }
      //fclose(table);
      //}
    }
  }
  
  //printf("\nReading STACKING interactions...\n");
  /* STACKING - BONDED */
  //s35
  table_nnB_0_s35=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //sprintf(tablename, "tab_energs/table_nnB_%d%d_0_s35.tab", i,j);
      //if((table=fopen(tablename, "r"))==NULL){
      //printf("No stacking interaction between %d and %d.\n", i,j);
      //table_nnB_N_0_s35[typ_ind][0]=-1;
      //table_nnB_N_0_s35[typ_ind][1]=-1;
      //table_nnB_N_0_s35[typ_ind][2]=-1;
      //}
      //else{
      //printf("[%d %d s35] ", i,j);
      //fscanf(table,"%d%d%d", &(table_nnB_N_0_s35[typ_ind][0]), &(table_nnB_N_0_s35[typ_ind][1]), &(table_nnB_N_0_s35[typ_ind][2])); //NX, NY, NZ
      fread(&(table_nnB_N_0_s35[typ_ind][0] ) , sizeof(int), 1, outfile);
      fread(&(table_nnB_N_0_s35[typ_ind][1] ) , sizeof(int), 1, outfile);
      fread(&(table_nnB_N_0_s35[typ_ind][2] ) , sizeof(int), 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s35[typ_ind][0][0]), &(table_nnB_params_0_s35[typ_ind][1][0]), &(table_nnB_params_0_s35[typ_ind][2][0]));//DX, DY, DZ
      fread(&(table_nnB_params_0_s35[typ_ind][0][0] ) , binnum, 1, outfile);
      fread(&(table_nnB_params_0_s35[typ_ind][1][0] ) , binnum, 1, outfile);
      fread(&(table_nnB_params_0_s35[typ_ind][2][0] ) , binnum, 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s35[typ_ind][0][1]), &(table_nnB_params_0_s35[typ_ind][1][1]), &(table_nnB_params_0_s35[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fread(&(table_nnB_params_0_s35[typ_ind][0][1] ) , binnum, 1, outfile);
      fread(&(table_nnB_params_0_s35[typ_ind][1][1] ) , binnum, 1, outfile);
      fread(&(table_nnB_params_0_s35[typ_ind][2][1] ) , binnum, 1, outfile);
      ntot=table_nnB_N_0_s35[typ_ind][0]*table_nnB_N_0_s35[typ_ind][1]*table_nnB_N_0_s35[typ_ind][2];
      table_nnB_0_s35[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	//fscanf(table, "%lf", &(table_nnB_0_s35[typ_ind][ener]));
	fread(&(table_nnB_0_s35[typ_ind][ener] ) , binnum, 1, outfile);
      }
      //fclose(table);
      //}
    }
  }
  //S53
  table_nnB_0_s53=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //sprintf(tablename, "tab_energs/table_nnB_%d%d_0_s53.tab", i,j);
      //if((table=fopen(tablename, "r"))==NULL){
      //printf("No stacking interaction between %d and %d.\n", i,j);
      //table_nnB_N_0_s53[typ_ind][0]=-1;
      //table_nnB_N_0_s53[typ_ind][1]=-1;
      //table_nnB_N_0_s53[typ_ind][2]=-1;
      //}
      //else{
      //printf("[%d %d s53] ", i,j);
      //fscanf(table,"%d%d%d", &(table_nnB_N_0_s53[typ_ind][0]), &(table_nnB_N_0_s53[typ_ind][1]), &(table_nnB_N_0_s53[typ_ind][2])); //NX, NY, NZ
      fread(&(table_nnB_N_0_s53[typ_ind][0] ) , sizeof(int), 1, outfile);
      fread(&(table_nnB_N_0_s53[typ_ind][1] ) , sizeof(int), 1, outfile);
      fread(&(table_nnB_N_0_s53[typ_ind][2] ) , sizeof(int), 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s53[typ_ind][0][0]), &(table_nnB_params_0_s53[typ_ind][1][0]), &(table_nnB_params_0_s53[typ_ind][2][0]));//DX, DY, DZ
      fread(&(table_nnB_params_0_s53[typ_ind][0][0] ) , binnum, 1, outfile);
      fread(&(table_nnB_params_0_s53[typ_ind][1][0] ) , binnum, 1, outfile);
      fread(&(table_nnB_params_0_s53[typ_ind][2][0] ) , binnum, 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s53[typ_ind][0][1]), &(table_nnB_params_0_s53[typ_ind][1][1]), &(table_nnB_params_0_s53[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fread(&(table_nnB_params_0_s53[typ_ind][0][1] ) , binnum, 1, outfile);
      fread(&(table_nnB_params_0_s53[typ_ind][1][1] ) , binnum, 1, outfile);
      fread(&(table_nnB_params_0_s53[typ_ind][2][1] ) , binnum, 1, outfile);
      ntot=table_nnB_N_0_s53[typ_ind][0]*table_nnB_N_0_s53[typ_ind][1]*table_nnB_N_0_s53[typ_ind][2];
      table_nnB_0_s53[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	//fscanf(table, "%lf", &(table_nnB_0_s53[typ_ind][ener]));
	fread(&(table_nnB_0_s53[typ_ind][ener] ) , binnum, 1, outfile);
      }
      //fclose(table);
      //}
    }
  }
  //S33
  table_nnB_0_s33=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //sprintf(tablename, "tab_energs/table_nnB_%d%d_0_s33.tab", i,j);
      //if((table=fopen(tablename, "r"))==NULL){
      //printf("No stacking interaction between %d and %d.\n", i,j);
      //table_nnB_N_0_s33[typ_ind][0]=-1;
      //table_nnB_N_0_s33[typ_ind][1]=-1;
      //table_nnB_N_0_s33[typ_ind][2]=-1;
      //}
      //      else{
      //printf("[%d %d s33] ", i,j);
      //fscanf(table,"%d%d%d", &(table_nnB_N_0_s33[typ_ind][0]), &(table_nnB_N_0_s33[typ_ind][1]), &(table_nnB_N_0_s33[typ_ind][2])); //NX, NY, NZ
      fread(&(table_nnB_N_0_s33[typ_ind][0] ) , sizeof(int), 1, outfile);
      fread(&(table_nnB_N_0_s33[typ_ind][1] ) , sizeof(int), 1, outfile);
      fread(&(table_nnB_N_0_s33[typ_ind][2] ) , sizeof(int), 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s33[typ_ind][0][0]), &(table_nnB_params_0_s33[typ_ind][1][0]), &(table_nnB_params_0_s33[typ_ind][2][0]));//DX, DY, DZ
      fread(&(table_nnB_params_0_s33[typ_ind][0][0] ) , binnum, 1, outfile);
      fread(&(table_nnB_params_0_s33[typ_ind][1][0] ) , binnum, 1, outfile);
      fread(&(table_nnB_params_0_s33[typ_ind][2][0] ) , binnum, 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s33[typ_ind][0][1]), &(table_nnB_params_0_s33[typ_ind][1][1]), &(table_nnB_params_0_s33[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fread(&(table_nnB_params_0_s33[typ_ind][0][1] ) , binnum, 1, outfile);
      fread(&(table_nnB_params_0_s33[typ_ind][1][1] ) , binnum, 1, outfile);
      fread(&(table_nnB_params_0_s33[typ_ind][2][1] ) , binnum, 1, outfile);
      ntot=table_nnB_N_0_s33[typ_ind][0]*table_nnB_N_0_s33[typ_ind][1]*table_nnB_N_0_s33[typ_ind][2];
      table_nnB_0_s33[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	//fscanf(table, "%lf", &(table_nnB_0_s33[typ_ind][ener]));
	fread(&(table_nnB_0_s33[typ_ind][ener] ) , binnum, 1, outfile);
      }
      //fclose(table);
      //}
    }
  }
  //S55
  table_nnB_0_s55=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //sprintf(tablename, "tab_energs/table_nnB_%d%d_0_s55.tab", i,j);
      //if((table=fopen(tablename, "r"))==NULL){
      //printf("No stacking interaction between %d and %d.\n", i,j);
      //    table_nnB_N_0_s55[typ_ind][0]=-1;
      //    table_nnB_N_0_s55[typ_ind][1]=-1;
      //    table_nnB_N_0_s55[typ_ind][2]=-1;
      //  }
      //else{
      //printf("[%d %d s55] ", i,j);
      //fscanf(table,"%d%d%d", &(table_nnB_N_0_s55[typ_ind][0]), &(table_nnB_N_0_s55[typ_ind][1]), &(table_nnB_N_0_s55[typ_ind][2])); //NX, NY, NZ
      fread(&(table_nnB_N_0_s55[typ_ind][0] ) , sizeof(int), 1, outfile);
      fread(&(table_nnB_N_0_s55[typ_ind][1] ) , sizeof(int), 1, outfile);
      fread(&(table_nnB_N_0_s55[typ_ind][2] ) , sizeof(int), 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s55[typ_ind][0][0]), &(table_nnB_params_0_s55[typ_ind][1][0]), &(table_nnB_params_0_s55[typ_ind][2][0]));//DX, DY, DZ
      fread(&(table_nnB_params_0_s55[typ_ind][0][0] ) , binnum, 1, outfile);
      fread(&(table_nnB_params_0_s55[typ_ind][1][0] ) , binnum, 1, outfile);
      fread(&(table_nnB_params_0_s55[typ_ind][2][0] ) , binnum, 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnB_params_0_s55[typ_ind][0][1]), &(table_nnB_params_0_s55[typ_ind][1][1]), &(table_nnB_params_0_s55[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fread(&(table_nnB_params_0_s55[typ_ind][0][1] ) , binnum, 1, outfile);
      fread(&(table_nnB_params_0_s55[typ_ind][1][1] ) , binnum, 1, outfile);
      fread(&(table_nnB_params_0_s55[typ_ind][2][1] ) , binnum, 1, outfile);
      ntot=table_nnB_N_0_s55[typ_ind][0]*table_nnB_N_0_s55[typ_ind][1]*table_nnB_N_0_s55[typ_ind][2];
      table_nnB_0_s55[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	//fscanf(table, "%lf", &(table_nnB_0_s55[typ_ind][ener]));
	fread(&(table_nnB_0_s55[typ_ind][ener] ) , binnum, 1, outfile);
      }
      //fclose(table);
      //}
    }
  }
  //////////////////////////////////
  //printf("\nReading STACKING-DIHEDRAL interactions...\n");
  table_nnB_1_s33=(double **)malloc(sizeof(double *)*ntypsq);
  /* for(i=0;i<N_BASES;i++){ */
  /*   for(j=0;j<N_BASES;j++){ */
  /*     typ_ind=N_BASES*i+j; */
  /*     //typ_ind2=N_BASES*j+i; */
  /*     //sprintf(tablename, "tab_energs/table_nnB_%d%d_1_s33.tab", i,j); */
  /*     //FILE *table; */
  /*     //if((table=fopen(tablename, "r"))==NULL){ */
  /*     //printf("[NO %d %d 33] ", i,j); */
  /*     //table_nnB_N_1_s33[typ_ind]=-1; */
  /*     //table_nnB_N_1[typ_ind2]=-1; */
  /*     //} */
  /*     //else{ */
  /*     //printf("[%d %d 33] ", i,j); */
  /*     //fscanf(table,"%d", &(table_nnB_N_1_s33[typ_ind])); //N_ETA */
  /*     fread(&(table_nnB_N_1_s33[typ_ind] ) , sizeof(int), 1, outfile); */
  /*     //fscanf(table,"%lf", &(table_nnB_params_1_s33[typ_ind][0]));//D_ETA */
  /*     fread(&(table_nnB_params_1_s33[typ_ind][0] ) , binnum, 1, outfile); */
  /*     //fscanf(table,"%lf", &(table_nnB_params_1_s33[typ_ind][1])); //ETAMIN */
  /*     fread(&(table_nnB_params_1_s33[typ_ind][1] ) , binnum, 1, outfile); */
  /*     ntot=table_nnB_N_1_s33[typ_ind]; */
  /*     table_nnB_1_s33[typ_ind]=(double *)malloc(sizeof(double)*ntot); */
  /*     for(ener=0;ener<ntot;ener++){ */
  /* 	//fscanf(table, "%lf", &(table_nnB_1_s33[typ_ind][ener])); */
  /* 	fread(&(table_nnB_1_s33[typ_ind][ener] ) , binnum, 1, outfile); */
  /*     } */
  /*     //fclose(table); */
  /*     //} */
  /*   } */
  /* } */
      
  /* table_nnB_1_s35=(double **)malloc(sizeof(double *)*ntypsq); */
  /* for(i=0;i<N_BASES;i++){ */
  /*   for(j=0;j<N_BASES;j++){ */
  /*     typ_ind=N_BASES*i+j; */
  /*     //typ_ind2=N_BASES*j+i; */
  /*     //sprintf(tablename, "tab_energs/table_nnB_%d%d_1_s35.tab", i,j); */
  /*     //FILE *table; */
  /*     //if((table=fopen(tablename, "r"))==NULL){ */
  /*     //printf("[NO %d %d 35] ", i,j); */
  /*     //table_nnB_N_1_s35[typ_ind]=-1; */
  /*     //table_nnB_N_1[typ_ind2]=-1; */
  /*     //} */
  /*     // else{ */
  /*     //printf("[%d %d 35] ", i,j); */
  /*     //fscanf(table,"%d", &(table_nnB_N_1_s35[typ_ind])); //N_ETA */
  /*     fread(&(table_nnB_N_1_s35[typ_ind] ) , sizeof(int), 1, outfile); */
  /*     //fscanf(table,"%lf", &(table_nnB_params_1_s35[typ_ind][0]));//D_ETA */
  /*     fread(&(table_nnB_params_1_s35[typ_ind][0] ) , binnum, 1, outfile); */
  /*     //fscanf(table,"%lf", &(table_nnB_params_1_s35[typ_ind][1])); //ETAMIN */
  /*     fread(&(table_nnB_params_1_s35[typ_ind][1] ) , binnum, 1, outfile); */
  /*     ntot=table_nnB_N_1_s35[typ_ind]; */
  /*     table_nnB_1_s35[typ_ind]=(double *)malloc(sizeof(double)*ntot); */
  /*     for(ener=0;ener<ntot;ener++){ */
  /* 	//fscanf(table, "%lf", &(table_nnB_1_s35[typ_ind][ener])); */
  /* 	fread(&(table_nnB_1_s35[typ_ind][ener] ) , binnum, 1, outfile); */
  /*     } */
  /*     //fclose(table); */
  /*     //} */
  /*   } */
  /* } */
  
  /* table_nnB_1_s53=(double **)malloc(sizeof(double *)*ntypsq); */
  /* for(i=0;i<N_BASES;i++){ */
  /*   for(j=0;j<N_BASES;j++){ */
  /*     typ_ind=N_BASES*i+j; */
  /*     //typ_ind2=N_BASES*j+i; */
  /*     //sprintf(tablename, "tab_energs/table_nnB_%d%d_1_s53.tab", i,j); */
  /*     //FILE *table; */
  /*     //if((table=fopen(tablename, "r"))==NULL){ */
  /*     //printf("[NO %d %d 53] ", i,j); */
  /*     //table_nnB_N_1_s53[typ_ind]=-1; */
  /*     //table_nnB_N_1[typ_ind2]=-1; */
  /*     //} */
  /*     // else{ */
  /*     //printf("[%d %d 53] ", i,j); */
  /*     //fscanf(table,"%d", &(table_nnB_N_1_s53[typ_ind])); //N_ETA */
  /*     fread(&(table_nnB_N_1_s53[typ_ind] ) , sizeof(int), 1, outfile); */
  /*     //fscanf(table,"%lf", &(table_nnB_params_1_s53[typ_ind][0]));//D_ETA */
  /*     fread(&(table_nnB_params_1_s53[typ_ind][0] ) , binnum, 1, outfile); */
  /*     //fscanf(table,"%lf", &(table_nnB_params_1_s53[typ_ind][1])); //ETAMIN */
  /*     fread(&(table_nnB_params_1_s53[typ_ind][1] ) , binnum, 1, outfile); */
  /*     ntot=table_nnB_N_1_s53[typ_ind]; */
  /*     table_nnB_1_s53[typ_ind]=(double *)malloc(sizeof(double)*ntot); */
  /*     for(ener=0;ener<ntot;ener++){ */
  /* 	//fscanf(table, "%lf", &(table_nnB_1_s53[typ_ind][ener])); */
  /* 	fread(&(table_nnB_1_s53[typ_ind][ener] ) , binnum, 1, outfile); */
  /*     } */
  /*     //fclose(table); */
  /*     //} */
  /*   } */
  /* } */
  
  /* table_nnB_1_s55=(double **)malloc(sizeof(double *)*ntypsq); */
  /* for(i=0;i<N_BASES;i++){ */
  /*   for(j=0;j<N_BASES;j++){ */
  /*     typ_ind=N_BASES*i+j; */
  /*     //typ_ind2=N_BASES*j+i; */
  /*     //sprintf(tablename, "tab_energs/table_nnB_%d%d_1_s55.tab", i,j); */
  /*     //FILE *table; */
  /*     //if((table=fopen(tablename, "r"))==NULL){ */
  /*     //printf("[NO %d %d 55] ", i,j); */
  /*     //    table_nnB_N_1_s55[typ_ind]=-1; */
  /*     //table_nnB_N_1[typ_ind2]=-1; */
  /*     //} */
  /*     //else{ */
  /*     //printf("[%d %d 55] ", i,j); */
  /*     //fscanf(table,"%d", &(table_nnB_N_1_s55[typ_ind])); //N_ETA */
  /*     fread(&(table_nnB_N_1_s55[typ_ind] ) , sizeof(int), 1, outfile); */
  /*     //fscanf(table,"%lf", &(table_nnB_params_1_s55[typ_ind][0]));//D_ETA */
  /*     fread(&(table_nnB_params_1_s55[typ_ind][0] ) , binnum, 1, outfile); */
  /*     //fscanf(table,"%lf", &(table_nnB_params_1_s55[typ_ind][1])); //ETAMIN */
  /*     fread(&(table_nnB_params_1_s55[typ_ind][1] ) , binnum, 1, outfile); */
  /*     ntot=table_nnB_N_1_s55[typ_ind]; */
  /*     table_nnB_1_s55[typ_ind]=(double *)malloc(sizeof(double)*ntot); */
  /*     for(ener=0;ener<ntot;ener++){ */
  /* 	//fscanf(table, "%lf", &(table_nnB_1_s55[typ_ind][ener])); */
  /* 	fread(&(table_nnB_1_s55[typ_ind][ener] ) , binnum, 1, outfile); */
  /*     } */
  /*     //fclose(table); */
  /*     //} */
  /*   } */
  /* } */
      
  //printf("\nReading NON-BONDED STACKING interactions...\n");
  /* STACKING -  NON-BONDED */
  //int nbtyp=1;
  table_nnN_0s3=(double **)malloc(sizeof(double *)*ntypsq);
  //typ_ind=0;
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      /*     //typ_ind2=N_BASES*j+i; */
      //sprintf(tablename, "tab_energs/table_nnN_%d%d_0s3.tab", i, j);
      //FILE *table;
      //if((table=fopen(tablename, "r"))==NULL){
      //printf("No generic non-bonded stacking interaction.\n");
      //table_nnN_N_0s3[typ_ind][0]=-1;
      //table_nnN_N_0s3[typ_ind][1]=-1;
      //table_nnN_N_0s3[typ_ind][2]=-1;
      //}
      // else{
      //printf("[%d %d s3] ", i,j);
      //fscanf(table,"%d%d%d", &(table_nnN_N_0s3[typ_ind][0]), &(table_nnN_N_0s3[typ_ind][1]), &(table_nnN_N_0s3[typ_ind][2])); //NX, NY, NZ
      fread(&(table_nnN_N_0s3[typ_ind][0] ) , sizeof(int), 1, outfile);
      fread(&(table_nnN_N_0s3[typ_ind][1] ) , sizeof(int), 1, outfile);
      fread(&(table_nnN_N_0s3[typ_ind][2] ) , sizeof(int), 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnN_params_0s3[typ_ind][0][0]), &(table_nnN_params_0s3[typ_ind][1][0]), &(table_nnN_params_0s3[typ_ind][2][0]));//DX, DY, DZ
      fread(&(table_nnN_params_0s3[typ_ind][0][0] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_0s3[typ_ind][1][0] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_0s3[typ_ind][2][0] ) , binnum, 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnN_params_0s3[typ_ind][0][1]), &(table_nnN_params_0s3[typ_ind][1][1]), &(table_nnN_params_0s3[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fread(&(table_nnN_params_0s3[typ_ind][0][1] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_0s3[typ_ind][1][1] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_0s3[typ_ind][2][1] ) , binnum, 1, outfile);
      ntot=table_nnN_N_0s3[typ_ind][0]*table_nnN_N_0s3[typ_ind][1]*table_nnN_N_0s3[typ_ind][2];
      table_nnN_0s3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	//fscanf(table, "%lf", &(table_nnN_0s3[typ_ind][ener]));
	fread(&(table_nnN_0s3[typ_ind][ener] ) , binnum, 1, outfile);
      }
      //fclose(table);
      //}
    }
  }
  
  table_nnN_0s5=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //typ_ind2=N_BASES*j+i;
      //typ_ind=0;
      //sprintf(tablename, "tab_energs/table_nnN_%d%d_0s5.tab", i, j);
      //if((table=fopen(tablename, "r"))==NULL){
      //printf("No generic non-bonded stacking_inv interaction.\n");
      //table_nnN_N_0s5[typ_ind][0]=-1;
      //table_nnN_N_0s5[typ_ind][1]=-1;
      //table_nnN_N_0s5[typ_ind][2]=-1;
      //}
      //else{
      //printf("[%d %d s5] ", i,j);
      //fscanf(table,"%d%d%d", &(table_nnN_N_0s5[typ_ind][0]), &(table_nnN_N_0s5[typ_ind][1]), &(table_nnN_N_0s5[typ_ind][2])); //NX, NY, NZ
      fread(&(table_nnN_N_0s5[typ_ind][0] ) , sizeof(int), 1, outfile);
      fread(&(table_nnN_N_0s5[typ_ind][1] ) , sizeof(int), 1, outfile);
      fread(&(table_nnN_N_0s5[typ_ind][2] ) , sizeof(int), 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnN_params_0s5[typ_ind][0][0]), &(table_nnN_params_0s5[typ_ind][1][0]), &(table_nnN_params_0s5[typ_ind][2][0]));//DX, DY, DZ
      fread(&(table_nnN_params_0s5[typ_ind][0][0] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_0s5[typ_ind][1][0] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_0s5[typ_ind][2][0] ) , binnum, 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnN_params_0s5[typ_ind][0][1]), &(table_nnN_params_0s5[typ_ind][1][1]), &(table_nnN_params_0s5[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fread(&(table_nnN_params_0s5[typ_ind][0][1] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_0s5[typ_ind][1][1] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_0s5[typ_ind][2][1] ) , binnum, 1, outfile);
      ntot=table_nnN_N_0s5[typ_ind][0]*table_nnN_N_0s5[typ_ind][1]*table_nnN_N_0s5[typ_ind][2];
      table_nnN_0s5[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	//fscanf(table, "%lf", &(table_nnN_0s5[typ_ind][ener]));
	fread(&(table_nnN_0s5[typ_ind][ener] ) , binnum, 1, outfile);
      }
      //fclose(table);
      //}
    }
  }
  
  /* WATSON-CRICK */
  //printf("\nReading BASE-PAIR interactions...\n");
  table_nnN_2=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //typ_ind2=N_BASES*j+i;
      //sprintf(tablename, "tab_energs/table_nnN_%d%d_2.tab", i,j);
      //FILE *table;
      //if((table=fopen(tablename, "r"))==NULL){
      //printf("[NO %d %d] ", i,j);
      //table_nnN_N_2[typ_ind][0]=-1;
      //table_nnN_N_2[typ_ind][1]=-1;
      //table_nnN_N_2[typ_ind][2]=-1;
      //table_nnN_N_2[typ_ind][3]=-1;
      //table_nnN_N_2[typ_ind2][0]=-1;
      //table_nnN_N_2[typ_ind2][1]=-1;
      //table_nnN_N_2[typ_ind2][2]=-1;
      //table_nnN_2[typ_ind]=(double *)malloc(sizeof(double)*10);
      //}
      //else{
      //printf("[%d %d] ", i,j);
      //fscanf(table,"%d%d%d", &(table_nnN_N_2[typ_ind][0]), &(table_nnN_N_2[typ_ind][1]), &(table_nnN_N_2[typ_ind][2])); //NX, NY, NZ
      fread(&(table_nnN_N_2[typ_ind][0] ) , sizeof(int), 1, outfile);
      fread(&(table_nnN_N_2[typ_ind][1] ) , sizeof(int), 1, outfile);
      fread(&(table_nnN_N_2[typ_ind][2] ) , sizeof(int), 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnN_params_2[typ_ind][0][0]), &(table_nnN_params_2[typ_ind][1][0]), &(table_nnN_params_2[typ_ind][2][0]));//DX, DY, DZ
      fread(&(table_nnN_params_2[typ_ind][0][0] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_2[typ_ind][1][0] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_2[typ_ind][2][0] ) , binnum, 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnN_params_2[typ_ind][0][1]), &(table_nnN_params_2[typ_ind][1][1]), &(table_nnN_params_2[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fread(&(table_nnN_params_2[typ_ind][0][1] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_2[typ_ind][1][1] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_2[typ_ind][2][1] ) , binnum, 1, outfile);
      ntot=table_nnN_N_2[typ_ind][0]*table_nnN_N_2[typ_ind][1]*table_nnN_N_2[typ_ind][2];
      table_nnN_N_2[typ_ind][3]=ntot;
      table_nnN_2[typ_ind]=(double *)malloc(sizeof(double)*ntot*WC_FACES);
      for(ener=0;ener<ntot*WC_FACES;ener++){
	//fscanf(table, "%lf", &(table_nnN_2[typ_ind][ener]));
	fread(&(table_nnN_2[typ_ind][ener] ) , binnum, 1, outfile);
	//table_nnN_2[typ_ind2][ener]=table_nnN_2[typ_ind][ener];
      }
      //fclose(table);
    }
  }
  
  table_nnN_2_inv=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //typ_ind2=N_BASES*j+i;
      //sprintf(tablename, "tab_energs/table_nnN_%d%d_2i.tab", i,j);
      //FILE *table;
      //if((table=fopen(tablename, "r"))==NULL){
      //printf("[NO inv %d %d] ", i,j);
      //table_nnN_N_2_inv[typ_ind][0]=-1;
      //table_nnN_N_2_inv[typ_ind][1]=-1;
      //table_nnN_N_2_inv[typ_ind][2]=-1;
      //table_nnN_N_2_inv[typ_ind][3]=-1;
      //table_nnN_N_2_inv[typ_ind2][0]=-1;
      //table_nnN_N_2_inv[typ_ind2][1]=-1;
      //table_nnN_N_2_inv[typ_ind2][2]=-1;
      //}
      //else{
      //printf("[inv %d %d] ", i,j);
      //fscanf(table,"%d%d%d", &(table_nnN_N_2_inv[typ_ind][0]), &(table_nnN_N_2_inv[typ_ind][1]), &(table_nnN_N_2_inv[typ_ind][2])); //NX, NY, NZ
      fread(&(table_nnN_N_2_inv[typ_ind][0] ) , sizeof(int), 1, outfile);
      fread(&(table_nnN_N_2_inv[typ_ind][1] ) , sizeof(int), 1, outfile);
      fread(&(table_nnN_N_2_inv[typ_ind][2] ) , sizeof(int), 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnN_params_2_inv[typ_ind][0][0]), &(table_nnN_params_2_inv[typ_ind][1][0]), &(table_nnN_params_2_inv[typ_ind][2][0]));//DX, DY, DZ
      fread(&(table_nnN_params_2_inv[typ_ind][0][0] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_2_inv[typ_ind][1][0] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_2_inv[typ_ind][2][0] ) , binnum, 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnN_params_2_inv[typ_ind][0][1]), &(table_nnN_params_2_inv[typ_ind][1][1]), &(table_nnN_params_2_inv[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fread(&(table_nnN_params_2_inv[typ_ind][0][1] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_2_inv[typ_ind][1][1] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_2_inv[typ_ind][2][1] ) , binnum, 1, outfile);
      ntot=table_nnN_N_2_inv[typ_ind][0]*table_nnN_N_2_inv[typ_ind][1]*table_nnN_N_2_inv[typ_ind][2];
      table_nnN_N_2_inv[typ_ind][3]=ntot;
      table_nnN_2_inv[typ_ind]=(double *)malloc(sizeof(double)*ntot*WC_FACES);
      for(ener=0;ener<ntot*WC_FACES;ener++){
	//fscanf(table, "%lf", &(table_nnN_2_inv[typ_ind][ener]));
	//table_nnN_2_inv[typ_ind2][ener]=table_nnN_2_inv[typ_ind][ener];
	fread(&(table_nnN_2_inv[typ_ind][ener] ) , binnum, 1, outfile);
      }
      //fclose(table);
      //}
    }
  }
  
  //printf("n max types : %d\n", N_MAX_TYPES);
  table_nnN_3=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //typ_ind2=N_BASES*j+i;
      //sprintf(tablename, "tab_energs/table_nnN_%d%d_3.tab", i,j);
      //FILE *table;
      //if((table=fopen(tablename, "r"))==NULL){
      //printf("[DIH %d %d] ", i,j);
      //table_nnN_N_3[typ_ind]=-1;
      //table_nnN_N_3[typ_ind2]=-1;
      //}
      //else{
      //printf("[DIH %d %d] ", i,j);
      //fscanf(table,"%d", &(table_nnN_N_3[typ_ind])); //N_THETA
      fread(&(table_nnN_N_3[typ_ind] ) , sizeof(int), 1, outfile);
      //fscanf(table,"%lf", &(table_nnN_params_3[typ_ind][0]));//D_THETA
      fread(&(table_nnN_params_3[typ_ind][0] ) , binnum, 1, outfile);
      //fscanf(table,"%lf", &(table_nnN_params_3[typ_ind][1])); //THETAMIN
      fread(&(table_nnN_params_3[typ_ind][1] ) , binnum, 1, outfile);
      ntot=table_nnN_N_3[typ_ind];
      table_nnN_3[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      //table_nnN_3[typ_ind2]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	//fscanf(table, "%lf", &(table_nnN_3[typ_ind][ener]));
	fread(&(table_nnN_3[typ_ind][ener] ) , binnum, 1, outfile);
	//table_nnN_3[typ_ind2][ener]=table_nnN_3[typ_ind][ener];
      }
      //fclose(table);
      //}
    }
  }
  
  /* for the case of parallel orientations */
  
  table_nnN_2_F=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //typ_ind2=N_BASES*j+i;
      //sprintf(tablename, "tab_energs/table_nnN_%d%d_2_F.tab", i,j);
      //FILE *table;
      //if((table=fopen(tablename, "r"))==NULL){
      //printf("[NO F %d %d] ", i,j);
      //table_nnN_N_2_F[typ_ind][0]=-1;
      //table_nnN_N_2_F[typ_ind][1]=-1;
      //table_nnN_N_2_F[typ_ind][2]=-1;
      //table_nnN_N_2_F[typ_ind][3]=-1;
      //}
      //else{
      //printf("[F %d %d] ", i,j);
      //fscanf(table,"%d%d%d", &(table_nnN_N_2_F[typ_ind][0]), &(table_nnN_N_2_F[typ_ind][1]), &(table_nnN_N_2_F[typ_ind][2])); //NX, NY, NZ
      fread(&(table_nnN_N_2_F[typ_ind][0] ) , sizeof(int), 1, outfile);
      fread(&(table_nnN_N_2_F[typ_ind][1] ) , sizeof(int), 1, outfile);
      fread(&(table_nnN_N_2_F[typ_ind][2] ) , sizeof(int), 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnN_params_2_F[typ_ind][0][0]), &(table_nnN_params_2_F[typ_ind][1][0]), &(table_nnN_params_2_F[typ_ind][2][0]));//DX, DY, DZ
      fread(&(table_nnN_params_2_F[typ_ind][0][0] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_2_F[typ_ind][1][0] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_2_F[typ_ind][2][0] ) , binnum, 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnN_params_2_F[typ_ind][0][1]), &(table_nnN_params_2_F[typ_ind][1][1]), &(table_nnN_params_2_F[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fread(&(table_nnN_params_2_F[typ_ind][0][1] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_2_F[typ_ind][1][1] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_2_F[typ_ind][2][1] ) , binnum, 1, outfile);
      ntot=table_nnN_N_2_F[typ_ind][0]*table_nnN_N_2_F[typ_ind][1]*table_nnN_N_2_F[typ_ind][2];
      table_nnN_N_2_F[typ_ind][3]=ntot;
      table_nnN_2_F[typ_ind]=(double *)malloc(sizeof(double)*ntot*WC_FACES);
      for(ener=0;ener<ntot*WC_FACES;ener++){
	//fscanf(table, "%lf", &(table_nnN_2_F[typ_ind][ener]));
	fread(&(table_nnN_2_F[typ_ind][ener] ) , binnum, 1, outfile);
      }
      //fclose(table);
      //}
    }
  }
  
  table_nnN_2_inv_F=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //typ_ind2=N_BASES*j+i;
      //sprintf(tablename, "tab_energs/table_nnN_%d%d_2i_F.tab", i,j);
      //FILE *table;
      //if((table=fopen(tablename, "r"))==NULL){
      //printf("[NO F inv %d %d] ", i,j);
      //table_nnN_N_2_inv_F[typ_ind][0]=-1;
      //table_nnN_N_2_inv_F[typ_ind][1]=-1;
      //table_nnN_N_2_inv_F[typ_ind][2]=-1;
      //table_nnN_N_2_inv_F[typ_ind][3]=-1;
      //}
      //else{
      //printf("[F inv %d %d] ", i,j);
      //fscanf(table,"%d%d%d", &(table_nnN_N_2_inv_F[typ_ind][0]), &(table_nnN_N_2_inv_F[typ_ind][1]), &(table_nnN_N_2_inv_F[typ_ind][2])); //NX, NY, NZ
      fread(&(table_nnN_N_2_inv_F[typ_ind][0] ) , sizeof(int), 1, outfile);
      fread(&(table_nnN_N_2_inv_F[typ_ind][1] ) , sizeof(int), 1, outfile);
      fread(&(table_nnN_N_2_inv_F[typ_ind][2] ) , sizeof(int), 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnN_params_2_inv_F[typ_ind][0][0]), &(table_nnN_params_2_inv_F[typ_ind][1][0]), &(table_nnN_params_2_inv_F[typ_ind][2][0]));//DX, DY, DZ
      fread(&(table_nnN_params_2_inv_F[typ_ind][0][0] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_2_inv_F[typ_ind][1][0] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_2_inv_F[typ_ind][2][0] ) , binnum, 1, outfile);
      //fscanf(table,"%lf%lf%lf", &(table_nnN_params_2_inv_F[typ_ind][0][1]), &(table_nnN_params_2_inv_F[typ_ind][1][1]), &(table_nnN_params_2_inv_F[typ_ind][2][1])); //XMIN, YMIN, ZMIN
      fread(&(table_nnN_params_2_inv_F[typ_ind][0][1] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_2_inv_F[typ_ind][1][1] ) , binnum, 1, outfile);
      fread(&(table_nnN_params_2_inv_F[typ_ind][2][1] ) , binnum, 1, outfile);
      ntot=table_nnN_N_2_inv_F[typ_ind][0]*table_nnN_N_2_inv_F[typ_ind][1]*table_nnN_N_2_inv_F[typ_ind][2];
      table_nnN_N_2_inv_F[typ_ind][3]=ntot;
      table_nnN_2_inv_F[typ_ind]=(double *)malloc(sizeof(double)*ntot*WC_FACES);
      for(ener=0;ener<ntot*WC_FACES;ener++){
	//fscanf(table, "%lf", &(table_nnN_2_inv_F[typ_ind][ener]));
	fread(&(table_nnN_2_inv_F[typ_ind][ener] ) , binnum, 1, outfile);
      }
      //fclose(table);
      //}
    }
  }
  
  table_nnN_3_F=(double **)malloc(sizeof(double *)*ntypsq);
  for(i=0;i<N_BASES;i++){
    for(j=0;j<N_BASES;j++){
      typ_ind=N_BASES*i+j;
      //sprintf(tablename, "tab_energs/table_nnN_%d%d_3_F.tab", i,j);
      //if((table=fopen(tablename, "r"))==NULL){
      //printf("[NO F DIH %d %d] ", i,j);
      //table_nnN_N_3_F[typ_ind]=-1;
      //}
      //else{
      //printf("[F DIH %d %d] ", i,j);
      //fscanf(table,"%d", &(table_nnN_N_3_F[typ_ind])); //N_THETA
      fread(&(table_nnN_N_3_F[typ_ind] ) , sizeof(int), 1, outfile);
      //fscanf(table,"%lf", &(table_nnN_params_3_F[typ_ind][0]));//D_THETA
      fread(&(table_nnN_params_3_F[typ_ind][0] ) , binnum, 1, outfile);
      //fscanf(table,"%lf", &(table_nnN_params_3_F[typ_ind][1])); //THETAMIN
      fread(&(table_nnN_params_3_F[typ_ind][1] ) , binnum, 1, outfile);
      ntot=table_nnN_N_3_F[typ_ind];
      table_nnN_3_F[typ_ind]=(double *)malloc(sizeof(double)*ntot);
      for(ener=0;ener<ntot;ener++){
	//fscanf(table, "%lf", &(table_nnN_3_F[typ_ind][ener]));
	fread(&(table_nnN_3_F[typ_ind][ener] ) , binnum, 1, outfile);
      }
      //fclose(table);
      //}
    }
  }
  /***********************/
  
  //printf("\nReading BASE-PHOSPHATE interactions...\n");
  // BASE-PHOSPHATE
  table_npN_0=(double **)malloc(sizeof(double *)*N_BASES);
  for(i=0;i<N_BASES;i++){
    //sprintf(tablename, "tab_energs/table_npN_%d_0.tab", i);
    //if((table=fopen(tablename, "r"))==NULL){
    //printf("[NO BPH %d] ", i);
    //table_npN_N_0[i][0]=-1;
    //table_npN_N_0[i][1]=-1;
    //table_npN_N_0[i][2]=-1;
    //}
    //else{
    //printf("[BPH %d] ", i);
    //fscanf(table,"%d%d%d", &(table_npN_N_0[i][0]), &(table_npN_N_0[i][1]), &(table_npN_N_0[i][2])); //NX, NY, NZ
    fread(&(table_npN_N_0[i][0] ) , sizeof(int), 1, outfile);
    fread(&(table_npN_N_0[i][1] ) , sizeof(int), 1, outfile);
    fread(&(table_npN_N_0[i][2] ) , sizeof(int), 1, outfile);
    //fscanf(table,"%lf%lf%lf", &(table_npN_params_0[i][0][0]), &(table_npN_params_0[i][1][0]), &(table_npN_params_0[i][2][0]));//DX, DY, DZ
    fread(&(table_npN_params_0[i][0][0] ) , binnum, 1, outfile);
    fread(&(table_npN_params_0[i][1][0] ) , binnum, 1, outfile);
    fread(&(table_npN_params_0[i][2][0] ) , binnum, 1, outfile);
    //fscanf(table,"%lf%lf%lf", &(table_npN_params_0[i][0][1]), &(table_npN_params_0[i][1][1]), &(table_npN_params_0[i][2][1])); //XMIN, YMIN, ZMIN
    fread(&(table_npN_params_0[i][0][1] ) , binnum, 1, outfile);
    fread(&(table_npN_params_0[i][1][1] ) , binnum, 1, outfile);
    fread(&(table_npN_params_0[i][2][1] ) , binnum, 1, outfile);
    ntot=table_npN_N_0[i][0]*table_npN_N_0[i][1]*table_npN_N_0[i][2];
    table_npN_0[i]=(double *)malloc(sizeof(double)*ntot);
    for(ener=0;ener<ntot;ener++){
      //fscanf(table, "%lf", &(table_npN_0[i][ener]));
      fread(&(table_npN_0[i][ener] ) , binnum, 1, outfile);
    }
    //fclose(table);
    //}
  }
  fclose(outfile);
  //printf("\n");
}


int main(){
  MC_read_write_energy_tables();
  //MC_read_bin_energy_tables();
  return 0;
}
