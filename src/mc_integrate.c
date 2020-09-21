#include "mc_integrate.h"

double MC_integrate(int mc_n, double **rx, double **ry, double **rz){
  int nt_n=mc_n/N_PARTS_PER_NT, i;
  int nt_c, trial;
  double old_energy, new_energy, wc_energ_diff=0.0, dwce=0; // wc energ is the difference between the new and the old : Ewc - Ewc_old
  int mc_trial_flag=0;
  int nt_to_select=nt_n;
  double d_energ=0.0;
  double oUCM[DIM],nUCM[DIM], oeUCM, neUCM;
  /** RG STUFF**/
  double oRG2=0,nRG2=0,oeRG,neRG,rgnats=0;
  double RCM[DIM];
  int rgat,d,rgnt;
#ifdef MCVLISTS
  for(i=0;i<nt_n;i++){
    if(vl_count[i]>=vl_ncrit-1)
      MC_build_verlet_lists(nt_n, *rx, *ry, *rz, i);
  }
#endif
  for(nt_c=0;nt_c<nt_n;nt_c++){
    mc_trial_flag=0;
    if(fr_is_mobile[nt_c]!=FR_MOB_FROZ){
/*       #ifdef ERMSDR */
/*       double tempermsd = get_first_ermsd(rx, ry, rz, nt_n, &ERMSD_SQ, &ERMSD_ENERG); */
/* #endif */
      MC_copy_nt(nt_c, *rx, *ry, *rz);

      /** CALC RADIUS OF GYRATION **/
      if(KRG>0 || UMBRELLA_TYPE==0){
	for(d=0;d<DIM;d++) RCM[d]=0;oRG2=0;
	rgnats=0;
	//for(rgat=0;rgat<mc_n;rgat++){
	for(rgnt=0;rgnt<nt_n;rgnt++){
	  if(fr_is_mobile[rgnt]!=FR_MOB_FROZ){
	    for(rgat=0;rgat<N_PARTS_PER_NT;rgat++){
	      RCM[0]+=(*rx)[rgnt*N_PARTS_PER_NT+rgat];
	      RCM[1]+=(*ry)[rgnt*N_PARTS_PER_NT+rgat];
	      RCM[2]+=(*rz)[rgnt*N_PARTS_PER_NT+rgat];
	      rgnats++;
	    }
	  }
	}
	for(d=0;d<DIM;d++) RCM[d]/=(double)rgnats;
	if(KRG>0){
	  for(rgnt=0;rgnt<nt_n;rgnt++){
	  if(fr_is_mobile[rgnt]!=FR_MOB_FROZ){
	    for(rgat=0;rgat<N_PARTS_PER_NT;rgat++){
	      if(rgnt==nt_c){
		oRG2+=SQ(mc_temp_x[rgnt*N_PARTS_PER_NT+rgat]-RCM[0]);
		oRG2+=SQ(mc_temp_y[rgnt*N_PARTS_PER_NT+rgat]-RCM[1]);
		oRG2+=SQ(mc_temp_z[rgnt*N_PARTS_PER_NT+rgat]-RCM[2]);
	      }
	      else{
		oRG2+=SQ((*rx)[rgnt*N_PARTS_PER_NT+rgat]-RCM[0]);
		oRG2+=SQ((*ry)[rgnt*N_PARTS_PER_NT+rgat]-RCM[1]);
		oRG2+=SQ((*rz)[rgnt*N_PARTS_PER_NT+rgat]-RCM[2]);
	      }
	    }
	  }
	  }
	oRG2/=(double)rgnats;
	}//printf("%lf %lf\t %lf %lf %lf   %lf\n", sqrt(oRG2), RG_target,RCM[0],RCM[1],RCM[2], rgnats);
	if(UMBRELLA_TYPE==0){
	  for(d=0;d<DIM;d++)
	    oUCM[d]=RCM[d];
	}
      }
      /******************************/
      
      mc_trial_flag=MC_calculate_local_energy(*rx, *ry, *rz, nt_c, &old_energy, nt_n, -1); //here we calculate the energy of the selected particle, from temp position
      
      if(mc_trial_flag!=0){ printf("THIS SHOULDNT HAPPEN\n flag = %d, selected %d\n", mc_trial_flag, nt_c);
	for(i=0;i<nt_n;i++){
	  printf("%d  GLYC %d  PUCK %d   -    temp: GLYC %d  PUCK  %d\n", i, mc_glyc[i], mc_puck[i], mc_temp_glyc[i], mc_temp_puck[i]);
	}
	MC_get_sp_bonded_energy(nt_c, nt_n, *rx, *ry, *rz);
	exit(1);
      }
      
      trial=MC_perform_trial_movement(nt_c); //here we update the position to temp position
      
      mc_trial_flag=MC_calculate_local_energy(*rx, *ry, *rz, nt_c, &new_energy, nt_n, trial); //here we calculate the energy of the selected particle, from temp position

      /** CALC RADIUS OF GYRATION AGAIN**/
      if(KRG>0 || UMBRELLA_TYPE==0){
	for(d=0;d<DIM;d++) RCM[d]=0;nRG2=0;
	rgnats=0;
	for(rgnt=0;rgnt<nt_n;rgnt++){
	  if(fr_is_mobile[rgnt]!=FR_MOB_FROZ){
	    for(rgat=0;rgat<N_PARTS_PER_NT;rgat++){
	      if(rgnt==nt_c){
		RCM[0]+=mc_temp_x[rgnt*N_PARTS_PER_NT+rgat];
		RCM[1]+=mc_temp_y[rgnt*N_PARTS_PER_NT+rgat];
		RCM[2]+=mc_temp_z[rgnt*N_PARTS_PER_NT+rgat];
		rgnats++;
	      }
	      else{
		RCM[0]+=(*rx)[rgnt*N_PARTS_PER_NT+rgat];
		RCM[1]+=(*ry)[rgnt*N_PARTS_PER_NT+rgat];
		RCM[2]+=(*rz)[rgnt*N_PARTS_PER_NT+rgat];
		rgnats++;
	      }
	    }
	  }
	}
	for(d=0;d<DIM;d++) RCM[d]/=(double)rgnats;
	if(KRG>0){
	  for(rgnt=0;rgnt<nt_n;rgnt++){
	    if(fr_is_mobile[rgnt]!=FR_MOB_FROZ){
	      for(rgat=0;rgat<N_PARTS_PER_NT;rgat++){
		if(rgnt==nt_c){
		  nRG2+=SQ(mc_temp_x[rgnt*N_PARTS_PER_NT+rgat]-RCM[0]);
		  nRG2+=SQ(mc_temp_y[rgnt*N_PARTS_PER_NT+rgat]-RCM[1]);
		  nRG2+=SQ(mc_temp_z[rgnt*N_PARTS_PER_NT+rgat]-RCM[2]);
		}
		else{
		  nRG2+=SQ((*rx)[rgnt*N_PARTS_PER_NT+rgat]-RCM[0]);
		  nRG2+=SQ((*ry)[rgnt*N_PARTS_PER_NT+rgat]-RCM[1]);
		  nRG2+=SQ((*rz)[rgnt*N_PARTS_PER_NT+rgat]-RCM[2]);
		}
	      }
	    }
	  }
	  nRG2/=(double)rgnats;
	  /******************************/
	  /** add RG energy **/
	  oeRG=KRG*SQ(sqrt(oRG2)-RG_target); 
	  neRG=KRG*SQ(sqrt(nRG2)-RG_target);
	  //printf("SECOND %lf %lf  \t %lf %lf %lf   %lf\n", sqrt(nRG2), RG_target,RCM[0],RCM[1],RCM[2], rgnats);
	}
	if(UMBRELLA_TYPE==0){
	  for(d=0;d<DIM;d++)
	    nUCM[d]=RCM[d];
	  
	  oeUCM=0.5*UCMK0*SQ(UCM[0]-oUCM[0])+0.5*UCMK1*SQ(UCM[1]-oUCM[1])+0.5*UCMK2*SQ(UCM[2]-oUCM[2]);
	  neUCM=0.5*UCMK0*SQ(UCM[0]-nUCM[0])+0.5*UCMK1*SQ(UCM[1]-nUCM[1])+0.5*UCMK2*SQ(UCM[2]-nUCM[2]);
	}
      }
#ifndef NOCTCS
      dwce=MC_calculate_local_wc_energy(nt_n, nt_c, *rx, *ry, *rz);
#endif
      if(mc_trial_flag==0){
	//printf("%lf %lf    %lf %lf\n", neRG, oeRG, new_energy, old_energy);
	//printf("  %lf\n", new_energy);
	if(KRG>0) {new_energy+=neRG-oeRG; }
	if(UMBRELLA_TYPE==0){new_energy+=neUCM-oeUCM;}
	d_energ+=MC_eval_displacement(nt_n,rx, ry, rz, nt_c, old_energy, new_energy, dwce); //if the movement is accepted, we update the position
      }
    }
  }
  return d_energ;
}

void MC_print_radius_of_gyration(int nt_n, double *rx, double *ry, double *rz){
  int mc_n=nt_n*N_PARTS_PER_NT,i,d,rgnt,rgat;
  double RG2, RCM[DIM];
  for(d=0;d<DIM;d++) RCM[d]=0;RG2=0;
  for(rgnt=0;rgnt<nt_n;rgnt++){
    for(rgat=0;rgat<N_PARTS_PER_NT;rgat++){
      RCM[0]+=rx[rgnt*N_PARTS_PER_NT+rgat];
      RCM[1]+=ry[rgnt*N_PARTS_PER_NT+rgat];
      RCM[2]+=rz[rgnt*N_PARTS_PER_NT+rgat];
    }
  }
  for(d=0;d<DIM;d++) RCM[d]/=(double)mc_n;
  for(rgnt=0;rgnt<nt_n;rgnt++){
    for(rgat=0;rgat<N_PARTS_PER_NT;rgat++){
      RG2+=SQ(rx[rgnt*N_PARTS_PER_NT+rgat]-RCM[0]);
      RG2+=SQ(ry[rgnt*N_PARTS_PER_NT+rgat]-RCM[1]);
      RG2+=SQ(rz[rgnt*N_PARTS_PER_NT+rgat]-RCM[2]);
    }
  }
  RG2/=mc_n;
  printf("RG: %lf\n",sqrt(RG2));
}

void MC_print_secondary_structure(int nt_n, double *rx, double *ry, double *rz){
  int nt_c, at_c, typ1, typ2;
  double b_st=0.0, b_sg=0.0, b_ev=0.0,nb_st=0.0, nb_wc=0.0;
  double tot_b_st=0.0, tot_b_sg=0.0, tot_b_ev=0.0,tot_nb_st=0.0, tot_nb_wc=0.0, tot_nb_ev=0.0;
  double d_vec[DIM], r_vec[DIM], r_vec_inv[DIM];
  double dist, sugdistsq;
  int typ_ind, typ_ind2;
  //int flag=0;
  int b, nt_neigh, ba, at_ne;
  double eta, theta, eta2;
  int temp_flag;
  double epsilon, sigma;
  int nt2,n;
  int wc_pair_array1[3],wc_pair_array2[3];
  int *secstr;
  secstr=(int*)malloc(sizeof(int)*nt_n);
  for(nt_c=0;nt_c<nt_n;nt_c++)
    secstr[nt_c]=-1;
  
  int temp_face1, temp_face2;
  double energ_wc=MC_calculate_total_wc_energy(nt_n,-1,rx, ry, rz);
  MC_update_wc_lists(nt_n); 
  for(nt_c=0;nt_c<nt_n;nt_c++){
    at_c=N_PARTS_PER_NT*nt_c;
    MC_copy_nt(nt_c, rx, ry, rz);
    typ1=mc_types[nt_c*N_PARTS_PER_NT];
    printf("%d ", nt_c);
    if(typ1==TYP_ADENINE)
      printf("(A)  :");
    if(typ1==TYP_URACIL)
      printf("(U)  :");
    if(typ1==TYP_GUANINE)
      printf("(G)  :");
    if(typ1==TYP_CYTOSINE)
      printf("(C)  :");
    
    //Base pairs
    for(n=0;n<vl_n_pairs[nt_c];n++){
      nt2=vl_neighbor_tables[nt_c][n];
      at_ne=N_PARTS_PER_NT*nt2;
      typ2=mc_types[nt2*N_PARTS_PER_NT];
      MC_is_wc_corresp(nt_c, nt2, nt_n, wc_pair_array1);
      MC_is_wc_corresp(nt2, nt_c, nt_n, wc_pair_array2);
      temp_face1=wc_pair_array1[0];
      temp_face2=wc_pair_array2[0];
      //if(temp_face1>=0 && temp_face2>=0) {
      if(temp_face1==WC_FACE_WATSCRICK && temp_face2==WC_FACE_WATSCRICK && is_canonical(typ1, typ2) == 1) {
	theta=calc_wc_psdihedr(mc_temp_x,mc_temp_y,mc_temp_z,rx, ry, rz,nt_c,nt2);
	//if(theta>=0 && theta <= M_PI/2)
	//printf(" pWC %d ", nt2);
	if(theta>M_PI/2 && theta <= M_PI){
	  printf(" %d ", nt2);
	  if(typ2==TYP_ADENINE)
	    printf("(A)");
	  if(typ2==TYP_URACIL)
	    printf("(U)");
	  if(typ2==TYP_GUANINE)
	    printf("(G)");
	  if(typ2==TYP_CYTOSINE)
	    printf("(C)");
	  if(nt_c<nt2) secstr[nt_c]=0;
	  if(nt_c>nt2) secstr[nt_c]=1;
	  
	}
      }
    }
    printf("\n");
  }
  for(nt_c=0;nt_c<nt_n;nt_c++){
     typ1=mc_types[nt_c*N_PARTS_PER_NT];
      if(typ1==TYP_ADENINE)
	printf("A");
      if(typ1==TYP_URACIL)
	printf("U");
      if(typ1==TYP_GUANINE)
	printf("G");
      if(typ1==TYP_CYTOSINE)
	printf("C");
      
  }
  printf("\n");
  for(nt_c=0;nt_c<nt_n;nt_c++){
    typ1=secstr[nt_c];

    if(typ1==-1)
      printf(".");
    if(typ1==0)
      printf("(");
    if(typ1==1)
    printf(")");
  }
  printf("\n");
  free(secstr);
}

int is_canonical(int typ1, int typ2){
  int ret=-1;
  if(typ1==TYP_ADENINE)
    if(typ2==TYP_URACIL)
      ret=1;
  if(typ1==TYP_URACIL)
    if(typ2==TYP_ADENINE)
      ret=1;
  if(typ1==TYP_CYTOSINE)
    if(typ2==TYP_GUANINE)
      ret=1;
  if(typ1==TYP_GUANINE)
    if(typ2==TYP_CYTOSINE)
      ret=1;
  return ret;
}


void MC_print_contact_list(int nt_n, double *rx, double *ry, double *rz){
  int nt_c, at_c;
  double b_st=0.0, b_sg=0.0, b_ev=0.0,nb_st=0.0, nb_wc=0.0, nb_ev=0.0;
  double tot_b_st=0.0, tot_b_sg=0.0, tot_b_ev=0.0,tot_nb_st=0.0, tot_nb_wc=0.0, tot_nb_ev=0.0;
  double d_vec[DIM], r_vec[DIM], r_vec_inv[DIM];
  double dist, sugdistsq;
  int typ_ind, typ_ind2;
  //int flag=0;
  int b, nt_neigh, ba, at_ne;
  double eta, theta, eta2;
  int temp_flag;
  double epsilon, sigma;
  int nt2,n;
  int wc_pair_array1[3],wc_pair_array2[3];
  
  int temp_face1, temp_face2;
  double energ_wc=MC_calculate_total_wc_energy(nt_n,-1,rx, ry, rz);
  MC_update_wc_lists(nt_n); 
  for(nt_c=0;nt_c<nt_n;nt_c++){
    at_c=N_PARTS_PER_NT*nt_c;
    MC_copy_nt(nt_c, rx, ry, rz);
    printf("\n%d ", nt_c);
    if(mc_types[nt_c*N_PARTS_PER_NT]==0)
      printf("(A)  :");
    if(mc_types[nt_c*N_PARTS_PER_NT]==1)
      printf("(U)  :");
    if(mc_types[nt_c*N_PARTS_PER_NT]==2)
      printf("(G)  :");
    if(mc_types[nt_c*N_PARTS_PER_NT]==3)
      printf("(C)  :");
    
    for(b=0;b<mc_nbonds[nt_c][0];b++){
      nt_neigh=mc_bondlist[nt_c][b];
      //first, we see if they are stacked!
      ba=nt_neigh*N_PARTS_PER_NT;
      typ_ind = mc_n_types*mc_types[at_c]     + mc_types[ba];
      typ_ind2= mc_n_types*mc_types[ba] + mc_types[at_c];
      d_vec[0]=get_unf_coo_x(rx, ba) - get_unf_coo_temp_x(at_c);
      d_vec[1]=get_unf_coo_y(ry, ba) - get_unf_coo_temp_y(at_c);
      d_vec[2]=get_unf_coo_z(rz, ba) - get_unf_coo_temp_z(at_c);
      
      dist=sqrt(d_vec[0]*d_vec[0]+d_vec[1]*d_vec[1]+d_vec[2]*d_vec[2]);
      proj_on_nt(d_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
      proj_on_nt_inv(d_vec, rx, ry, rz, nt_neigh, r_vec_inv);
      //calc_rel_pos(d_vec, r_vec);
      //calc_rel_pos_inv(d_vec, nt_neigh,rx, ry,rz, r_vec_inv);
      
      eta=calc_st_psdihedr(mc_temp_x, mc_temp_y, mc_temp_z,rx,ry,rz,nt_c,nt_neigh); 
      temp_flag=0;
      eta2=calc_wc_psdihedr(mc_temp_x, mc_temp_y, mc_temp_z, rx,ry,rz,nt_c,nt_neigh); 
      b_st=0;
      if(fabs(eta2)<STANG || fabs(eta2) > M_PI-STANG)
	b_st=MC_calc_nnB_stacking(typ_ind, typ_ind2, r_vec, r_vec_inv, eta, &temp_flag);
      else temp_flag=1;
      //printf("flag : %d \n", *mc_flag);
      if(temp_flag==0){
	printf(" stb  %d  ", nt_neigh);// if(cflag==-1)printf(" %lf ", b_st);
      }
    }
    
    //NON BONDED
    //if(cflag!=1) printf("\nnonbonded : ");
    
    //for(n=0;n<vl_n_pairs[nt_c];n++){
    for(nt2=0;nt2<nt_n;nt2++){
      //nt2=vl_neighbor_tables[nt_c][n];
      at_ne=N_PARTS_PER_NT*nt2;
      typ_ind = mc_n_types*mc_types[at_c]  + mc_types[at_ne];
	typ_ind2= mc_n_types*mc_types[at_ne] + mc_types[at_c];
	sugdistsq=calc_min_dist_sq(rx[at_ne+ISUG], ry[at_ne+ISUG], rz[at_ne+ISUG], mc_temp_x[at_c+ISUG], mc_temp_y[at_c+ISUG], mc_temp_z[at_c+ISUG]);
	
	
	if(sugdistsq<mc_r_cut_sq){
	  calc_min_vec(rx[at_ne], ry[at_ne], rz[at_ne], mc_temp_x[at_c], mc_temp_y[at_c], mc_temp_z[at_c], d_vec, &dist);
	  proj_on_nt(d_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
	  proj_on_nt_inv(d_vec, rx, ry, rz, nt2, r_vec_inv);
	  
	  MC_is_wc_corresp(nt_c, nt2, nt_n, wc_pair_array1);
	  MC_is_wc_corresp(nt2, nt_c, nt_n, wc_pair_array2);
	  temp_face1=wc_pair_array1[0];
	  temp_face2=wc_pair_array2[0];
	  //printf("%d   %d %d %d\t%d %d %d\n", nt_c,wc_pair_array1[0], wc_pair_array1[1], wc_pair_array1[2], wc_pair_array2[0], wc_pair_array2[1], wc_pair_array2[2]);
	  if(temp_face1>=0 && temp_face2>=0) {
	    printf(" wcp %d (", nt2);
	    if(temp_face2==WC_FACE_SUGAR)
	      printf("S");
	    else if(temp_face2==WC_FACE_WATSCRICK)
	      printf("W");
	    else if(temp_face2==WC_FACE_HOOGSTEEN)
	      printf("H");
	    if(temp_face1==WC_FACE_SUGAR)
	      printf("S");
	    else if(temp_face1==WC_FACE_WATSCRICK)
	      printf("W");
	    else if(temp_face1==WC_FACE_HOOGSTEEN)
	      printf("H");
	    printf(")");
	  }
	  
	  temp_face1=wc_pair_array1[1];
	  temp_face2=wc_pair_array2[2];
	  if(temp_face1>=0 && temp_face2>=0) {
	    printf(" bph %d (", nt2);
	    if(temp_face2==WC_FACE_SUGAR)
	      printf("S");
	    else if(temp_face2==WC_FACE_WATSCRICK)
	      printf("W");
	    else if(temp_face2==WC_FACE_HOOGSTEEN)
	      printf("H");
	    printf("P)");
	  }
	  
	  temp_face1=wc_pair_array1[2];
	  temp_face2=wc_pair_array2[1];
	  if(temp_face1>=0 && temp_face2>=0) {
	    printf(" bph %d (P", nt2);
	    if(temp_face1==WC_FACE_SUGAR)
	    printf("S");
	  else if(temp_face1==WC_FACE_WATSCRICK)
	    printf("W");
	  else if(temp_face1==WC_FACE_HOOGSTEEN)
	    printf("H");
	  printf(")");
	  }
	
	  if(!MC_are_neighbors(nt_c,nt2)){
	    temp_flag=0;
	    eta2=calc_wc_psdihedr(mc_temp_x, mc_temp_y, mc_temp_z, rx,ry,rz,nt_c,nt2);
	    nb_st=0;
	    if(fabs(eta2)<STANG || fabs(eta2) > M_PI-STANG)
	      nb_st = MC_calc_nnN_stacking(typ_ind,typ_ind2, r_vec, r_vec_inv, &temp_flag); // NO PAIR SPECIFIC!!
	    else temp_flag=1;
	    
	    if(temp_flag==0){
	      if(nb_st<0)printf(" nst %d ", nt2);// if(cflag==-1)printf( " %lf ", nb_st);
	    }
	  }
	}
    }
        
    printf("\n");
    
  
  }
}


  
void MC_write_contact_list(int nt_n, double *rx, double *ry, double *rz, int init){
  FILE * mcconf;
  char filename[256];
  
  sprintf(filename, "configs/min_energ.p%02d.map", init);
  
  int nt_c, at_c;
  double b_st=0.0, b_sg=0.0, b_ev=0.0,nb_st=0.0, nb_wc=0.0, nb_ev=0.0;
  double tot_b_st=0.0, tot_b_sg=0.0, tot_b_ev=0.0,tot_nb_st=0.0, tot_nb_wc=0.0, tot_nb_ev=0.0;
  double d_vec[DIM], r_vec[DIM], r_vec_inv[DIM];
  double dist, sugdistsq;
  int typ_ind, typ_ind2;
  //int flag=0;
  int b, nt_neigh, ba;
  double eta, theta, eta2;
  int temp_flag;
  double epsilon, sigma;
  int nt2,n;
  int wc_pair_array1[3],wc_pair_array2[3];
  //we calculate wc energy to leave the wc arrays initialized
  double energ_wc=MC_calculate_total_wc_energy(nt_n,-1,rx, ry, rz);
  MC_update_wc_lists(nt_n);
  int temp_face1, temp_face2;
  
  if(!(mcconf=fopen(filename , "w"))){
    fprintf(stderr, "Unable to write minimum energy MC configuration file.\n");
    //exit(ERR_INPUT);
  } else {
    for(nt_c=0;nt_c<nt_n;nt_c++){
      at_c=N_PARTS_PER_NT*nt_c;
      MC_copy_nt(nt_c, rx, ry, rz);
      fprintf(mcconf, "\n%d ", nt_c);
      if(mc_types[nt_c*N_PARTS_PER_NT]==0)
	fprintf(mcconf, "(A)  :");
      if(mc_types[nt_c*N_PARTS_PER_NT]==1)
	fprintf(mcconf, "(U)  :");
      if(mc_types[nt_c*N_PARTS_PER_NT]==2)
	fprintf(mcconf, "(G)  :");
      if(mc_types[nt_c*N_PARTS_PER_NT]==3)
	fprintf(mcconf, "(C)  :");
      
      for(b=0;b<mc_nbonds[nt_c][0];b++){
	nt_neigh=mc_bondlist[nt_c][b];
	//first, we see if they are stacked!
	ba=nt_neigh*N_PARTS_PER_NT;
	typ_ind = mc_n_types*mc_types[at_c]     + mc_types[N_PARTS_PER_NT*nt_neigh];
	typ_ind2= mc_n_types*mc_types[N_PARTS_PER_NT*nt_neigh] + mc_types[at_c];
	d_vec[0]=get_unf_coo_x(rx, ba) - get_unf_coo_temp_x(at_c+0);
	//x[ba]+mc_pbox[ba][0]*box_l[0]) - (mc_temp_x[0]+mc_temp_pbox[0][0]*box_l[0]);
	d_vec[1]=get_unf_coo_y(ry, ba) - get_unf_coo_temp_y(at_c+0);
	d_vec[2]=get_unf_coo_z(rz, ba) - get_unf_coo_temp_z(at_c+0);
	
	dist=sqrt(d_vec[0]*d_vec[0]+d_vec[1]*d_vec[1]+d_vec[2]*d_vec[2]);
	proj_on_nt(d_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
	proj_on_nt_inv(d_vec, rx, ry, rz, nt_neigh, r_vec_inv);
	//calc_rel_pos(d_vec, r_vec);
	//calc_rel_pos_inv(d_vec, nt_neigh,rx, ry,rz, r_vec_inv);
	eta=calc_st_psdihedr(mc_temp_x, mc_temp_y, mc_temp_z,rx,ry,rz,nt_c,nt_neigh); 
	temp_flag=0;
	b_st=0;
	eta2=calc_wc_psdihedr(mc_temp_x, mc_temp_y, mc_temp_z, rx,ry,rz,nt_c,nt_neigh); 
	if(fabs(eta2)<STANG || fabs(eta2) > M_PI-STANG)
	  b_st=MC_calc_nnB_stacking(typ_ind, typ_ind2, r_vec, r_vec_inv, eta, &temp_flag);
	else
	  temp_flag=1;
	if(temp_flag==0){
	  fprintf(mcconf, " stb  %d  ", nt_neigh);// if(cflag==-1)printf(" %lf ", b_st);
	}
      }
      
      //NON BONDED
      //if(cflag!=1) printf("\nnonbonded : ");
      for(n=0;n<vl_n_pairs[nt_c];n++){
	//for(nt2=0;nt2<nt_n;nt2++){
	//if(nt_c != nt2 && nt_c-1 != nt2 && nt_c+1 !=nt2){
	nt2=vl_neighbor_tables[nt_c][n];
	typ_ind = mc_n_types*mc_types[at_c]     + mc_types[N_PARTS_PER_NT*nt2];
	typ_ind2= mc_n_types*mc_types[N_PARTS_PER_NT*nt2] + mc_types[at_c];
	
	sugdistsq=calc_min_dist_sq(rx[N_PARTS_PER_NT*nt2+ISUG], ry[N_PARTS_PER_NT*nt2+ISUG], rz[N_PARTS_PER_NT*nt2+ISUG], mc_temp_x[at_c+ISUG], mc_temp_y[at_c+ISUG], mc_temp_z[at_c+ISUG]);
	if(sugdistsq<mc_r_cut_sq){
	  calc_min_vec(rx[N_PARTS_PER_NT*nt2], ry[N_PARTS_PER_NT*nt2], rz[N_PARTS_PER_NT*nt2], mc_temp_x[at_c], mc_temp_y[at_c], mc_temp_z[at_c], d_vec, &dist);
	  proj_on_nt(d_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
	  proj_on_nt_inv(d_vec, rx, ry, rz, nt2, r_vec_inv);
	  //calc_rel_pos(d_vec, r_vec);
	  //calc_rel_pos_inv(d_vec, nt2,rx, ry,rz, r_vec_inv);
	  //theta=calc_wc_psdihedr(mc_temp_x, mc_temp_y, mc_temp_z, rx,ry,rz,nt_c,nt2); 
	  //temp_flag=0;
	  //nb_wc = 0;//MC_calc_nnN_watscric(typ_ind, r_vec, r_vec_inv, theta, &temp_flag);
	  //if(temp_flag==0){
	  //if(nb_wc<0)fprintf(mcconf," wcp %d ", nt2);//if(cflag==-1)printf( " %lf ", nb_wc);
	  //}
	  //if(MC_is_wc_corresp(nt_c, nt2)>=0) fprintf(mcconf," wcp %d ", nt2);
	  
	  MC_is_wc_corresp(nt_c, nt2, nt_n, wc_pair_array1);
	  MC_is_wc_corresp(nt2, nt_c, nt_n, wc_pair_array2);
	  temp_face1=wc_pair_array1[0];
	  temp_face2=wc_pair_array2[0];
	  
	  if(temp_face1>=0 && temp_face2>=0) {
	    fprintf(mcconf, " wcp %d (", nt2);
	    if(temp_face2==WC_FACE_SUGAR)
	      fprintf(mcconf, "S");
	    else if(temp_face2==WC_FACE_WATSCRICK)
	      fprintf(mcconf, "W");
	    else if(temp_face2==WC_FACE_HOOGSTEEN)
	      fprintf(mcconf, "H");
	    if(temp_face1==WC_FACE_SUGAR)
	      fprintf(mcconf, "S");
	    else if(temp_face1==WC_FACE_WATSCRICK)
	      fprintf(mcconf, "W");
	    else if(temp_face1==WC_FACE_HOOGSTEEN)
	      fprintf(mcconf, "H");
	    fprintf(mcconf, ")");
	  }
	  
	  temp_face1=wc_pair_array1[1];
	  temp_face2=wc_pair_array2[2];
	  if(temp_face1>=0 && temp_face2>=0) {
	    fprintf(mcconf, " bph %d (", nt2);
	    if(temp_face2==WC_FACE_SUGAR)
	      fprintf(mcconf, "S");
	    else if(temp_face2==WC_FACE_WATSCRICK)
	      fprintf(mcconf, "W");
	    else if(temp_face2==WC_FACE_HOOGSTEEN)
	      fprintf(mcconf, "H");
	    fprintf(mcconf, "P)");
	  }
	  
	  temp_face1=wc_pair_array1[2];
	  temp_face2=wc_pair_array2[1];
	  if(temp_face1>=0 && temp_face2>=0) {
	    fprintf(mcconf, " bph %d (P", nt2);
	    if(temp_face1==WC_FACE_SUGAR)
	      fprintf(mcconf, "S");
	    else if(temp_face1==WC_FACE_WATSCRICK)
	      fprintf(mcconf, "W");
	    else if(temp_face1==WC_FACE_HOOGSTEEN)
	      fprintf(mcconf, "H");
	    fprintf(mcconf, ")");
	  }
	  
	  temp_flag=0;
	  nb_st=0;
	  eta2=calc_wc_psdihedr(mc_temp_x, mc_temp_y, mc_temp_z, rx,ry,rz,nt_c,nt2); 
	  if(fabs(eta2)<STANG || fabs(eta2) > M_PI-STANG)
	    nb_st = MC_calc_nnN_stacking(typ_ind,typ_ind2, r_vec, r_vec_inv, &temp_flag); // PAIR SPECIFIC!!
	  else temp_flag=1;
	  if(temp_flag==0){
	    if(nb_st<0)fprintf(mcconf," nst %d ", nt2);// if(cflag==-1)printf( " %lf ", nb_st);
	  }
	}
      }
      
      //THE ADDITIONAL PAIR
      fprintf(mcconf,"\n");
    }
    fclose(mcconf);
  }
}

void MC_get_sp_bonded_energy(int nt_c, int nt_n, double *rx, double *ry, double *rz){
  //central
  double energ_c=0, energ_pr=0, energ_po=0;
  int gly_flag=0;
  MC_copy_nt(nt_c, rx, ry, rz);
  int tflag=0;
  //energ_c+=MC_calc_bonded_energy(nt_c, rx, ry, rz, &tflag,0,1,0,0,&gly_flag)/2.0;//we stick to the current GLYC conformation : don't allow the change between A and H since this is not in the integrator  // energ+=MC_calc_bonded_energy(nt_c, rx, ry, rz, &tflag, tflag_G1, tflag_G2, &gly_flag);
  //prev
  if(nt_c>0){
    mc_nbonds[nt_c-1][0]=1;
    mc_bondlist[nt_c-1][0]=nt_c;
    MC_copy_nt(nt_c-1, rx, ry, rz);
    tflag=0;
    //energ_pr+=MC_calc_bonded_energy(nt_c-1, rx, ry, rz, &tflag,0,1,0,0,&gly_flag)/2.0;
  }
  //post
  if(nt_c<nt_n){
    mc_nbonds[nt_c+1][0]=1;
    mc_bondlist[nt_c+1][0]=nt_c;
    MC_copy_nt(nt_c+1, rx, ry, rz);
    tflag=0;
    //energ_po+=MC_calc_bonded_energy(nt_c+1, rx, ry, rz, &tflag,0,1,0,0,&gly_flag)/2.0;
  }
  printf("%d  :  %lf \t%d  :  %lf \t%d  :  %lf\n", nt_c, energ_c, nt_c-1, energ_pr, nt_c+1, energ_po);
  
}

double MC_get_energy(int nt_n, double *rx, double *ry, double *rz, int index){
  if(index !=0 && index !=1 && index !=2 && index !=3 && index !=4) return 0;
  int nt_c, at_c, at_ne;
  double b_st=0.0, b_sg=0.0, b_ev=0.0,nb_st=0.0, nb_wc=0.0, nb_ev=0.0;
  double tot_b_st=0.0, tot_b_sg=0.0, tot_b_ev=0.0,tot_nb_st=0.0, tot_nb_wc=0.0, tot_nb_ev=0.0;
  double d_vec[DIM], r_vec[DIM], r_vec_inv[DIM];
  double dist, sugdistsq;
  int typ_ind, typ_ind2;
  int tflag=0, tflag_G1=0, tflag_G2=0, tflag_G1_pre=0, tflag_G2_pre=0;
  double dumm, dumm2;
  
  int b, nt_neigh, ba;
  double eta, theta;
  int temp_flag;
  double sigma, r;
  int n, nt2;
  double b_ev_temp1, b_ev_temp2, t_vec[DIM], nb_ev_temp1, nb_ev_temp2;
  double EintraG1=0, EintraG2=0;
  
  double etemp;
 
  double energ_c=0;
  int gly_flag=0, gly_flag_pre=0;
  for(nt_c=0;nt_c<nt_n;nt_c++){

    at_c=N_PARTS_PER_NT*nt_c;
    MC_copy_nt(nt_c, rx, ry, rz);
    tflag=0; tflag_G1=0; tflag_G2=0;
    EintraG1=0; EintraG2=0;

    /* self interaction (or intra-nt) */
    if(index==0 || index ==3){
      if(fr_is_mobile[nt_c]!=FR_MOB_FROZ){
	etemp=MC_calc_intra_energy(nt_c,&tflag_G1, &tflag_G2, &EintraG1, &EintraG2);
	energ_c+=etemp;//MC_calc_intra_energy(nt_c,&tflag);
      }
    }
    
    if(tflag_G1!=0 && (index==0 || index == 3))
      {  printf("Intra: flag = %d, glp = %d   %d   %d\n", tflag, mc_glyc[nt_c], mc_puck[nt_c], glp_is_flippable[nt_c]);exit(0);	}//   return;	}
    /* bonded interactions */
    if(index==1 || index==3){
      //froz check is done insidle the MC_calc_bonded_energy function
	energ_c+=MC_calc_bonded_energy(nt_c, rx, ry, rz, &tflag, 0, 1, 0,0,
				       0, 1, 0, 0, &gly_flag, &gly_flag_pre, &dumm,&dumm2)/2.0;

	
    }
    if(tflag!=0 && index==1)
      {  printf("Bonded flag = %d\n", tflag);exit(0);	}
    
    /* non bonded loop */
    if(index==2 || index==3){
      for(n=0;n<vl_n_pairs[nt_c];n++) {
	nt2=vl_neighbor_tables[nt_c][n];
	at_ne=N_PARTS_PER_NT*nt2;
#ifdef FROZEN
	if(fr_is_mobile[nt_c]!=FR_MOB_FROZ || fr_is_mobile[nt2]!=FR_MOB_FROZ)
#endif
	  {
	    sugdistsq=calc_min_dist_sq(rx[at_ne+ISUG], ry[at_ne+ISUG], rz[at_ne+ISUG], mc_temp_x[at_c+ISUG], mc_temp_y[at_c+ISUG], mc_temp_z[at_c+ISUG]);
	    if(sugdistsq<mc_nb_rcut_sq) {
	      calc_min_vec(rx[at_ne], ry[at_ne], rz[at_ne], mc_temp_x[at_c], mc_temp_y[at_c], mc_temp_z[at_c], r_vec, &r);
	      tflag=0;
	      //printf("%d  and %d\n", nt_c, nt2);
	      energ_c+=MC_calc_non_bonded_energy(nt_c, rx, ry, rz, nt2, r_vec, r, &tflag)/2.0;
	      if(tflag!=0 && index==2)
		{  printf("Nonbonded flag = %d between %d and %d\n", tflag, nt_c, nt2);exit(0);	   	}
	    }
	  }
      }
    }
    
}
  //ADD WATSON-CRICK!
#ifndef NOCTCS
if(index==3 || index==4) energ_c+=MC_calculate_total_wc_energy(nt_n,-1,rx, ry, rz);
#endif
  
  return energ_c;
}




void MC_update_positions(double **rx, double **ry, double **rz, int nt){
  int i;
  for(i=0;i<N_PARTS_PER_NT;i++){
    (*rx)[N_PARTS_PER_NT*nt+i]=mc_temp_x[N_PARTS_PER_NT*nt+i];
    (*ry)[N_PARTS_PER_NT*nt+i]=mc_temp_y[N_PARTS_PER_NT*nt+i];
    (*rz)[N_PARTS_PER_NT*nt+i]=mc_temp_z[N_PARTS_PER_NT*nt+i];
#ifdef PBC
    mc_pbox[N_PARTS_PER_NT*nt+i][0]=mc_temp_pbox[N_PARTS_PER_NT*nt+i][0];
    mc_pbox[N_PARTS_PER_NT*nt+i][1]=mc_temp_pbox[N_PARTS_PER_NT*nt+i][1];
    mc_pbox[N_PARTS_PER_NT*nt+i][2]=mc_temp_pbox[N_PARTS_PER_NT*nt+i][2];
#endif
  }
}

void MC_update_glypuck(int nt){
  mc_glyc[nt]=mc_temp_glyc[nt];
  mc_puck[nt]=mc_temp_puck[nt];
}

double MC_eval_displacement(int nt_n, double **rx, double **ry, double **rz, int nt, double o_ene, double n_ene, double wc_ene_nmo){
  
  double d_energ=0;
  //double pref=1.0;
  double ranf;
  double expfac=0;
  double sign=1; if(n_ene-o_ene+wc_ene_nmo<0) sign=-1; //sign of DELTA_E
  int accflag=0;
  if(mc_target_temp>0.05){
    if(sign*(n_ene-o_ene+wc_ene_nmo)/mc_target_temp < 100){
      expfac=exp(-(n_ene-o_ene+wc_ene_nmo)/mc_target_temp);
      ranf=rand_d(1.0);
      if(ranf<expfac) accflag=1;
    }
    else{
      //value overflowed - accept of reject depending on sign
      if(sign>0) accflag=0;
      else accflag=1;
    }
  }
  else{
    if(sign>0) accflag=0;
    else accflag=1;
  }
  if(accflag==1){
    //if(ranf <= exp(-(n_ene-o_ene+wc_ene_nmo)/mc_target_temp )){
    //printf("accepted\n");
    d_energ=n_ene-o_ene+wc_ene_nmo;
    MC_update_positions(rx, ry, rz, nt);
    MC_update_glypuck(nt);
    MC_add_vl_count(nt);
    MC_update_wc_lists(nt_n);
    
    if(nt>0){
      MC_update_positions(rx, ry, rz, nt-1);
      MC_update_glypuck(nt-1);
      //printf("updating %d and %d\n", nt, nt-1);
    }
#ifdef ERMSDR
    ERMSD_SQ=ERMSD_SQ+DELTA_ERMSD_SQ;
    ERMSD_ENERG=ERMSD_ENERG+DELTA_ERMSD_ENERG;
    WALL_ENERG=WALL_ENERG+DELTA_WALL_ENERG;
#endif
  }
  return d_energ;
}

void MC_perform_bb_rotation(int nt_c){
  /* rotational axis */
  int at_c=N_PARTS_PER_NT*nt_c;
  int i,d;
  double sign=1.0;
  if(rand_d(1)<0.5) sign=-1.0;
  double xi[DIM], epsilon[DIM], eta[DIM];
  double vdotR;
  double Xvec[DIM], Yvec[DIM];
  double rot_axis[DIM];
   
  double to_rot[DIM];
  double v_parall[DIM], v_perpen[DIM];
  double CM[DIM];
  int CM_pbox[DIM];
  double MASSES[N_PARTS_PER_NT], TOT_MASS=0.0;
  
  //HERE WE SET THE ROTATION POINT - NOW ITS THE SUGAR
  for(i=0;i<N_PARTS_PER_NT-1;i++)
    MASSES[i]=0.0;
  MASSES[ISUG]=1.0;
  for(i=0;i<N_PARTS_PER_NT-1;i++){
    TOT_MASS+=MASSES[i];
  }
  for(d=0;d<DIM;d++){
    CM[d]=0.0;
#ifdef PBC
    CM_pbox[d]=0;
#endif
  }
  for(i=0;i<N_PARTS_PER_NT-1;i++){
    CM[0]+=MASSES[i]*get_unf_coo_temp_x(at_c+i);
    CM[1]+=MASSES[i]*get_unf_coo_temp_y(at_c+i);
    CM[2]+=MASSES[i]*get_unf_coo_temp_z(at_c+i);
    //mc_temp_x[i]+box_l[0]*mc_temp_pbox[i][0]);
  }
  for(d=0;d<DIM;d++)
    CM[d]/=TOT_MASS;
  
#ifdef PBC
  for(d=0;d<DIM;d++)
    if(CM[d]<0) { CM_pbox[d]-- ; CM[d]+=box_l[d];}
  for(d=0;d<DIM;d++)
    if(CM[d]>=box_l[d]) { CM_pbox[d]++ ; CM[d]-=box_l[d];}
#endif
  
  //if(CM_pbox[0]!=mc_temp_pbox[3][0] || CM_pbox[1]!=mc_temp_pbox[3][1] || CM_pbox[2]!=mc_temp_pbox[3][2]){
  //printf("AAAAAAAAAAAAAAAAAA\n");
  //exit(1);}
  //printf("CM   %lf %lf %lf\nSUG  %lf %lf %lf\n", CM[0], CM[1], CM[2], mc_temp_x[3], mc_temp_y[3], mc_temp_z[3]);
  /* rot axis is perpendicular to the base plane  - X x Y*/
  Xvec[0]=dist_1d(mc_temp_x[at_c+1], mc_temp_x[at_c],0);
  Xvec[1]=dist_1d(mc_temp_y[at_c+1], mc_temp_y[at_c],1);
  Xvec[2]=dist_1d(mc_temp_z[at_c+1], mc_temp_z[at_c],2);
  Yvec[0]=dist_1d(mc_temp_x[at_c+2], mc_temp_x[at_c],0);
  Yvec[1]=dist_1d(mc_temp_y[at_c+2], mc_temp_y[at_c],1);
  Yvec[2]=dist_1d(mc_temp_z[at_c+2], mc_temp_z[at_c],2);
  
  rot_axis[0]=sign*(Xvec[1]*Yvec[2]-Xvec[2]*Yvec[1]);
  rot_axis[1]=sign*(Xvec[2]*Yvec[0]-Xvec[0]*Yvec[2]);
  rot_axis[2]=sign*(Xvec[0]*Yvec[1]-Xvec[1]*Yvec[0]);
  double m_axis=sqrt(dot_prod(rot_axis,rot_axis));
  rot_axis[0]/=m_axis;rot_axis[1]/=m_axis;rot_axis[2]/=m_axis;
  
  //printf("BB   %lf %lf %lf   %lf %lf %lf        %lf %lf %lf\n",mc_temp_x[0], mc_temp_y[0], mc_temp_z[0],mc_temp_x[1], mc_temp_y[1], mc_temp_z[1], Xvec[0], Xvec[1], Xvec[2]);
  //printf("BB   %lf %lf %lf %lf %lf\n", dot_prod(Xvec,Xvec), dot_prod(Yvec,Yvec),rot_axis[0],rot_axis[1],rot_axis[2]);
  //intf("%lf %lf %lf\n", MC_BB_ANGLE, MC_BB_ANGLE_COS, MC_BB_ANGLE_SIN);
  for(i=0;i<N_PARTS_PER_NT-2;i++){
    to_rot[0]=dist_1d(mc_temp_x[at_c+i],CM[0]     ,0);
    to_rot[1]=dist_1d(mc_temp_y[at_c+i],CM[1]     ,1);
    to_rot[2]=dist_1d(mc_temp_z[at_c+i],CM[2]     ,2);
    vdotR=dot_prod(to_rot, rot_axis);
    
    for(d=0;d<DIM;d++){
      v_parall[d]=vdotR*rot_axis[d];
      v_perpen[d]=to_rot[d]-v_parall[d];
    }
    vec_prod(v_perpen, rot_axis, eta);

    for(d=0;d<DIM;d++){
      epsilon[d]=v_perpen[d]*MC_BB_ANGLE_COS+eta[d]*MC_BB_ANGLE_SIN;
      xi[d]=v_parall[d]+epsilon[d];
    }
    mc_temp_x[at_c+i]=xi[0]+CM[0];
    mc_temp_y[at_c+i]=xi[1]+CM[1];
    mc_temp_z[at_c+i]=xi[2]+CM[2];
    
#ifdef PBC
    mc_temp_x[at_c+i]+=(CM_pbox[0])*box_l[0];
    mc_temp_y[at_c+i]+=(CM_pbox[1])*box_l[1];
    mc_temp_z[at_c+i]+=(CM_pbox[2])*box_l[2];
    
    mc_temp_pbox[at_c+i][0]=0;mc_temp_pbox[at_c+i][1]=0;mc_temp_pbox[at_c+i][2]=0;
    if(mc_temp_x[at_c+i]<0) { mc_temp_pbox[at_c+i][0]-- ; mc_temp_x[at_c+i]+=box_l[0];}
    if(mc_temp_y[at_c+i]<0) { mc_temp_pbox[at_c+i][1]-- ; mc_temp_y[at_c+i]+=box_l[1];}
    if(mc_temp_z[at_c+i]<0) { mc_temp_pbox[at_c+i][2]-- ; mc_temp_z[at_c+i]+=box_l[2];}
    
    if(mc_temp_x[at_c+i]>=box_l[0]) { mc_temp_pbox[at_c+i][0]++; mc_temp_x[at_c+i]-=box_l[0];}
    if(mc_temp_y[at_c+i]>=box_l[1]) { mc_temp_pbox[at_c+i][1]++; mc_temp_y[at_c+i]-=box_l[1];}
    if(mc_temp_z[at_c+i]>=box_l[2]) { mc_temp_pbox[at_c+i][2]++; mc_temp_z[at_c+i]-=box_l[2];}
#endif
  } 
}

void MC_perform_nt_rotation(int nt_c, int rot_type){
  /* rotational axis */
  int at_c=N_PARTS_PER_NT*nt_c;
  int i,d;
  double xi[DIM], epsilon[DIM], eta[DIM];
  double vdotR;
  double rot_axis[DIM];
  double R1=rand_d(1.0);
  double R2=rand_d(1.0);
  double phi=2*M_PI*R1;
  double rho=2*R2-1;
  double to_rot[DIM];
  double v_parall[DIM], v_perpen[DIM];
  double CM[DIM];
  int CM_pbox[DIM];
  double MASSES[N_PARTS_PER_NT], TOT_MASS=0.0;
  double Xvec[DIM], Yvec[DIM];
  //HERE WE SET THE ROTATION POINT - NOW ITS THE SUGAR
  double rang_sc=rand_d(1.0);
  double cos_angle=cos(MC_NT_ANGLE*rang_sc),sin_angle=sin(MC_NT_ANGLE*rang_sc);
  switch(rot_type){
  case 0:
    CM[0]=get_unf_coo_temp_x(at_c+ISUG);
    CM[1]=get_unf_coo_temp_y(at_c+ISUG);
    CM[2]=get_unf_coo_temp_z(at_c+ISUG);
    break;
  case 1:
    CM[0]=get_unf_coo_temp_x(at_c+IBAS);
    CM[1]=get_unf_coo_temp_y(at_c+IBAS);
    CM[2]=get_unf_coo_temp_z(at_c+IBAS);
    break;
      
  case 2:
    CM[0]=0.5*(get_unf_coo_temp_x(at_c+IBAS)+get_unf_coo_temp_x(at_c+ISUG));
    CM[1]=0.5*(get_unf_coo_temp_y(at_c+IBAS)+get_unf_coo_temp_y(at_c+ISUG));
    CM[2]=0.5*(get_unf_coo_temp_z(at_c+IBAS)+get_unf_coo_temp_z(at_c+ISUG));
    break;
      
  default:
    printf("Invalid rotation type  %d!\n", rot_type); exit(ERR_INTEG);
    break;
  }
  /* if(rot_type==0){ */
  /*   for(i=0;i<N_PARTS_PER_NT-1;i++) */
  /*     MASSES[i]=0.0; */
  /*   MASSES[ISUG]=1.0; */
  /* } else if(rot_type==1){ //NOW ITS THE BASE CENTER */
  /*   for(i=0;i<N_PARTS_PER_NT-1;i++) */
  /*     MASSES[i]=1.0; */
  /* } */
  /* else {printf("Invalid rotation type  %d!\n", rot_type); exit(1);} */
  /* if(rot_type==0){ */
  /*   for(i=0;i<N_PARTS_PER_NT-1;i++) */
  /*     MASSES[i]=0.0; */
  /*   MASSES[ISUG]=1.0; */
  /* } else if(rot_type==1){ //NOW ITS THE BASE CENTER */
  /*   for(i=0;i<N_PARTS_PER_NT-1;i++) */
  /*     MASSES[i]=1.0; */
  /* } */
  /* else {printf("Invalid rotation type!\n"); exit(1);}  */
 /*  for(i=0;i<N_PARTS_PER_NT-1;i++){ */
/*     TOT_MASS+=MASSES[i]; */
/*   } */
/*   for(d=0;d<DIM;d++){ */
/*     CM[d]=0.0; */
/* #ifdef PBC */
/*     CM_pbox[d]=0.0; */
/* #endif */
/*   } */
  
/*   for(i=0;i<N_PARTS_PER_NT-1;i++){ */
/*     CM[0]+=MASSES[i]*get_unf_coo_temp_x(at_c+i); */
/*       //(mc_temp_x[i]+box_l[0]*mc_temp_pbox[i][0]); */
/*     CM[1]+=MASSES[i]*get_unf_coo_temp_y(at_c+i); */
/*     CM[2]+=MASSES[i]*get_unf_coo_temp_z(at_c+i); */
/*   } */
/*   for(d=0;d<DIM;d++) */
/*     CM[d]/=TOT_MASS; */
  
 #ifdef PBC 
  for(d=0;d<DIM;d++) 
    if(CM[d]<0) { CM_pbox[d]-- ; CM[d]+=box_l[d];} 
  for(d=0;d<DIM;d++) 
    if(CM[d]>=box_l[d]) { CM_pbox[d]++ ; CM[d]-=box_l[d];} 
#endif 
  
  rot_axis[0]=cos(phi)*sqrt(1-rho*rho);
  rot_axis[1]=sin(phi)*sqrt(1-rho*rho);
  rot_axis[2]=rho;
  
  Xvec[0]=dist_1d(mc_temp_x[at_c+1], mc_temp_x[at_c],0);
  Xvec[1]=dist_1d(mc_temp_y[at_c+1], mc_temp_y[at_c],1);
  Xvec[2]=dist_1d(mc_temp_z[at_c+1], mc_temp_z[at_c],2);
  Yvec[0]=dist_1d(mc_temp_x[at_c+2], mc_temp_x[at_c],0);
  Yvec[1]=dist_1d(mc_temp_y[at_c+2], mc_temp_y[at_c],1);
  Yvec[2]=dist_1d(mc_temp_z[at_c+2], mc_temp_z[at_c],2);
  
  //printf("NT   %lf %lf %lf   %lf %lf %lf        %lf %lf %lf\n",mc_temp_x[0], mc_temp_y[0], mc_temp_z[0],mc_temp_x[1], mc_temp_y[1], mc_temp_z[1], Xvec[0], Xvec[1], Xvec[2]);
  //printf("NT   %lf %lf %lf\n", dot_prod(Xvec,Xvec), dot_prod(Yvec,Yvec),rot_axis[0]*rot_axis[0]+ rot_axis[1]*rot_axis[1]+ rot_axis[2]*rot_axis[2]);
  
  
  for(i=0;i<N_PARTS_PER_NT-1;i++){
    to_rot[0]=dist_1d(mc_temp_x[at_c+i],CM[0]     ,0);
    to_rot[1]=dist_1d(mc_temp_y[at_c+i],CM[1]     ,1);
    to_rot[2]=dist_1d(mc_temp_z[at_c+i],CM[2]     ,2);
    vdotR=dot_prod(to_rot, rot_axis);
    
    for(d=0;d<DIM;d++){
      v_parall[d]=vdotR*rot_axis[d];
      v_perpen[d]=to_rot[d]-v_parall[d];
    }
    vec_prod(v_perpen, rot_axis, eta);
    
    for(d=0;d<DIM;d++){
      epsilon[d]=v_perpen[d]*cos_angle+eta[d]*sin_angle;
      xi[d]=v_parall[d]+epsilon[d];
    }
    mc_temp_x[at_c+i]=xi[0]+CM[0];
    mc_temp_y[at_c+i]=xi[1]+CM[1];
    mc_temp_z[at_c+i]=xi[2]+CM[2];
    
#ifdef PBC
    mc_temp_x[at_c+i]+=(CM_pbox[0]-mc_temp_pbox[at_c+i][0])*box_l[0];
    mc_temp_y[at_c+i]+=(CM_pbox[1]-mc_temp_pbox[at_c+i][1])*box_l[1];
    mc_temp_z[at_c+i]+=(CM_pbox[2]-mc_temp_pbox[at_c+i][2])*box_l[2];
    if(mc_temp_x[at_c+i]<0) { mc_temp_pbox[at_c+i][0]-- ; mc_temp_x[at_c+i]+=box_l[0];}
    if(mc_temp_y[at_c+i]<0) { mc_temp_pbox[at_c+i][1]-- ; mc_temp_y[at_c+i]+=box_l[1];}
    if(mc_temp_z[at_c+i]<0) { mc_temp_pbox[at_c+i][2]-- ; mc_temp_z[at_c+i]+=box_l[2];}
    
    if(mc_temp_x[at_c+i]>=box_l[0]) { mc_temp_pbox[at_c+i][0]++; mc_temp_x[at_c+i]-=box_l[0];}
    if(mc_temp_y[at_c+i]>=box_l[1]) { mc_temp_pbox[at_c+i][1]++; mc_temp_y[at_c+i]-=box_l[1];}
    if(mc_temp_z[at_c+i]>=box_l[2]) { mc_temp_pbox[at_c+i][2]++; mc_temp_z[at_c+i]-=box_l[2];}
#endif
  }
}


void MC_perform_nt_translation(int nt_c){
  int i;
  int at_c=N_PARTS_PER_NT*nt_c;
  double ranx=rand_d(1.0);
  double rany=rand_d(1.0);
  double ranz=rand_d(1.0);
  for(i=0;i<N_PARTS_PER_NT-1;i++){
    mc_temp_x[at_c+i]+=(ranx-0.5)*MC_NT_XYZ;
    mc_temp_y[at_c+i]+=(rany-0.5)*MC_NT_XYZ;
    mc_temp_z[at_c+i]+=(ranz-0.5)*MC_NT_XYZ;

#ifdef PBC
    if(mc_temp_x[at_c+i]<0) { mc_temp_pbox[at_c+i][0]-- ; mc_temp_x[at_c+i]+=box_l[0];}
    if(mc_temp_y[at_c+i]<0) { mc_temp_pbox[at_c+i][1]-- ; mc_temp_y[at_c+i]+=box_l[1];}
    if(mc_temp_z[at_c+i]<0) { mc_temp_pbox[at_c+i][2]-- ; mc_temp_z[at_c+i]+=box_l[2];}
    
    if(mc_temp_x[at_c+i]>=box_l[0]) { mc_temp_pbox[at_c+i][0]++; mc_temp_x[at_c+i]-=box_l[0];}
    if(mc_temp_y[at_c+i]>=box_l[1]) { mc_temp_pbox[at_c+i][1]++; mc_temp_y[at_c+i]-=box_l[1];}
    if(mc_temp_z[at_c+i]>=box_l[2]) { mc_temp_pbox[at_c+i][2]++; mc_temp_z[at_c+i]-=box_l[2];}
#endif
  }
}

void MC_perform_p_translation(int nt_c){
  int i;
  int at_c=N_PARTS_PER_NT*nt_c;
  double ranx=rand_d(1.0);
  double rany=rand_d(1.0);
  double ranz=rand_d(1.0);
  mc_temp_x[at_c+IPHO]+=(ranx-0.5)*MC_PH_XYZ;
  mc_temp_y[at_c+IPHO]+=(rany-0.5)*MC_PH_XYZ;
  mc_temp_z[at_c+IPHO]+=(ranz-0.5)*MC_PH_XYZ;
  
#ifdef PBC
  if(mc_temp_x[at_c+IPHO]<0) { mc_temp_pbox[at_c+IPHO][0]-- ; mc_temp_x[at_c+IPHO]+=box_l[0];}
  if(mc_temp_y[at_c+IPHO]<0) { mc_temp_pbox[at_c+IPHO][1]-- ; mc_temp_y[at_c+IPHO]+=box_l[1];}
  if(mc_temp_z[at_c+IPHO]<0) { mc_temp_pbox[at_c+IPHO][2]-- ; mc_temp_z[at_c+IPHO]+=box_l[2];}
  
  if(mc_temp_x[at_c+IPHO]>=box_l[0]) { mc_temp_pbox[at_c+IPHO][0]++; mc_temp_x[at_c+IPHO]-=box_l[0];}
  if(mc_temp_y[at_c+IPHO]>=box_l[1]) { mc_temp_pbox[at_c+IPHO][1]++; mc_temp_y[at_c+IPHO]-=box_l[1];}
  if(mc_temp_z[at_c+IPHO]>=box_l[2]) { mc_temp_pbox[at_c+IPHO][2]++; mc_temp_z[at_c+IPHO]-=box_l[2];}
#endif
}

int MC_perform_trial_movement(int nt_c){
  int maxmoves=9;
  //printf("performing trial mov\n");
  //printf("glycs =  %d %d %d\n", mc_glyc[0],mc_glyc[1],mc_glyc[2]);
  //int maxmoves=7;//THIS IS TEMPORARY !!! THIS LINE SHOULD NOT BE UNCOMMENTED!!!
  
  //if(glp_is_flippable[nt_c]==GLP_FIXED || mc_temp_glyc[nt_c]==GLYC_H)
  if(glp_is_flippable[nt_c]==GLP_FIXED)
    maxmoves=maxmoves-2;
  int rand1=rand_i(maxmoves);
#ifdef FROZEN
  if(fr_is_mobile[nt_c]==FR_MOB_BASE) //only NT is mobile - the flag has to be > 0 
    rand1=rand_i(maxmoves-1)+1;
  else if(fr_is_mobile[nt_c]==FR_MOB_PHOS)
    rand1=0;
#endif
  if(rand1==0 || rand1==maxmoves){
    MC_perform_p_translation(nt_c);
    //MC_perform_nt_translation(nt_c);
  }
  else if(rand1==1){
    MC_perform_nt_rotation(nt_c,0); // ROTATION AROUND SUGAR
  }
  else if(rand1==2){
    MC_perform_nt_rotation(nt_c,1); // ROTATION AROUND BASE
  }
  else if(rand1==3){
    MC_perform_nt_rotation(nt_c,2); // ROTATION AROUND (SUG+BAS)/2
  }
  else if(rand1==4 || rand1==5 || rand1==6){
    MC_perform_nt_translation(nt_c); // TRANSLATION OF NT
  }
  else if(rand1==7){
    if(glp_is_flippable[nt_c]==GLP_PUCK){
      MC_perform_base_flip(nt_c,FLIP_PUCK);
    }
    else{
      if(is_pyrimidine(nt_c)){
	if(glp_is_flippable[nt_c]==GLP_BOTH)
	  MC_perform_base_flip(nt_c,FLIP_PUCK);
	else
	  MC_perform_nt_translation(nt_c); // do something else when glyc is the only change you can make
      }
      else
	MC_perform_base_flip(nt_c,FLIP_GLYC);
    }
  }
  else if(rand1==8){
    if(glp_is_flippable[nt_c]==GLP_GLYC){
      if(is_purine(nt_c))
	MC_perform_base_flip(nt_c,FLIP_GLYC);
      else
	MC_perform_nt_translation(nt_c); // do something else when glyc is the only change you can make
      
    }
    else{
      MC_perform_base_flip(nt_c,FLIP_PUCK);
    }
  }
  return rand1;
}

void MC_perform_base_flip(int nt_c, int flip_type){
  //we displace the base
  //then flip it
  //then update the glycosidic index
  //then remap the sugar
  int i,j,d,m,n;
  double xvec[DIM], yvec[DIM], zvec[DIM];
  double newxvec[DIM], newyvec[DIM], newzvec[DIM];
  double mat_U[DIMSQ], mat_Utr[DIMSQ], mat_A[DIMSQ], mat_B[DIMSQ];
  double vec_a[DIM], vec_b[DIM];
  //double **mat_A, **mat_U, **mat_Utr, *vec_a, *vec_b;
  //choose parameters
  int at_c=nt_c*N_PARTS_PER_NT;
  //printf("Flipping %d\t", nt_c);
  int newconf=MC_init_flip(mc_types[at_c], mc_temp_glyc[nt_c],mc_temp_puck[nt_c], flip_type, mat_A, vec_a); // and we have set the rotation matrix and the displacement vector
  xvec[0]=mc_temp_x[at_c+IX]-mc_temp_x[at_c];  xvec[1]=mc_temp_y[at_c+IX]-mc_temp_y[at_c];  xvec[2]=mc_temp_z[at_c+IX]-mc_temp_z[at_c];
  yvec[0]=mc_temp_x[at_c+IY]-mc_temp_x[at_c];  yvec[1]=mc_temp_y[at_c+IY]-mc_temp_y[at_c];  yvec[2]=mc_temp_z[at_c+IY]-mc_temp_z[at_c];
  vec_prod(xvec, yvec, zvec);
  mat_U[0*DIM+0]=xvec[0];  mat_U[0*DIM+1]=xvec[1];  mat_U[0*DIM+2]=xvec[2];
  mat_U[1*DIM+0]=yvec[0];  mat_U[1*DIM+1]=yvec[1];  mat_U[1*DIM+2]=yvec[2];
  mat_U[2*DIM+0]=zvec[0];  mat_U[2*DIM+1]=zvec[1];  mat_U[2*DIM+2]=zvec[2];
  for(i=0;i<DIM;i++){
    for(j=0;j<DIM;j++){
      mat_Utr[i*DIM+j]=mat_U[j*DIM+i];
      mat_B[i*DIM+j]=0;
    }
    vec_b[i]=0;
    newxvec[i]=0;
    newyvec[i]=0;
    newzvec[i]=0;
  }
  for(i=0;i<DIM;i++)
    for(j=0;j<DIM;j++)
      for(m=0;m<DIM;m++)
	for(n=0;n<DIM;n++)
	  mat_B[i*DIM+j]+=mat_Utr[i*DIM+m]*mat_A[m*DIM+n]*mat_U[n*DIM+j];
  for(d=0;d<DIM;d++)
    for(i=0;i<DIM;i++)
      for(j=0;j<DIM;j++)
      //vec_b[d]+=mat_B[d][i]*vec_a[i];
      //vec_b[d]+=mat_B[d*DIM+i]*mat_A[i*DIM+j]*vec_a[j];
	vec_b[d]+=mat_Utr[d*DIM+i]*mat_A[i*DIM+j]*vec_a[j];
  
  //so we have the new nucleoside orientation
  for(d=0;d<DIM;d++)
    for(i=0;i<DIM;i++){
      newxvec[d]+=mat_B[d*DIM+i]*xvec[i];
      newyvec[d]+=mat_B[d*DIM+i]*yvec[i];
      //newzvec[d]+=mat_B[d][i]*zvec[i];
    }
  //normalize these vectors
  double normx=sqrt(SQ(newxvec[0])+SQ(newxvec[1])+SQ(newxvec[2]));
  double normy=sqrt(SQ(newyvec[0])+SQ(newyvec[1])+SQ(newyvec[2]));
  for(d=0;d<DIM;d++){
    newxvec[d]/=normx;
    newyvec[d]/=normy;
  }
  double xyproj=newxvec[0]*newyvec[0]+newxvec[1]*newyvec[1]+newxvec[2]*newyvec[2];
  for(d=0;d<DIM;d++)
    newyvec[d]=newyvec[d]-newxvec[d]*xyproj;
  normy=sqrt(SQ(newyvec[0])+SQ(newyvec[1])+SQ(newyvec[2]));
  for(d=0;d<DIM;d++)
    newyvec[d]/=normy;
  
  
  
  //and the new base position
  mc_temp_x[at_c]+=vec_b[0];  mc_temp_y[at_c]+=vec_b[1];  mc_temp_z[at_c]+=vec_b[2];
  mc_temp_x[at_c+IX]=mc_temp_x[at_c]+newxvec[0];  mc_temp_y[at_c+IX]=mc_temp_y[at_c]+newxvec[1];  mc_temp_z[at_c+IX]=mc_temp_z[at_c]+newxvec[2];
  mc_temp_x[at_c+IY]=mc_temp_x[at_c]+newyvec[0];  mc_temp_y[at_c+IY]=mc_temp_y[at_c]+newyvec[1];  mc_temp_z[at_c+IY]=mc_temp_z[at_c]+newyvec[2];
  
  //we update the glyc bond/sugar pucker
  MC_assign_temp_glp(nt_c, newconf);
  //we remap the sugar with the "new" pucker/glyc
  MC_map_sugar_temp(nt_c);
  
}

int MC_init_flip(int bastyp, int oldglyc, int oldpuck, int flip_type, double *matrix, double *vec){
  int newglyc, newpuck;
  //this returns the new glyc-puck conformation and selects the proper matrices and displacement vectors
  
  //TEMPORARY, WE treat GLUC_H AS GLYC_A. THE JUMP SYN->HIGH ANTI IS NOT CONSIDERED, SINCE HIGH ANTI IS A "SUBSPACE" OF ANTI 
  /* if(oldglyc==GLYC_H){ */
  /*   printf("Trying to flip high-anti.\n"); */
  /*   exit(1); */
  /* } */
  //printf("In init flip, from %d , %d\n", oldglyc, oldpuck); 
  if(flip_type==FLIP_GLYC){
    newpuck=oldpuck;
    //now, we allow the A-H -> S change, but with the same parameters of A->S change
    if(oldpuck==PUCK_3){
      if(oldglyc==GLYC_A || oldglyc==GLYC_H){
	newglyc=GLYC_S;
	if(bastyp==TYP_ADENINE || bastyp==TYP_GUANINE){
	  //{{-0.452135, 0.0978822, -0.886562}, {0.372873,   0.923683, -0.0881798}, {0.810271, -0.370445, -0.454127}}
	  //{-0.645402, -1.32813, 1.45666}
	  matrix[0*DIM+0]=-0.452135;	matrix[0*DIM+1]= 0.0978822;     matrix[0*DIM+2]=-0.886562;
	  matrix[1*DIM+0]= 0.372873;	matrix[1*DIM+1]= 0.923683;	matrix[1*DIM+2]=-0.0881798;
	  matrix[2*DIM+0]= 0.810271;	matrix[2*DIM+1]=-0.370445;	matrix[2*DIM+2]=-0.454127;
	  
	  vec[0]=-0.645402;	vec[1]=-1.32813;	vec[2]= 1.45666;
	}
	else if(bastyp==TYP_CYTOSINE || bastyp==TYP_URACIL){
	  printf("i shouldnt be here!\n");
	  exit(1);
	}
      }
      else if(oldglyc==GLYC_S){
	newglyc=GLYC_A;
	if(bastyp==TYP_ADENINE || bastyp==TYP_GUANINE){
	  //{{-0.452135, 0.372873, 0.810271}, {0.0978822,   0.923683, -0.370445}, {-0.886562, -0.0881798, -0.454127}}
	  //{1.12961, 1.59587, 0.692458}
	  matrix[0*DIM+0]=-0.452135;	matrix[0*DIM+1]=0.372873;	matrix[0*DIM+2]= 0.810271;
	  matrix[1*DIM+0]=0.0978822;	matrix[1*DIM+1]=0.923683;	matrix[1*DIM+2]=-0.370445;
	  matrix[2*DIM+0]=-0.886562;	matrix[2*DIM+1]=-0.0881798;     matrix[2*DIM+2]= -0.454127;
	  
	  vec[0]= 1.12961;	vec[1]= 1.59587;	vec[2]= 0.692458;
	}
	else if(bastyp==TYP_CYTOSINE || bastyp==TYP_URACIL){
	  printf("i shouldnt be here!\n");
	  exit(1);
	}
      }
    }
    else if(oldpuck==PUCK_2){
      if(oldglyc==GLYC_A || oldglyc==GLYC_H){
	newglyc=GLYC_S;
	if(bastyp==TYP_ADENINE || bastyp==TYP_GUANINE){
	  //{{-0.452135, 0.0978822, -0.886562}, {0.372873,   0.923683, -0.0881798}, {0.810271, -0.370445, -0.454127}}
	  //{-0.645402, -1.32813, 1.45666}
	  //puck 2{{-0.917397, -0.199475, -0.344373}, {-0.162731,  0.977691, -0.13281}, {0.363183, -0.0657992, -0.929392}}
	  //{2.1465, 0.16011, 1.19316}
	  matrix[0*DIM+0]=-0.917397 ;	matrix[0*DIM+1]=-0.199475;      matrix[0*DIM+2]=-0.344373;
	  matrix[1*DIM+0]=-0.162731 ;	matrix[1*DIM+1]= 0.977691;	matrix[1*DIM+2]= -0.13281;
	  matrix[2*DIM+0]= 0.363183 ;	matrix[2*DIM+1]=-0.0657992;	matrix[2*DIM+2]=-0.929392;
	  
	  vec[0]=2.1465;	vec[1]=0.16011;	      vec[2]= 1.19316;
	}
	else if(bastyp==TYP_CYTOSINE || bastyp==TYP_URACIL){
	  printf("i shouldnt be here!\n");
	  exit(1);
	}
      }
      else if(oldglyc==GLYC_S){
	newglyc=GLYC_A;
	if(bastyp==TYP_ADENINE || bastyp==TYP_GUANINE){
	  //{{-0.452135, 0.372873, 0.810271}, {0.0978822,   0.923683, -0.370445}, {-0.886562, -0.0881798, -0.454127}}
	  //{1.12961, 1.59587, 0.692458}
	  //puck2 {{-0.917397, -0.162731, 0.363183}, {-0.199475,   0.977691, -0.0657992}, {-0.344373, -0.13281, -0.929392}}
	  //{2.41203, 0.351229, 0.339875}
	  matrix[0*DIM+0]=-0.917397;	matrix[0*DIM+1]=-0.162731;	matrix[0*DIM+2]=0.363183;
	  matrix[1*DIM+0]=-0.199475;	matrix[1*DIM+1]= 0.977691;	matrix[1*DIM+2]=-0.0657992;
	  matrix[2*DIM+0]=-0.344373;	matrix[2*DIM+1]=-0.13281;       matrix[2*DIM+2]=-0.929392;
	  
	  vec[0]= 2.41203;	vec[1]= 0.351229;	vec[2]= 0.339875;
	}
	else if(bastyp==TYP_CYTOSINE || bastyp==TYP_URACIL){
	  printf("i shouldnt be here!\n");
	  exit(1);
	}
      }
      
    
    }
  }
  else if(flip_type==FLIP_PUCK){
    newglyc=oldglyc;
    if(oldglyc==GLYC_A){//WE ASUMME THAT THE GLYCOSIDIC BOND ANGLE CONFORMATION IS ANTI!!
      if(oldpuck==PUCK_3){
	newpuck=PUCK_2;
	if(bastyp==TYP_ADENINE || bastyp==TYP_GUANINE){
	  //{{0.921656, -0.0221719, -0.387374}, {0.232272, 0.831246,   0.505054}, {0.310805, -0.555462, 0.771273}}
	  //{-1.03938, 0.973513, 0.355191}
	  matrix[0*DIM+0]=0.921656;	matrix[0*DIM+1]=-0.0221719;	matrix[0*DIM+2]=-0.387374;
	  matrix[1*DIM+0]=0.232272;	matrix[1*DIM+1]=0.831246;	matrix[1*DIM+2]=0.505054;
	  matrix[2*DIM+0]=0.310805;	matrix[2*DIM+1]=-0.555462;	matrix[2*DIM+2]=0.771273;
	  
	  vec[0]=-1.03938;	vec[1]=0.973513;	vec[2]=0.355191;
	}
	else if(bastyp==TYP_CYTOSINE || bastyp==TYP_URACIL){
	  //{{0.859784, 0.20379, -0.468233}, {0.119701, 0.810944,   0.572748}, {0.496431, -0.548488, 0.672843}}
	  //{-0.860058, 0.391386, -1.00444}
	  matrix[0*DIM+0]=0.859784;	matrix[0*DIM+1]=0.20379;	matrix[0*DIM+2]= -0.468233;
	  matrix[1*DIM+0]=0.119701;	matrix[1*DIM+1]=0.810944;	matrix[1*DIM+2]=0.572748;
	  matrix[2*DIM+0]=0.496431;	matrix[2*DIM+1]=-0.548488;	matrix[2*DIM+2]=0.672843;
	  
	  vec[0]=-0.860058;	vec[1]=0.391386;	vec[2]=-1.00444;
	}
      }
      else if(oldpuck==PUCK_2){
	newpuck=PUCK_3;
	if(bastyp==TYP_ADENINE || bastyp==TYP_GUANINE){
	  //{{0.921656, 0.232272, 0.310805}, {-0.0221719,  0.831246, -0.555462}, {-0.387374, 0.505054, 0.771273}}
	  //{1.11713, -0.7472, 0.589845}
	  matrix[0*DIM+0]=0.921656;	matrix[0*DIM+1]=0.232272;	matrix[0*DIM+2]=0.310805;
	  matrix[1*DIM+0]=-0.0221719;	matrix[1*DIM+1]=0.831246;	matrix[1*DIM+2]=-0.555462;
	  matrix[2*DIM+0]=-0.387374;	matrix[2*DIM+1]=0.505054;	matrix[2*DIM+2]=0.771273;
	  
	  vec[0]=1.11713;	vec[1]=-0.7472;	vec[2]=0.589845;
	}
	else if(bastyp==TYP_CYTOSINE || bastyp==TYP_URACIL){
	  //{{0.859784, 0.119701, 0.496431}, {0.20379,   0.810944, -0.548488}, {-0.468233, 0.572748, 0.672843}}
	  //{0.189393, 0.360848, 1.31746}
	  matrix[0*DIM+0]=0.859784;	 matrix[0*DIM+1]=0.119701;	 matrix[0*DIM+2]=0.496431;
	  matrix[1*DIM+0]=0.20379;	 matrix[1*DIM+1]=0.810944;	 matrix[1*DIM+2]=-0.548488;
	  matrix[2*DIM+0]=-0.468233;	 matrix[2*DIM+1]=0.572748;	 matrix[2*DIM+2]= 0.672843;
	  
	  vec[0]=0.189393;	 vec[1]=0.360848;	 vec[2]= 1.31746;
	}
      }
    }
    else if(oldglyc==GLYC_H){
       if(oldpuck==PUCK_3){
	newpuck=PUCK_2;
	if(bastyp==TYP_ADENINE || bastyp==TYP_GUANINE){
	  //{{0.726084, 0.494585, -0.477691}, {-0.398939, 0.868843,   0.293189}, {0.560045, -0.0223101, 0.828162}}
	  //{1.07218, -0.41053, -0.231989}
	  matrix[0*DIM+0]=0.726084;	matrix[0*DIM+1]=0.494585;	matrix[0*DIM+2]=-0.477691;
	  matrix[1*DIM+0]=-0.398939;	matrix[1*DIM+1]=0.868843;	matrix[1*DIM+2]=0.293189;
	  matrix[2*DIM+0]=0.560045;	matrix[2*DIM+1]=-0.0223101;	matrix[2*DIM+2]=0.828162;
	  
	  vec[0]=1.07218;	vec[1]=-0.41053;	vec[2]=-0.231989;
	}
	else if(bastyp==TYP_CYTOSINE || bastyp==TYP_URACIL){
	  //{{0.608478, 0.116058, -0.785038}, {-0.268734,  0.960935, -0.0662318}, {0.746684, 0.251267, 0.615896}}
	  //{-0.612358, -0.0141759, 0.503744}
	  matrix[0*DIM+0]=0.608478;	matrix[0*DIM+1]=0.116058;	matrix[0*DIM+2]=-0.785038;
	  matrix[1*DIM+0]=-0.268734;	matrix[1*DIM+1]=0.960935;	matrix[1*DIM+2]=-0.0662318;
	  matrix[2*DIM+0]=0.746684;	matrix[2*DIM+1]=0.251267;	matrix[2*DIM+2]=0.615896;
	  
	  vec[0]=-0.612358;	vec[1]=-0.0141759;	vec[2]=0.503744;
	}
      }
      else if(oldpuck==PUCK_2){
	newpuck=PUCK_3;
	if(bastyp==TYP_ADENINE || bastyp==TYP_GUANINE){
	  //{{0.726084, -0.398939, 0.560045}, {0.494585,   0.868843, -0.0223101}, {-0.477691, 0.293189, 0.828162}}
	  //{-0.686272, 0.852439, -0.417506}
	  matrix[0*DIM+0]=0.726084;	matrix[0*DIM+1]=-0.398939;	matrix[0*DIM+2]=0.560045;
	  matrix[1*DIM+0]=0.494585;	matrix[1*DIM+1]=0.868843;	matrix[1*DIM+2]=-0.0223101;
	  matrix[2*DIM+0]=-0.477691;	matrix[2*DIM+1]=0.293189;	matrix[2*DIM+2]=0.828162;
	  
	  vec[0]=-0.686272;	vec[1]= 0.852439;	vec[2]=-0.417506;
	}
	else if(bastyp==TYP_CYTOSINE || bastyp==TYP_URACIL){
	  //{{0.608478, -0.268734, 0.746684}, {0.116058, 0.960935,   0.251267}, {-0.785038, -0.0662318, 0.615896}}
	  //{0.76971, -0.117576, 0.150545}
	  matrix[0*DIM+0]=0.608478;	 matrix[0*DIM+1]=-0.268734;	 matrix[0*DIM+2]=0.746684;
	  matrix[1*DIM+0]=0.116058;	 matrix[1*DIM+1]=0.960935;	 matrix[1*DIM+2]=0.251267;
	  matrix[2*DIM+0]=-0.785038;	 matrix[2*DIM+1]=-0.0662318;	 matrix[2*DIM+2]=0.615896;
	  
	  vec[0]=0.76971;	 vec[1]=-0.117576;	 vec[2]=0.150545;
	}
      }
    }
    else if(oldglyc==GLYC_S){
      if(oldpuck==PUCK_3){
	newpuck=PUCK_2;
	if(bastyp==TYP_ADENINE || bastyp==TYP_GUANINE){
	  //{{0.588908, -0.806364, -0.0544508}, {0.456939,   0.387771, -0.800525}, {0.666629, 0.446555, 0.59682}}
	  //{-1.34632, 2.41762, 4.21156}
	  matrix[0*DIM+0]=0.588908;	matrix[0*DIM+1]=-0.806364;	matrix[0*DIM+2]=-0.0544508;
	  matrix[1*DIM+0]=0.456939;	matrix[1*DIM+1]=0.387771;	matrix[1*DIM+2]=-0.800525;
	  matrix[2*DIM+0]=0.666629;	matrix[2*DIM+1]=0.446555;	matrix[2*DIM+2]=0.59682;
	  
	  vec[0]=-1.34632;	vec[1]=2.41762;	vec[2]=4.21156;
	}
      }
      else if(oldpuck==PUCK_2){
	newpuck=PUCK_3;
	if(bastyp==TYP_ADENINE || bastyp==TYP_GUANINE){
	  //{{0.588908, 0.456939, 0.666629}, {-0.806364, 0.387771,   0.446555}, {-0.0544508, -0.800525, 0.59682}}
	  //{2.97166, 3.04916, -2.69565}
	  matrix[0*DIM+0]=0.588908;	matrix[0*DIM+1]=0.456939;	matrix[0*DIM+2]=0.666629;
	  matrix[1*DIM+0]=-0.806364;	matrix[1*DIM+1]= 0.387771;	matrix[1*DIM+2]=0.446555;
	  matrix[2*DIM+0]=-0.0544508;	matrix[2*DIM+1]=-0.800525;	matrix[2*DIM+2]=0.59682;
	  
	  vec[0]=2.97166;	vec[1]=3.04916;	vec[2]=-2.69565;
	}
      }
    }
  }
  else {
    printf("Invalid FLIP type.\n");
    exit(1);
  }
  
  return newpuck*N_GLYC_STATES+newglyc;
}


int MC_get_temp_glp(int nt, int conf){
  return mc_temp_puck[nt]*N_GLYC_STATES + mc_temp_glyc[nt];
}

void MC_assign_temp_glp(int nt, int conf){
  int temp=conf/N_GLYC_STATES;
  //if(temp==0) mc_temp_puck[nt]=PUCK_2; else mc_temp_puck[nt]=PUCK_3;
  if(temp==0) mc_temp_puck[nt]=PUCK_3; else mc_temp_puck[nt]=PUCK_2;
  temp=conf-temp*N_GLYC_STATES;
  if(temp==GLYC_A) mc_temp_glyc[nt]=GLYC_A;  else if(temp==GLYC_H) mc_temp_glyc[nt]=GLYC_H;  else if(temp==GLYC_S) mc_temp_glyc[nt]=GLYC_S; else{printf("ERROR: GLYC TYPE NOT RECOGNIZED!\n"); exit(1);}
}

int MC_calculate_local_energy(double *rx, double *ry, double *rz, int nt_c, double *energ_calc, int nt_n, int trial){
   //NO CELL STRUCTURE IMPLEMENTED YET
   /* for three dimensions! */
  int n, nt2, at_c=N_PARTS_PER_NT*nt_c, at_ne,  nt_pre=nt_c-1, at_pre=(nt_c-1)*N_PARTS_PER_NT;
  double r, r_vec[DIM], centdistsq;
  double TEMP_ERMSD_SQ=0,TEMP_ERMSD_ENERG=0;
#ifdef ERMSDR
  double e_vec[DIM], e_vec_inv[DIM];
  double termsd_p, termsd_q;
  //double ermsd_energ;
  double wall_energy=MC_wall_energy(mc_temp_x[at_c+IPHO],mc_temp_y[at_c+IPHO],mc_temp_z[at_c+IPHO]);
#endif  
  double energ=0.0;
  int flag=0, tflag=0;
  double EintraG1, EintraG2, etemp, etemp_PRE;
  double EintraG1_PRE, EintraG2_PRE;
  int iflag_G1=0, iflag_G2=0, gly_flag=-1, eval, temp_flag_G1, temp_flag_G2;
  int iflag_G1_PRE=0, iflag_G2_PRE=0, gly_flag_PRE=-1,  temp_flag_G1_PRE, temp_flag_G2_PRE;
  double next_intra=0, next_inter=0, t_vec[DIM], mc_ev_glob_rcut_sq, bond_bp_G1, bond_bp_G2, bond_bp_G1_PRE, bond_bp_G2_PRE;
  double pr_vec[DIM], pr_vec_inv[DIM], cr_vec[DIM], cr_vec_inv[DIM], cd_vec[DIM];
  double tnb_energ;
  /*********************/
  etemp=MC_calc_intra_energy(nt_c,&iflag_G1, &iflag_G2, &EintraG1, &EintraG2); //WE SELECT LATER THE CORRECT ENERGY
  if(nt_c>0){
    MC_copy_nt(nt_pre, rx, ry, rz);
    etemp_PRE=MC_calc_intra_energy(nt_pre,&iflag_G1_PRE, &iflag_G2_PRE, &EintraG1_PRE, &EintraG2_PRE); //WE SELECT LATER THE CORRECT ENERGY, for the PREVIOUS nucleotide
    if(fr_is_mobile[nt_pre]!=FR_MOB_FULL){
      //if the previous nt is frozen, we do not modify its glp state
      iflag_G2_PRE=1;
      EintraG2_PRE=0;
    }
  }
  //we dont add the intra energy yet - we have to see if there was jump between A and H
  if(iflag_G1!=0 && iflag_G2!=0) // we dont check PRE bc it didnt move
    return iflag_G1;
  
  if((glp_is_flippable[nt_c]==GLP_FIXED || glp_is_flippable[nt_c]==GLP_PUCK) && iflag_G1!=0) return iflag_G1;
  /* bonded interactions */
  energ+=MC_calc_bonded_energy(nt_c, rx, ry, rz, &tflag, iflag_G1, iflag_G2, EintraG1, EintraG2, iflag_G1_PRE, iflag_G2_PRE, EintraG1_PRE, EintraG2_PRE, &gly_flag, &gly_flag_PRE, &TEMP_ERMSD_SQ, &TEMP_ERMSD_ENERG);
  if(tflag!=0)
    return tflag;
  
  //note that the gly_flag has already been updated, and the bad cases have been filtered
  
  //this must happen uniquely and exclusively if the energy has been taken from the "different" glyc state
  //since it is the minimum energy state, 
#ifndef WARMUP
  if(gly_flag!=mc_temp_glyc[nt_c] && (glp_is_flippable[nt_c]==GLP_BOTH || glp_is_flippable[nt_c]==GLP_GLYC) && gly_flag!=-1) { 
    mc_temp_glyc[nt_c]=gly_flag;
    energ+=EintraG2; //now, we add the correct intra energy
  }
  else
#endif
    energ+=EintraG1; // if not, we add the one of the original configuration
  
  if(nt_c>0){  
#ifndef WARMUP
    if(fr_is_mobile[nt_pre]==FR_MOB_FULL && gly_flag_PRE!=mc_temp_glyc[nt_pre] && (glp_is_flippable[nt_pre]==GLP_BOTH || glp_is_flippable[nt_pre]==GLP_GLYC) && gly_flag_PRE!=-1) { 
      //printf("MODIFYING NEIGHBOR GLYC %d  from %d to %d\n", nt_c-1, mc_temp_glyc[nt_pre], gly_flag_PRE);
      mc_temp_glyc[nt_pre]=gly_flag_PRE;
      energ+=EintraG2_PRE; //now, we add the correct intra energy
    }
    else
#endif
      energ+=EintraG1_PRE; // if not, we add the one of the original configuration
  }
#ifdef LNKRMV
  //double phenerg=calc_phpull_energ(nt_c,get_unf_coo_temp_x(at_c+IPHO),get_unf_coo_temp_y(at_c+IPHO),get_unf_coo_temp_z(at_c+IPHO));
  double phenerg=0;
  //if(my_link[nt_c]>-1){
  if(in_link[nt_c]>-1){
    phenerg=calc_link_energy(nt_c, rx, ry, rz);
    energ+=phenerg;
    //printf("%d  %lf\n", nt_c, phenerg);
  }
#endif
  /* non bonded loop */
  //printf("beg\n");
  for(n=0;n<vl_n_pairs[nt_c];n++){
    nt2=vl_neighbor_tables[nt_c][n];
#ifdef FROZEN
    if(fr_is_mobile[nt_c]!=FR_MOB_FROZ || fr_is_mobile[nt2]!=FR_MOB_FROZ)
#endif

      {
	at_ne=N_PARTS_PER_NT*nt2;
	at_c=N_PARTS_PER_NT*nt_c;
#ifdef LNKRMV
	if(in_link[nt_c]<0 || in_link[nt2]<0)
	  if(nts_in_same_link_but_different_loops(nt_c, nt2)==0)
#endif
#ifdef EBUILD
	    //this is for making nts inside ermsd group invisible. For building models from the scratch
	    if(G_groups[nt_c][nt2]<0)
#endif
	    {
	      centdistsq=calc_min_dist_sq(rx[at_ne+ISUG], ry[at_ne+ISUG], rz[at_ne+ISUG], mc_temp_x[at_c+ISUG], mc_temp_y[at_c+ISUG], mc_temp_z[at_c+ISUG]);
	      if(centdistsq<mc_nb_rcut_sq){
		calc_min_vec(rx[at_ne], ry[at_ne], rz[at_ne], mc_temp_x[at_c], mc_temp_y[at_c], mc_temp_z[at_c], r_vec, &r);
		tnb_energ=MC_calc_non_bonded_energy(nt_c, rx, ry, rz, nt2, r_vec, r, &flag);
		energ+=tnb_energ;
		/* if(nt_c==1){ */
		/* printf("%d %lf\n", nt2, tnb_energ); */
		/* } */
		if(flag!=0)
		  return flag;
	      }
	    }
#ifdef ERMSDR
	if(G_groups[nt_c][nt2]>-1){
	  calc_min_vec(rx[at_ne], ry[at_ne], rz[at_ne], mc_temp_x[at_c], mc_temp_y[at_c], mc_temp_z[at_c], r_vec, &r);
	  //if(r<ERMSD_CUTOFF*ERMSDX){
	  proj_on_nt(r_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, pr_vec);
	  proj_on_nt_inv(r_vec, rx, ry,rz, nt2, pr_vec_inv);
	  termsd_p=MC_get_pair_ermsd(pr_vec[0]/ERMSDX    , pr_vec[1]/ERMSDY    , pr_vec[2]/ERMSDZ    , G_ref[nt_c][nt2][0], G_ref[nt_c][nt2][1], G_ref[nt_c][nt2][2], G_ref[nt_c][nt2][3]);
	    termsd_q=MC_get_pair_ermsd(pr_vec_inv[0]/ERMSDX, pr_vec_inv[1]/ERMSDY, pr_vec_inv[2]/ERMSDZ, G_ref[nt2][nt_c][0], G_ref[nt2][nt_c][1], G_ref[nt2][nt_c][2], G_ref[nt2][nt_c][3]);
	    TEMP_ERMSD_SQ+=(termsd_p+termsd_q);
	    TEMP_ERMSD_ENERG+=(0.5*ERMSD_PREF[G_groups[nt_c][nt2]]*(termsd_p+termsd_q)*ERMSD_SSTRUCT[nt_c][nt2]);
	    //}
	  }
#endif  
      }
    }
    *energ_calc=energ;
#ifdef ERMSDR
    *energ_calc+=wall_energy;
    TEMP_ERMSD_SQ/=((double)ERMSD_NNT);
    if(trial==-1){
      DELTA_ERMSD_SQ=-TEMP_ERMSD_SQ;
      DELTA_ERMSD_ENERG=-TEMP_ERMSD_ENERG;
      DELTA_WALL_ENERG=-wall_energy;
      //DELTA_ERMSD_SQ=0;
      //DELTA_ERMSD_ENERG=0;
      //printf("trial is minus 1! %lf  %lf\n", DELTA_ERMSD_SQ, DELTA_ERMSD_ENERG);
    }
    else if(trial==-2){
      //this case is called in the initialization when checking the consistency of gkypuck. it must not affect the values of DELTA_ERMSD_SQ not DELTA_ERMSD_ENERG
      DELTA_ERMSD_SQ=0;
      DELTA_ERMSD_ENERG=0;
      DELTA_WALL_ENERG=0;
      //printf("trial is minus 2! %lf  %lf\n", DELTA_ERMSD_SQ, DELTA_ERMSD_ENERG);
    }
    else{
      DELTA_ERMSD_SQ+=TEMP_ERMSD_SQ;
      DELTA_ERMSD_ENERG+=TEMP_ERMSD_ENERG;
      DELTA_WALL_ENERG+=wall_energy;
      //in this manner, if the step is completely evaluated, we have that DELTA ERMSD = ERMSD_trial - ERMSD_curr. Once this is calculated, we add the corresponding difference of energy
      //ermsd_energ=0.5*ERMSD_PREF*DELTA_ERMSD_SQ;
      //*energ_calc+=ermsd_energ;
      *energ_calc+=DELTA_ERMSD_ENERG;
      
    }
#endif
    //printf("energ = %lf phenerg=%lf\n",energ, phenerg);
    return flag;
}

void MC_select_rand_nt(int nt_n, double *rx, double *ry, double *rz,  int *nt){
  int thnt=rand_i(nt_n);
#ifdef FROZEN
  while(fr_is_mobile[thnt]==FR_MOB_FROZ)
    thnt=rand_i(nt_n);
#endif
  *nt=thnt;
  MC_copy_nt(thnt, rx, ry, rz);
}

void MC_copy_nt(int nt_c, double *rx, double *ry, double *rz){
  int i;
  int at_c=N_PARTS_PER_NT*nt_c;
  //printf("selected %d\n", at_c);
  for(i=0;i<N_PARTS_PER_NT;i++){
    mc_temp_x[at_c+i]=rx[at_c+i];
    mc_temp_y[at_c+i]=ry[at_c+i];
    mc_temp_z[at_c+i]=rz[at_c+i];
    //printf("COPY %lf %lf %lf\n", mc_temp_x[at_c+i], mc_temp_y[at_c+i], mc_temp_z[at_c+i]);
#ifdef PBC
    mc_temp_pbox[at_c+i][0]=mc_pbox[at_c+i][0];
    mc_temp_pbox[at_c+i][1]=mc_pbox[at_c+i][1];
    mc_temp_pbox[at_c+i][2]=mc_pbox[at_c+i][2];
#endif
  }
  mc_temp_glyc[nt_c]=mc_glyc[nt_c];
  mc_temp_puck[nt_c]=mc_puck[nt_c];
}


void MC_map_sugar_temp(int nt_c){
  double Xv[DIM], Yv[DIM], Zv[DIM], msug[DIM];
  int d, at_c=nt_c*N_PARTS_PER_NT;
  Xv[0]=get_unf_coo_temp_x(at_c+IX)-get_unf_coo_temp_x(at_c+IBAS);
  Xv[1]=get_unf_coo_temp_y(at_c+IX)-get_unf_coo_temp_y(at_c+IBAS);
  Xv[2]=get_unf_coo_temp_z(at_c+IX)-get_unf_coo_temp_z(at_c+IBAS);
  Yv[0]=get_unf_coo_temp_x(at_c+IY)-get_unf_coo_temp_x(at_c+IBAS);
  Yv[1]=get_unf_coo_temp_y(at_c+IY)-get_unf_coo_temp_y(at_c+IBAS);
  Yv[2]=get_unf_coo_temp_z(at_c+IY)-get_unf_coo_temp_z(at_c+IBAS);
  
  vec_prod(Xv, Yv, Zv);
  
  msug[0]=MAP_SUG_GL_X[mc_types[at_c]][mc_temp_glyc[nt_c]] * Xv[0] + MAP_SUG_GL_Y[mc_types[at_c]][mc_temp_glyc[nt_c]] * Yv[0] + MAP_SUG_GL_Z[mc_types[at_c]][mc_temp_glyc[nt_c]] * Zv[0] ;
  msug[1]=MAP_SUG_GL_X[mc_types[at_c]][mc_temp_glyc[nt_c]] * Xv[1] + MAP_SUG_GL_Y[mc_types[at_c]][mc_temp_glyc[nt_c]] * Yv[1] + MAP_SUG_GL_Z[mc_types[at_c]][mc_temp_glyc[nt_c]] * Zv[1] ;
  msug[2]=MAP_SUG_GL_X[mc_types[at_c]][mc_temp_glyc[nt_c]] * Xv[2] + MAP_SUG_GL_Y[mc_types[at_c]][mc_temp_glyc[nt_c]] * Yv[2] + MAP_SUG_GL_Z[mc_types[at_c]][mc_temp_glyc[nt_c]] * Zv[2] ;
  
  msug[0]+=get_unf_coo_temp_x(at_c+IBAS);
  msug[1]+=get_unf_coo_temp_y(at_c+IBAS);
  msug[2]+=get_unf_coo_temp_z(at_c+IBAS);
  
  mc_temp_x[at_c+ISUG]=msug[0];
  mc_temp_y[at_c+ISUG]=msug[1];
  mc_temp_z[at_c+ISUG]=msug[2];

  //NEED TO DEFINE PBC!
}

void MC_map_sugar(int nt_c, double **rx, double **ry, double **rz){
  double Xv[DIM], Yv[DIM], Zv[DIM], msug[DIM];
  int d, at_c=nt_c*N_PARTS_PER_NT;
  Xv[0]=get_unf_coo_x(*rx, at_c+IX)-get_unf_coo_x(*rx, at_c+IBAS);
  Xv[1]=get_unf_coo_y(*ry, at_c+IX)-get_unf_coo_y(*ry, at_c+IBAS);
  Xv[2]=get_unf_coo_z(*rz, at_c+IX)-get_unf_coo_z(*rz, at_c+IBAS);
  Yv[0]=get_unf_coo_x(*rx, at_c+IY)-get_unf_coo_x(*rx, at_c+IBAS);
  Yv[1]=get_unf_coo_y(*ry, at_c+IY)-get_unf_coo_y(*ry, at_c+IBAS);
  Yv[2]=get_unf_coo_z(*rz, at_c+IY)-get_unf_coo_z(*rz, at_c+IBAS);
  
  vec_prod(Xv, Yv, Zv);
  
  msug[0]=MAP_SUG_GL_X[mc_types[at_c]][mc_glyc[nt_c]] * Xv[0] + MAP_SUG_GL_Y[mc_types[at_c]][mc_glyc[nt_c]] * Yv[0] + MAP_SUG_GL_Z[mc_types[at_c]][mc_glyc[nt_c]] * Zv[0] ;
  msug[1]=MAP_SUG_GL_X[mc_types[at_c]][mc_glyc[nt_c]] * Xv[1] + MAP_SUG_GL_Y[mc_types[at_c]][mc_glyc[nt_c]] * Yv[1] + MAP_SUG_GL_Z[mc_types[at_c]][mc_glyc[nt_c]] * Zv[1] ;
  msug[2]=MAP_SUG_GL_X[mc_types[at_c]][mc_glyc[nt_c]] * Xv[2] + MAP_SUG_GL_Y[mc_types[at_c]][mc_glyc[nt_c]] * Yv[2] + MAP_SUG_GL_Z[mc_types[at_c]][mc_glyc[nt_c]] * Zv[2] ;
  
  
  msug[0]+=get_unf_coo_x(*rx, at_c+IBAS);
  msug[1]+=get_unf_coo_y(*ry, at_c+IBAS);
  msug[2]+=get_unf_coo_z(*rz, at_c+IBAS);
  
  (*rx)[at_c+ISUG]=msug[0];
  (*ry)[at_c+ISUG]=msug[1];
  (*rz)[at_c+ISUG]=msug[2];
  //NEED TO DEFINE PBC!!
}

int is_purine(int nt){
  if(mc_types[nt*N_PARTS_PER_NT]==TYP_ADENINE || mc_types[nt*N_PARTS_PER_NT]==TYP_GUANINE) return 1;
  else return 0;
}

int is_pyrimidine(int nt){
  
  if(mc_types[nt*N_PARTS_PER_NT]==TYP_URACIL || mc_types[nt*N_PARTS_PER_NT]==TYP_CYTOSINE) return 1;
  else return 0;
}
