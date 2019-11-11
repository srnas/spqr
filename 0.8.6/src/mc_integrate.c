#include "mc_integrate.h"

void MC_integrate(int mc_n, double **rx, double **ry, double **rz){
  int nt_n=mc_n/N_PARTS_PER_NT, i;
  /* int a,b; */
  /* for(a=0;a<nt_n;a++){ */
  /*   printf("%d    %d : ", a, vl_n_pairs[a]); */
  /*   for(b=0;b<vl_n_pairs[a];b++) */
  /*     printf("%d  ",vl_neighbor_tables[a][b]); */
  /*   printf("\n"); */
  /* } */
  //  nt2=vl_neighbor_tables[the_nt][n];

  int nt_c;
  double old_energy, new_energy, wc_energ_diff=0.0; // wc energ is the difference between the new and the old : Ewc - Ewc_old
  int mc_trial_flag=0;
  int nt_to_select=nt_n;
  
#ifdef MCVLISTS
  for(i=0;i<nt_n;i++){
    if(vl_count[i]>=vl_ncrit-1)
      MC_build_verlet_lists(nt_n, *rx, *ry, *rz, i);
  }
  //printf("Updating Verlet list - vlncrit = %d\n", vl_ncrit);
#endif
  //MC_eval_total_energy(nt_n, *rx, *ry, *rz);
  
  MC_select_rand_nt(nt_to_select, *rx, *ry, *rz, &nt_c); //here we copy the position to the temp position (the whole nucleotide). nt_c is between 0 and mc_n/N_PARTS_PER_NT
  //printf("*************************   SELECTED  %d   ****************************\n", nt_c);
  //printf("tot pre\n");
  //double wce_init=MC_calculate_total_wc_energy(nt_n, -1, *rx, *ry, *rz);
  

  mc_trial_flag=MC_calculate_local_energy(*rx, *ry, *rz, nt_c, &old_energy); //here we calculate the energy of the selected particle, from temp position
  if(mc_trial_flag!=0){ printf("THIS SHOULDNT HAPPEN\n flag = %d, selected %d\n", mc_trial_flag, nt_c);
    MC_get_sp_bonded_energy(nt_c, nt_n, *rx, *ry, *rz);
    //MC_get_energy(nt_n, *rx, *ry, *rz,0);MC_get_energy(nt_n, *rx, *ry, *rz,1);MC_get_energy(nt_n, *rx, *ry, *rz,2);
    exit(1);}
  int trial=MC_perform_trial_movement(nt_c); //here we update the position to temp position
  mc_trial_flag=MC_calculate_local_energy(*rx, *ry, *rz, nt_c, &new_energy); //here we calculate the energy of the selected particle, from temp position
  //printf("tot post\n");
  //double wce_trial=MC_calculate_total_wc_energy(nt_n, nt_c, *rx, *ry, *rz);
  
  double dwce=MC_calculate_local_wc_energy(nt_n, nt_c, *rx, *ry, *rz);
  //wc_energ_diff=wce_trial-wce_init;
  /* fflush(stdout); */
  /* if(fabs(wc_energ_diff - dwce) > 0.000001){ */
  /*   printf("ways of calculating energy differ!\nLOCAL=%lf\nTOTAL=%lf\n%lf\n", dwce, wc_energ_diff,fabs( wc_energ_diff-dwce)); */
  /*   exit(1); */
  /* } */
  if(mc_trial_flag==0){
    MC_eval_displacement(nt_n,rx, ry, rz, nt_c, old_energy, new_energy, dwce); //if the movement is accepted, we update the position
    //printf("trial evaluated!  ENERG OLD : %lf\tENERG NEW %lf \t ENERG WC: %lf\n", old_energy, new_energy, wc_energ_diff);
  }
  //else printf("trial rejected (not evaluated)!\n");
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
      //(rx[ba]+mc_pbox[ba][0]*box_l[0]) - (mc_temp_x[0]+mc_temp_pbox[0][0]*box_l[0]);
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
    
    for(n=0;n<vl_n_pairs[nt_c];n++){
      nt2=vl_neighbor_tables[nt_c][n];
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
	//nb_wc =0;// MC_calc_nnN_watscric(typ_ind, r_vec, r_vec_inv, theta, &temp_flag);
	//if(temp_flag==0){
	// if(nb_wc<0)printf(" wcp %d ", nt2);//if(cflag==-1)printf( " %lf ", nb_wc);
	//}
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
    //THE ADDITIONAL PAIR
    nt2=nt_c-1;
    for(nt2=nt_c-1;nt2<=nt_c+1;nt2+=2){
      if(nt2>=0 && nt2<nt_n){
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
      for(nt2=nt_c-1;nt2<=nt_c+1;nt2+=2){
      	if(nt2>=0 && nt2<nt_n){
      	  ba=N_PARTS_PER_NT*nt2;
      	  typ_ind = mc_n_types*mc_types[at_c]  + mc_types[ba];
      	  typ_ind2= mc_n_types*mc_types[ba] + mc_types[at_c];
      	  sugdistsq=calc_min_dist_sq(rx[ba+ISUG], ry[ba+ISUG], rz[ba+ISUG], mc_temp_x[at_c+ISUG], mc_temp_y[at_c+ISUG], mc_temp_z[at_c+ISUG]);
      	  if(sugdistsq<mc_r_cut_sq){
      	    calc_min_vec(rx[ba], ry[ba], rz[ba], mc_temp_x[at_c], mc_temp_y[at_c], mc_temp_z[at_c], d_vec, &dist);
      	    proj_on_nt(d_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
      	    proj_on_nt_inv(d_vec, rx, ry, rz, nt2, r_vec_inv);
	    
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
      		fprintf(mcconf,"S");
      	      else if(temp_face2==WC_FACE_WATSCRICK)
      		fprintf(mcconf,"W");
      	      else if(temp_face2==WC_FACE_HOOGSTEEN)
      		fprintf(mcconf,"H");
      	      fprintf(mcconf,"P)");
      	    }
      	    temp_face1=wc_pair_array1[2];
      	    temp_face2=wc_pair_array2[1];
      	    if(temp_face1>=0 && temp_face2>=0) {
      	      fprintf(mcconf," bph %d (P", nt2);
      	      if(temp_face1==WC_FACE_SUGAR)
      		fprintf(mcconf,"S");
      	      else if(temp_face1==WC_FACE_WATSCRICK)
      		fprintf(mcconf,"W");
      	      else if(temp_face1==WC_FACE_HOOGSTEEN)
      		fprintf(mcconf,"H");
      	      fprintf(mcconf,")");
      	    }
      	  }
      	}
      }
      fprintf(mcconf,"\n");
    }
    fclose(mcconf);
  }
}

void MC_get_sp_bonded_energy(int nt_c, int nt_n, double *rx, double *ry, double *rz){
  //central
  double energ_c=0, energ_pr=0, energ_po=0;
  MC_copy_nt(nt_c, rx, ry, rz);
  int tflag=0;
  energ_c+=MC_calc_bonded_energy(nt_c, rx, ry, rz, &tflag)/2.0;
  //prev
  if(nt_c>0){
    mc_nbonds[nt_c-1][0]=1;
    mc_bondlist[nt_c-1][0]=nt_c;
    MC_copy_nt(nt_c-1, rx, ry, rz);
    tflag=0;
    energ_pr+=MC_calc_bonded_energy(nt_c-1, rx, ry, rz, &tflag)/2.0;
  }
  //post
  if(nt_c<nt_n){
    mc_nbonds[nt_c+1][0]=1;
    mc_bondlist[nt_c+1][0]=nt_c;
    MC_copy_nt(nt_c+1, rx, ry, rz);
    tflag=0;
    energ_po+=MC_calc_bonded_energy(nt_c+1, rx, ry, rz, &tflag)/2.0;
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
  int tflag=0;
  int b, nt_neigh, ba;
  double eta, theta;
  int temp_flag;
  double sigma, r;
  int n, nt2;
  double b_ev_temp1, b_ev_temp2, t_vec[DIM], nb_ev_temp1, nb_ev_temp2;
  //int a;//,b;
  /* for(a=0;a<nt_n;a++){ */
  /*   printf("%d    %d : ", a, vl_n_pairs[a]); */
  /*   for(b=0;b<vl_n_pairs[a];b++) */
  /*     printf("%d  ",vl_neighbor_tables[a][b]); */
  /*   printf("\n"); */
  /* } */
  //printf("Evaluating energy\n");
  double energ_c=0;
  for(nt_c=0;nt_c<nt_n;nt_c++){
    //printf("Part %d (of %d)\n", nt_c, nt_n);
    at_c=N_PARTS_PER_NT*nt_c;
    MC_copy_nt(nt_c, rx, ry, rz);
    tflag=0;
    
    /* self interaction (or intra-nt) */
    if(index==0 || index ==3)
      energ_c+=MC_calc_intra_energy(nt_c,&tflag);
    
    if(tflag!=0 && index==0)
      {  printf("Intra: flag = %d\n", tflag);exit(0);	}//   return;	}
    /* bonded interactions */
    
     //tflag=0;
     if(index==1 || index==3)
       energ_c+=MC_calc_bonded_energy(nt_c, rx, ry, rz, &tflag)/2.0;
     if(tflag!=0 && index==1)
       {  printf("Bonded flag = %d\n", tflag);exit(0);	}
     
     
     /* non bonded loop */
     
     //printf("i survived the bonded loop\n");
     if(index==2 || index==3){
       for(n=0;n<vl_n_pairs[nt_c];n++) {
	 nt2=vl_neighbor_tables[nt_c][n];
	 at_ne=N_PARTS_PER_NT*nt2;
	  sugdistsq=calc_min_dist_sq(rx[at_ne+ISUG], ry[at_ne+ISUG], rz[at_ne+ISUG], mc_temp_x[at_c+ISUG], mc_temp_y[at_c+ISUG], mc_temp_z[at_c+ISUG]);                                     
	  if(sugdistsq<mc_r_cut_sq) {
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
   //ADD WATSON-CRICK!
  if(index==3 || index==4) energ_c+=MC_calculate_total_wc_energy(nt_n,-1,rx, ry, rz);
   
   return energ_c;
}




void MC_update_positions(double **rx, double **ry, double **rz, int nt){
  int i;
  //printf("from %lf %lf %lf\to to %lf %lf %lf\n", (*rx)[N_PARTS_PER_NT*nt], (*ry)[N_PARTS_PER_NT*nt],(*rz)[N_PARTS_PER_NT*nt], mc_temp_x[nt_c], mc_temp_y[0],mc_temp_z[0]);
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
  mc_glyc[nt]=mc_temp_glyc[nt];
  mc_puck[nt]=mc_temp_puck[nt];
}

void MC_eval_displacement(int nt_n, double **rx, double **ry, double **rz, int nt, double o_ene, double n_ene, double wc_ene_nmo){
  double ranf=rand_d(1.0);
  if(ranf <= exp(-(n_ene-o_ene+wc_ene_nmo)/mc_target_temp )){
    MC_update_positions(rx, ry, rz, nt);
    //printf("trial accepted!\t %lf %lf     %lf\t %lf %lf  %lf\n", o_ene, n_ene, wc_ene_nmo, ranf,exp(-(n_ene-o_ene+wc_ene_nmo)/mc_target_temp ), mc_target_temp );
    MC_add_vl_count(nt);
    MC_update_wc_lists(nt_n);
  }
  //else printf("trial rejected!\t %lf %lf     %lf\t%lf %lf  %lf\n", o_ene, n_ene, wc_ene_nmo,ranf,exp(-(n_ene-o_ene+wc_ene_nmo)/mc_target_temp ), mc_target_temp);
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
  int rand1=rand_i(6);
#ifdef FROZEN
  if(fr_is_mobile[nt_c]==FR_MOB_BASE) //only NT is mobile - the flag has to be > 0 
    rand1=rand_i(5)+1;
  else if(fr_is_mobile[nt_c]==FR_MOB_PHOS)
    rand1=0;
#endif
  if(N_PARTS_PER_NT>1){
    if(rand1==0){
      MC_perform_nt_rotation(nt_c,0); // ROTATION AROUND SUGAR
      //printf(" NT rot       \n");
    }
    else if(rand1==1){
      MC_perform_nt_rotation(nt_c,1); // ROTATION AROUND BASE
      //printf(" NT rot       \n");
    }
    else if(rand1==2){
      MC_perform_nt_rotation(nt_c,2); // ROTATION AROUND (SUG+BAS)/2
      //printf(" NT rot       \n");
    }
    else if(rand1==3){
        MC_perform_p_translation(nt_c);
    }
    //else if(rand1==4){
    else{
      MC_perform_nt_translation(nt_c);
    }

    
    /* else if (rand1==3){ */
    /*   MC_perform_nt_rotation(nt_c,1); */
    /* } */
      /* else */
    /* 	{ */
    /* 	  MC_perform_nt_rotation(1); // ROTATION AROUND BASE CENTER */
    /* 	  //printf("NT rot 1 "); */
    /* 	} */
    /* #ifdef GLYC_TRIAL */
    /* else if(rand1==3){ */
    /*   MC_perform_base_flip(nt_c); */
    /* } */
    /* #endif */
    /* else{ */
    /*   MC_perform_bb_rotation(nt_c); */
    /*   //printf(" BB rot        \n"); */
    /* } */
  }
  else
    MC_perform_nt_translation(nt_c);
  return rand1;
}

/* void MC_perform_base_flip(int nt_c){ */
/*   // we displace the base */
/*   //then flip it */
/*   //then update the glycosidic index */
/*   //then remap the sugar */
  
/* } */

int MC_calculate_local_energy(double *rx, double *ry, double *rz, int nt_c, double *energ_calc){
  //NO CELL STRUCTURE IMPLEMENTED YET
   /* for three dimensions! */
  int n, nt2, at_c, at_ne;
  double r, r_vec[DIM], centdistsq;
  //HERE WE HAVE TO CALCULATE THE ENERGIES OF THE AUXILIAR PARTICLES TOO!!
  double energ=0.0;
  int flag=0;
  
  /* self interaction (or intra-nt) */
  energ+=MC_calc_intra_energy(nt_c,&flag);
  if(flag!=0)
    return flag;
  /* bonded interactions */
  energ+=MC_calc_bonded_energy(nt_c, rx, ry, rz, &flag);
  if(flag!=0)
    return flag;
  /* non bonded loop */
  //printf("i survived the bonded loop\n");
  for(n=0;n<vl_n_pairs[nt_c];n++){
    nt2=vl_neighbor_tables[nt_c][n];
    at_ne=N_PARTS_PER_NT*nt2;
    at_c=N_PARTS_PER_NT*nt_c;
    centdistsq=calc_min_dist_sq(rx[at_ne+ISUG], ry[at_ne+ISUG], rz[at_ne+ISUG], mc_temp_x[at_c+ISUG], mc_temp_y[at_c+ISUG], mc_temp_z[at_c+ISUG]);
    
    
    //centdistsq=calc_min_dist_sq(0.5*(rx[at_ne+IPHO]+rx[at_ne]), 0.5*(ry[at_ne+IPHO]+ry[at_ne]), 0.5*(rz[at_ne+IPHO]+rz[at_ne]), 				0.5*(mc_temp_x[at_c+IPHO] + mc_temp_x[at_c]), 0.5*(mc_temp_y[at_c+IPHO]+mc_temp_y[at_c]), 0.5*(mc_temp_z[at_c+IPHO]+mc_temp_z[at_c]));
    
    if(centdistsq<mc_nb_rcut_sq){
      calc_min_vec(rx[at_ne], ry[at_ne], rz[at_ne], mc_temp_x[at_c], mc_temp_y[at_c], mc_temp_z[at_c], r_vec, &r);
      energ+=MC_calc_non_bonded_energy(nt_c, rx, ry, rz, nt2, r_vec, r, &flag);
      if(flag!=0)
  	return flag;
    }
  }
  //HERE WE HAVE TO CALCULATE THE ENERGIES OF THE AUXILIAR PARTICLES TOO!!
  *energ_calc=energ;
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
