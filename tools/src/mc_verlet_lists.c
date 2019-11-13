#include "mc_verlet_lists.h"

int offset_array[N_OFFSET][DIM]= {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{-1,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1},{-1,1,1},{-1,0,1},{-1,-1,1},{0,-1,1},{1,-1,1}};

void MC_initialize_verlet_lists(int nt_n, double *rx, double *ry, double *rz){
  int i;
  //int nt_n=mc_n/N_PARTS_PER_NT;
  /* EVERYTHING IS DEFINED WITH RESPECT TO THE NUMBER OF NUCLEOTIDES */
  double maxdist=MC_NT_XYZ;
  if(MC_PH_XYZ > maxdist)
    maxdist=MC_PH_XYZ;
  vl_rmax=mc_r_cut+vl_skin;
  vl_rmax_sq=vl_rmax*vl_rmax;
  vl_max_displ = 2*sqrt(3.0*maxdist);
  vl_ncrit = (int) (vl_skin/vl_max_displ);
  
  vl_n_pairs=(int*)malloc(sizeof(int)*nt_n);
  vl_count=(int *)malloc(sizeof(int)*nt_n);
  //vl_disp=(double *)malloc(sizeof(double)*nt_n);
  
  for(i=0;i<nt_n;i++){
    vl_n_pairs[i]=0;
    vl_count[i]=vl_ncrit-1;
    //vl_disp[i]=0.0;
  }
  
  vl_neighbors_per_atom=nt_n;
  //#ifdef MCVLISTS
  // vl_neighbors_per_atom=6;
  //#endif
  
  //vl_max_vel=0.0;
  vl_n_pairs_max=nt_n;
  
  //vl_neighbor_table=(int*)malloc(sizeof(int)*2*vl_n_pairs_max);
  vl_neighbor_tables=(int **)malloc(sizeof(int*)*nt_n);
  for(i=0;i<nt_n;i++){
    vl_neighbor_tables[i]=(int *)malloc(sizeof(int)*vl_neighbors_per_atom);
  }
  //if(vl_rmax>mc_linked_cell_l){
  //fprintf(stderr, "Side of the linked cells %f incompatible with skin of neighbor list %f.\n", mc_linked_cell_l, vl_rmax);
  //exit(ERR_NEIGHBOR);
  //}
#ifdef MCVLISTS
  for(i=0;i<nt_n;i++){
    MC_build_verlet_lists(nt_n, rx, ry, rz, i);
  }
#else
  int cvl,j;
  for(i=0;i<nt_n;i++){
    cvl=0;
    for(j=0;j<nt_n;j++){
      if(j!=i && !MC_are_neighbors(i,j)){
	//if(j!=i){
	vl_neighbor_tables[i][cvl]=j;
	cvl++;
      }
    }
    vl_n_pairs[i]=cvl;
  }
#endif
}

void MC_build_verlet_lists(int nt_n, double *rx, double *ry, double *rz, int np){
  int i,j,p,jind, ni;
  //mci, mcj
  //int nt_n=mc_n/N_PARTS_PER_NT;
  //int offset;
  //int ic[DIM], jc[DIM];
  double r, r_vec[DIM];
  
  //#ifndef MCDILUTE
  //MC_update_linked_lists(nt_n, rx, ry, rz);
  //#endif
  vl_n_pairs[np]=0;
  if(vl_rmax>0){
#ifdef MCDILUTE
    //LOOP IS OVER NUCLEOTIDES!! VERY IMPORTANT!!!
    p=np*N_PARTS_PER_NT;
    for(ni=0;ni<nt_n;ni++){
      //for(nj=i+1;nj<nt_n;nj++){
#endif
      i=ni*N_PARTS_PER_NT;
      //  j=nj*N_PARTS_PER_NT;
      //IF the nonbonded interaction between such pairs exists...
      if(i!=p){
	//if(mc_lj_eps[mc_types[i]][mc_types[p]]!=0 && i!=p){
	calc_min_vec(rx[i], ry[i], rz[i], rx[p], ry[p], rz[p], r_vec, &r);
#ifdef EXCLUSIONS
	if(!MC_are_neighbors(ni,np))
#endif
	  {
#ifdef MCVLISTS
	    if(r<=vl_rmax)
#endif		
	      {
		//printf("%d %d   %lf\n", i,j,r);
		if(vl_n_pairs[np]>=vl_n_pairs_max || vl_n_pairs[ni]>=vl_n_pairs_max){
		  fprintf(stderr, "Too many neighbors (%d , %d) for neighbor list (max %d) of nucleotides %d and %d!\n", vl_n_pairs[np], vl_n_pairs[ni], vl_n_pairs_max, np, ni);
		  exit(ERR_NEIGHBOR);
		}
		vl_neighbor_tables[np][vl_n_pairs[np]]=ni;
		vl_n_pairs[np]++;
		
		jind=0;
		for(j=0;j<vl_n_pairs[ni];j++){
		  if(vl_neighbor_tables[ni][j]==np)
		    jind=1;
		}
		if(jind==0){
		  vl_neighbor_tables[ni][vl_n_pairs[ni]]=np;
		  vl_n_pairs[ni]++;
		}
	      }
	  }
      }
#ifdef MCDILUTE
      //}
    }
#endif 
  }
  
  //vl_max_vel=0.0;
  vl_count[np]=0;
}

void MC_add_vl_count(int np){
  vl_count[np]++;
}

void MC_free_verlet_lists(){
  free(vl_neighbor_tables);

}

#ifdef EXCLUSIONS
int MC_are_neighbors(int i, int j){
  int b;
  int n=0;
  for(b=0;b<mc_nbonds[i][0];b++){
    if(j==mc_bondlist[i][b])
      n=1;
  }
  for(b=0;b<mc_nbonds[j][0];b++){
    if(i==mc_bondlist[j][b])
      n=1;
  }
  return n;
}
#endif
