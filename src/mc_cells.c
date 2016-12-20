#include "mc_cells.h"

void add_particle_to_mc_cell(int mc_n, int cell_index, int part){
  mc_cells[part]=mc_cells[mc_n+cell_index];
  mc_cells[mc_n+cell_index]=part;
}

void remove_particle_from_mc_cell(int prev){
  mc_cells[prev]=mc_cells[mc_cells[prev]];
}

void MC_update_linked_lists(int mc_n, double *rx, double *ry, double *rz){
  int j, tmpi, cell;
  int ic[DIM], kc[DIM];
  for(ic[0]=0;ic[0]<mc_nc[0];ic[0]++) 
    for(ic[1]=0;ic[1]<mc_nc[1];ic[1]++)
      for(ic[2]=0;ic[2]<mc_nc[2];ic[2]++){
	cell=index(ic,mc_nc);
	j=mc_cells[mc_n+cell];
	if(j != -1){      /* if cell is not empty...*/
	  /* ...we check all the particles but the first... */
	  if(mc_cells[j] != -1){
	    while(mc_cells[j] != -1){
	      tmpi=mc_cells[j];
	      calculate_cell_index(kc, rx[mc_cells[j]], ry[mc_cells[j]], rz[mc_cells[j]]);
	      //mc_particles[mc_cells[j]], kc);
	      //if(index != cell){
	      if(kc[0]!=ic[0] || kc[1]!=ic[1] || kc[2]!=ic[2]){
		if(kc[0]>=mc_nc[0] || kc[1] >= mc_nc[1] || kc[2]>=mc_nc[2])
		  //		   index(kc,mc_nc)>mc_n_linked_cells)
		  printf("\n\nHORROR 2!!!!%d %d\t%d %d %d  of  %d %d %d\t%lf %lf %lf\n\n", index(kc,mc_nc), mc_n_linked_cells,kc[0], kc[1], kc[2], mc_nc[0], mc_nc[1], mc_nc[2], rx[mc_cells[j]], ry[mc_cells[j]], rz[mc_cells[j]]);
		
		remove_particle_from_mc_cell(j);
		//printf("Particle %d   %lf %lf %lf   old cell %d %d %d  new cell %d %d %d\n", tmpi, rx[tmpi], ry[tmpi], rz[tmpi], ic[0], ic[1], ic[2], kc[0], kc[1], kc[2]);
		add_particle_to_mc_cell(mc_n,index(kc,mc_nc), tmpi);
	      }
	      /* if not, jump to the next */
	      else j=mc_cells[j];
	    }
	  }
	  /* ...and at the end, the first */
	  j=mc_cells[mc_n+cell];
	  calculate_cell_index(kc, rx[j], ry[j], rz[j]);
	  //mc_particles[j], kc);
	  //if(index != cell){
	  if(kc[0]!=ic[0] || kc[1]!=ic[1] || kc[2]!=ic[2]){
	    tmpi=mc_cells[j];
	    mc_cells[mc_n+cell]=mc_cells[j];
	    add_particle_to_mc_cell(mc_n,index(kc,mc_nc), j);
	  }
	}
      }

#ifdef VIRTUAL_SITES
  VS_update_linked_lists(mc_n, rx, ry, rz);
#endif
}

void calculate_cell_index(int *kc, double x, double y, double z){
  kc[0]=(int)floor(x*((double)mc_nc[0])/box_l[0]);
  kc[1]=(int)floor(y*((double)mc_nc[1])/box_l[1]);
  kc[2]=(int)floor(z*((double)mc_nc[2])/box_l[2]);
}
