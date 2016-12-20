#include "mc_checkpoints.h"

void MC_initialize_save_configs(int mc_n, int mc_ini, int mpi_id){
  int i;
  char confname[256];
  if(mc_ini==0){
    sprintf(confname, "configs/confs.p%02d.mc", mpi_id);
    printf("Saving configurations to file %s.\n", confname);
  }
  else {
    sprintf(confname, "configs/confs.%08d.p%02d.mc", mc_ini, mpi_id);
  }
  if((mc_configs=fopen(confname,"ab"))!=NULL){
    //the simulation begins! write the basic data
    fwrite(&mc_n, sizeof(int), 1, mc_configs);
    for(i=0;i<mc_n;i++){
      fwrite(&(mc_types[i]), sizeof(int), 1, mc_configs);
    }
  }
}

void MC_close_configs(){
  fclose(mc_configs);
}

void MC_save_configuration(int mc_n, double *rx, double *ry, double *rz){
  int i;
  double tempf;
  
  for(i=0;i<mc_n;i++){
    tempf=(double)(get_unf_coo_x(rx,i));
      //rx[i]+mc_pbox[i][0]*box_l[0]);
    fwrite(&tempf,sizeof(double),1,mc_configs);
    tempf=(double)(get_unf_coo_y(ry,i));
    fwrite(&tempf,sizeof(double),1,mc_configs);
    tempf=(double)(get_unf_coo_z(rz,i));
    fwrite(&tempf,sizeof(double),1,mc_configs);
  }
}

void MC_save_xyz_configuration(int mc_n, double *rx, double *ry, double *rz, double final_iter, int mpi_id){
  FILE * mcconf;
  char filename[256];
  sprintf(filename, "configs/lastmc.p%02d.xyz", mpi_id);
  
  if(!(mcconf=fopen(filename , "w"))){
    fprintf(stderr, "Unable to write final MC configuration file.\n");
    //exit(ERR_INPUT);
  } else {
    int i;
    //double unfolded_pos[DIM];
    fprintf(mcconf, "%d\n", mc_n);
    fprintf(mcconf, "Final configuration at step %f\n", final_iter);
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
}

void MC_save_current_configuration(int mc_n, double *rx, double *ry, double *rz, double iter, int index, int mpi_id){
  FILE * mcconf;
  char filename[256];
  sprintf(filename, "configs/cmc.%08d.p%02d.xyz",index, mpi_id);
  //  char filename[11]="lastmc.xyz";
  
  if(!(mcconf=fopen(filename , "w"))){
    fprintf(stderr, "Unable to write final MC configuration file.\n");
    //exit(ERR_INPUT);
  } else {
    int i;
    //double unfolded_pos[DIM];
    fprintf(mcconf, "%d\n", mc_n);
    fprintf(mcconf, "Configuration at step %f\n", iter);
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
}


void MC_min_energ_xyz_configuration(int mc_n, double *rx, double *ry, double *rz, double energ, int init){
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
    fprintf(mcconf, "Energy %lf\n", energ);
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
