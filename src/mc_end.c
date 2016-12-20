#include "mc_end.h"

void MC_end(int mc_n, double *rx, double *ry, double *rz, double final_time, int mpi_id){
  MC_save_xyz_configuration(mc_n, rx, ry, rz, final_time, mpi_id);
}
