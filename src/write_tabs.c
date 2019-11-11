#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mc.h"
#include "mc_energies.h"
//#include "mc_base_pairing.h"
//#include "mc_checkpoints.h"

int main(){
  MC_read_write_energy_tables();
  //MC_read_bin_energy_tables();
  return 0;
}
