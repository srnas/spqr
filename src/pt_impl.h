#ifndef PTIMPL_H
#define PTIMPL_H

#include "mc_global.h"
#include "mc_utils.h"

void PT_read_params(int, int, FILE**);
void PT_write(int, double, double, FILE**);
void PT_end(FILE **);

#endif
