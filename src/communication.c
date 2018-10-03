#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "constant.h"
#include "prototypes.h"
#include "morton.h"

#ifdef WMPI
// =================================================
void init_MPI(struct CPU *cpu){
  MPI_Comm_size(MPI_COMM_WORLD,&(cpu->nproc));
  MPI_Comm_rank(MPI_COMM_WORLD,&(cpu->rank));
}
// =================================================
#endif WMPI
