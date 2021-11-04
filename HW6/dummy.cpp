#include "mpi.h"

int main(int argc, char* argv[]){
  MPI_Init(&argc, &argv);
  int rank = 0, nprocs = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  printf("Rank:%d\tNprocs:%d\n", rank, nprocs);
  MPI_Finalize();
  return 0;
}
