#include "stdio.h"
#include "stdlib.h"

const int ny = 200;
const int nx = 500;
const int nhalo = 1;

double **malloc2D(int jmax, int imax);

int main(int argc, char *argv[]) {
  printf("Running Shallow Water equations");
  int gdims[] = {ny + 2 * nhalo, nx + 2 * nhalo};
  double **restrict H = malloc2D(gdims[0], gdims[1]);
  double **restrict U = malloc2D(gdims[0], gdims[1]);
  double **restrict V = malloc2D(gdims[0], gdims[1]);
  free(H);
  free(U);
  free(V);
  return 0;
}
