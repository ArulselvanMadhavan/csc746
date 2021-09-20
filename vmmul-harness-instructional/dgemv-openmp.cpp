#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const char *dgemv_desc = "OpenMP dgemv.";

/*
 * This routine performs a dgemv operation
 * Y :=  A * X + Y
 * where A is n-by-n matrix stored in row-major format, and X and Y are n by 1
 * vectors. On exit, A and X maintain their input values.
 */

void my_dgemv(int n, double *A, double *x, double *y) {
  int i, j;
#pragma omp parallel for default(none) shared(n, A, x, y) private(i, j)
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      y[i] += A[i * n + j] * x[j];
    }
  }
}
