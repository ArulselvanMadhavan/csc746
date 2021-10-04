#include "likwid-stuff.h"
#include <iostream>
#include <omp.h>
const char *dgemm_desc =
    "Basic implementation, OpenMP-enabled, three-loop dgemm.";

int vec_idx(int r, int c, int n) { return (n * c) + r; }

int col_iter(int col_id, int k, int n) { return (n * col_id) + k; }
/*
 * This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in column-major format.
 * On exit, A and B maintain their input values.
 */
void square_dgemm(int n, double *A, double *B, double *C) {
  int i, j;
#pragma omp parallel default(none) shared(n, A, B, C) private(i, j)
  {
#ifdef LIKWID_PERFMON
    LIKWID_MARKER_START(MY_MARKER_REGION_NAME);
#endif
#pragma omp for collapse(2)
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
          C[n * j + i] += A[n * k + i] * B[n * j + k];
        }
      }
    }
#ifdef LIKWID_PERFMON
    LIKWID_MARKER_STOP(MY_MARKER_REGION_NAME);
#endif
  }
}
