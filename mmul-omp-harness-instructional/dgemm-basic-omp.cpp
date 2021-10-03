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

  //   int row_id, col_id;
  // #pragma omp parallel for default(none) \
//     shared(n, A, B, C, std::cout) private(row_id, col_id) collapse(2)
  //   for (row_id = 0; row_id < n; row_id++) {
  //     for (col_id = 0; col_id < n; col_id++) {
  //       int out_idx = n * col_id + row_id;
  //       int thread_id = omp_get_thread_num();
  // #pragma omp critical
  //       { std::cout << thread_id << "\t" << out_idx << "\n"; }
  //       // for (int k = 0; k < n; k++) {

  //       //   // int left_idx = vec_idx(row_id, k, n);
  //       //   // int right_idx = col_iter(col_id, k, n);
  //       //   // C[out_idx] = C[out_idx] + (A[left_idx] * B[right_idx]);
  //       // }
  //     }
  //   }
}
