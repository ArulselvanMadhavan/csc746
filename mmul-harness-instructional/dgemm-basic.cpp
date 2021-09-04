#include <iostream>
using std::cout;

const char *dgemm_desc = "Basic implementation, three-loop dgemm.";

int vec_idx(int r, int c, int n) { return (n * c) + r; }

int col_iter(int col_id, int k, int n) { return (n * col_id) + k; }

/*
 * This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in column-major format.
 * On exit, A and B maintain their input values.
 */
void square_dgemm(int n, double *A, double *B, double *C) {
  // insert your code here: implementation of basic matrix multiple
  for (int row_id = 0; row_id < n; row_id++) {
    for (int col_id = 0; col_id < n; col_id++) {
      int out_idx = vec_idx(row_id, col_id, n);      
      for (int k = 0; k < n; k++) {
        int left_idx = vec_idx(row_id, k, n);
        int right_idx = col_iter(col_id, k, n);
        C[out_idx] = C[out_idx] + (A[left_idx] * B[right_idx]);
      }
    }
  }
}
