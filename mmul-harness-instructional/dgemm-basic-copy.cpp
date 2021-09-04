#include "vector"
#include <iostream>

using std::cout;

const char *dgemm_desc = "Basic implementation Copy opt, three-loop dgemm.";

int vec_idx(int r, int c, int n) { return (n * c) + r; }

int col_iter(int col_id, int k, int n) { return (n * col_id) + k; }

void copy_row(double *src, double *dest, int r, int n) {
  for (int c = 0; c < n; c++) {
    dest[c] = src[vec_idx(r, c, n)];
  }
}

void copy_col(double *src, double *dest, int c, int n) {
  for (int k = 0; k < n; k++) {
    dest[k] = src[col_iter(c, k, n)];
  }
}
/*
 * This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in column-major format.
 * On exit, A and B maintain their input values.
 */
void square_dgemm(int n, double *A, double *B, double *C) {
  // insert your code here: implementation of basic matrix multiple
  std::vector<double> row_buf(n);
  std::vector<double> col_buf(n);
  for (int row_id = 0; row_id < n; row_id++) {
    copy_row(A, row_buf.data(), row_id, n);
    for (int col_id = 0; col_id < n; col_id++) {
      double c_out = 0;
      for (int k = 0; k < n; k++) {
        c_out = c_out + (row_buf[k] * B[(n * col_id) + k]);
      }
      C[vec_idx(row_id, col_id, n)] = c_out;
    }
  }
}