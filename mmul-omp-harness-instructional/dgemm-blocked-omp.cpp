#include "vector"
#include <iostream>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
// #include <vector.h>
#include "likwid-stuff.h"

const char *dgemm_desc = "Blocked dgemm, OpenMP-enabled";

int get_start(int b_r, int b_c, int b, int n) {
  return (b_r * b) + ((b_c * b) * n);
}

void copy_block(int start, double *src, double *dest, int b, int n) {
  for (int c = 0; c < b; c++) {
    for (int r = 0; r < b; r++) {
      dest[(c * b) + r] = src[start + (c * n) + r];
    }
  }
}

void write_block(int start, double *src, double *dest, int b, int n) {
  for (int c = 0; c < b; c++) {
    for (int r = 0; r < b; r++) {
      dest[start + (c * n) + r] = src[(c * b) + r];
    }
  }
}

int vec_idx(int r, int c, int n) { return (n * c) + r; }

int col_iter(int col_id, int k, int n) { return (n * col_id) + k; }

void square_dgemm_blocked(int n, int block_size, double *A, double *B,
                          double *C) {
  int Nb = n / block_size;
  int total_threads = omp_get_num_threads();
  int bi, bj;
#pragma omp parallel default(none) shared(                                     \
    n, block_size, Nb, A, B, C, total_threads, std::cout) private(bi, bj)
  {
#ifdef LIKWID_PERFMON
    LIKWID_MARKER_START(MY_MARKER_REGION_NAME);
#endif

    std::vector<double> block_buf(3 * block_size * block_size);
    double *AA = block_buf.data() + 0;
    double *BB = AA + block_size * block_size;
    double *CC = BB + block_size * block_size;

#pragma omp for collapse(2)
    for (bi = 0; bi < Nb; bi++) {
      for (bj = 0; bj < Nb; bj++) {
        int cpos = get_start(bi, bj, block_size, n);
        copy_block(cpos, C, CC, block_size, n);
        for (int k = 0; k < Nb; k++) {
          int apos = get_start(bi, k, block_size, n);
          int bpos = get_start(k, bj, block_size, n);

          copy_block(apos, A, AA, block_size, n); // Nb^3
          copy_block(bpos, B, BB, block_size, n); // Nb^3
                                                  // Basic matrix mul on block
          for (int row_id = 0; row_id < block_size; row_id++) {
            for (int col_id = 0; col_id < block_size; col_id++) {
              int out_idx = vec_idx(row_id, col_id, block_size); // 2*Nb^3*b^2
              double c_cout = CC[out_idx];
              for (int m = 0; m < block_size; m++) {
                int left_idx = vec_idx(row_id, m, block_size);    // 2*Nb^3*b^3
                int right_idx = col_iter(col_id, m, block_size);  // 2*Nb^3*b^3
                c_cout = c_cout + (AA[left_idx] * BB[right_idx]); // 2*Nb^3*b^3
              }
              CC[out_idx] = c_cout;
            }
          }
        }
        write_block(cpos, CC, C, block_size, n);
#pragma omp critical
        {
          std::cout << omp_get_thread_num() << "\t" << AA << "\t" << BB << "\n";
        }
      }
    }
#ifdef LIKWID_PERFMON
    LIKWID_MARKER_STOP(MY_MARKER_REGION_NAME);
#endif
  }
}
