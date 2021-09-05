#include "iostream"
#include "vector"

const char *dgemm_desc = "Blocked dgemm.";

// What happens when N%b != 0

int get_start(int b_r, int b_c, int b, int n) {
  return (b_r * b) + ((b_c * b) * n);
}

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
void square_dgemm(int n, double *A, double *B, double *C) {
  // insert your code here: implementation of basic matrix multiple
  std::vector<double> row_buf(n);
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

void square_dgemm_blocked(int n, int block_size, double *A, double *B,
                          double *C) {
  // insert your code here
  int Nb = n / block_size;
  // TODO: Single block_buf
  std::vector<double> a_block_buf(block_size * block_size);
  std::vector<double> b_block_buf(block_size * block_size);
  std::vector<double> c_block_buf(block_size * block_size);
  for (int i = 0; i < Nb; i++) {
    for (int j = 0; j < Nb; j++) {

      int cpos = get_start(i, j, block_size, n);
      copy_block(cpos, C, c_block_buf.data(), block_size, n);
      for (int k = 0; k < Nb; k++) {
        // std::cout << i << "\t" << k << "block\n";
        // std::cout << k << "\t" << j << "block\n";

        int apos = get_start(i, k, block_size, n);
        int bpos = get_start(k, j, block_size, n);
        // std::cout << apos << "\t" << bpos << "\n";
        copy_block(apos, A, a_block_buf.data(), block_size, n);
        copy_block(bpos, B, b_block_buf.data(), block_size, n);

        // Basic matrix mul on block
        for (int row_id = 0; row_id < block_size; row_id++) {
          for (int col_id = 0; col_id < block_size; col_id++) {
            int out_idx = vec_idx(row_id, col_id, block_size);
            for (int m = 0; m < block_size; m++) {
              int left_idx = vec_idx(row_id, m, block_size);
              int right_idx = col_iter(col_id, m, block_size);
              // std::cout << a_block_buf[left_idx] << " * "
              //           << b_block_buf[right_idx] << "\n";
              c_block_buf[out_idx] =
                  c_block_buf[out_idx] +
                  (a_block_buf[left_idx] * b_block_buf[right_idx]);
            }
          }
        }
      }
      // Write result back
      write_block(cpos, c_block_buf.data(), C, block_size, n);
    }
  }
}
