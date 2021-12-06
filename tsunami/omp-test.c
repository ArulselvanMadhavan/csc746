#include "malloc.h"
#include "memory.h"
#include "omp.h"
#include "stdio.h"
#include "timer.h"

#define N (1024)

void test() {

  double *restrict a = malloc(sizeof(double) * N * N);
  double *restrict b = malloc(sizeof(double) * N * N);
  double *restrict c = malloc(sizeof(double) * N * N);
  struct timespec starttime;
  cpu_timer_start(&starttime);

  // clang-format off
#pragma omp target data map(alloc : a[:N*N], b [:N*N], c[:N*N])
  // clang-format on
  {
#pragma omp target teams distribute parallel for map(always, from : c [0:5])
    for (long long int i = 0; i < N*N; i++) {
      a[i] = b[i] = 2.0;
      c[i] = a[i] * b[i];
    }
  }

  double totaltime = cpu_timer_stop(starttime);
  printf("Print:%f\t%f\n", c[3], totaltime);
}
