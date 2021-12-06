#include "malloc.h"
#include "omp.h"
#include "stdio.h"
#include "timer.h"

#define N (100000000)

void test() {

  /* float *a; */
  /* float *b; */
  float *a = malloc(sizeof(float) * N);
  float *b = malloc(sizeof(float) * N);
  float *c = malloc(sizeof(float) * 5);
  /* long long int idx; */
  struct timespec starttime;
  cpu_timer_start(&starttime);

#pragma omp target data map(alloc : a [0:N - 1], b [0:N - 1], c [0:N - 1])
  {
#pragma omp target teams distribute parallel for map(always, from:c[0:5])
    for (long long int i = 0; i < N; i++) {
      a[i] = b[i] = 2.0;
      c[i] = a[i] * b[i];
    }
  }

  double totaltime = cpu_timer_stop(starttime);
  printf("%f\t%f\n", c[3], totaltime);
}
