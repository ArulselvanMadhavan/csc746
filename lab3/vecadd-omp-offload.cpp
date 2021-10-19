#include <chrono>
#include <iostream>
#include <math.h>

void add(int N, float *X, float *Y) {
#pragma omp target data map(to : X [0:N], N) map(tofrom : Y [0:N])
  {
#pragma omp target teams distribute parallel for
    for (int i = 0; i < N; i++) {
      Y[i] = X[i] + Y[i];
    }
  }
}

int main(void) {
  const int N = 1 << 28;

  float *X = new float[N];
  float *Y = new float[N];

  // Initialize
  for (int i = 0; i < N; i++) {
    X[i] = 1.0f;
    Y[i] = 2.0f;
  }

  std::chrono::time_point<std::chrono::high_resolution_clock> start_time =
      std::chrono::high_resolution_clock::now();
  add(N, X, Y);
  std::chrono::time_point<std::chrono::high_resolution_clock> end_time =
      std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = end_time - start_time;
  std::cout << " Elapsed time is : " << elapsed.count() << " " << std::endl;

  float maxError = 0.0f;
  for (int i = 0; i < N; i++) {
    maxError = fmax(maxError, Y[i] - 3.0f);
  }

  std::cout << "MaxError:" << maxError << std::endl;

  delete[] X;
  delete[] Y;
  return 0;
}
