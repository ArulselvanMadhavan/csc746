#include <math.h>
#include <iostream>

void add(int N, float* X, float* Y){
  for (int i = 0; i < N; i++) {
    Y[i] = X[i]+Y[i];
  }
}

int main(void){
  const int N = 1 << 30;

  float* X = new float[N];
  float* Y = new float[N];
  
  // Initialize
  for(int i = 0; i < N; i++){
    X[i] = 1.0f;
    Y[i] = 2.0f;
  }

  add(N, X, Y);

  float maxError = 0.0f;
  for(int i = 0; i < N; i++){
    maxError = fmax(maxError, Y[i] - 3.0f);
  }

  std::cout << "MaxError:" << maxError << std::endl;

  delete []X;
  delete []Y;
  return 0;
}
