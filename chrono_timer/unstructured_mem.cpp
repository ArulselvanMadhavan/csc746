//
// (C) 2021, E. Wes Bethel
// chrono_timer: simple example showing how to use the std::chrono library
//  	to do elapsed time measurement
// usage:
// 	chrono_timer [secs]
//

#include <chrono>
#include <iostream>
#include <unistd.h>
#include <vector>

int main(int ac, char *av[]) {
  int problem_size = 100000;
  unsigned long long int sum = 0;
  int ptr = 0;
  if (ac > 1) {
    problem_size = std::atoi(av[1]);
    std::cout << "String value = " << av[1] << std::endl;
    std::cout << "Problem size value = " << problem_size << std::endl;
  }

  std::vector<unsigned long long> vect =
      std::vector<unsigned long long>(problem_size);

  // Initialize to rand integers
  for (int i = 0; i < problem_size - 1; i++) {
    vect[i] = lrand48() % problem_size;
  }

  std::chrono::time_point<std::chrono::high_resolution_clock> start_time =
      std::chrono::high_resolution_clock::now();

  for (int i = 0; i < problem_size; i++) {
    sum += vect[ptr];
    ptr = vect[ptr];
  }

  std::chrono::time_point<std::chrono::high_resolution_clock> end_time =
      std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = end_time - start_time;
  std::cout << " Sum is: " << sum << " " << std::endl;
  std::cout << " Elapsed time is : " << elapsed.count() << " " << std::endl;
  return 0;
}
