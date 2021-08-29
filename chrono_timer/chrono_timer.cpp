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

int main(int ac, char *av[]) {
  int sleep_secs = 2; // default sleep for 2 secs
  int problem_size = 100000;

  // if (ac > 1) // if user entered a value after the prog name, parse it
  //   sleep_secs = std::atoi(av[1]);

  if (ac > 1)
    problem_size = std::atoi(av[1]);
  std::cout << "String value = " << av[1] << std::endl;
  std::cout << "Problem size value = " << problem_size << std::endl;
  
  std::cout << "Going to sleep for " << sleep_secs << " seconds " << std::endl;

  std::chrono::time_point<std::chrono::high_resolution_clock> start_time =
      std::chrono::high_resolution_clock::now();

  sleep(sleep_secs);

  std::chrono::time_point<std::chrono::high_resolution_clock> end_time =
      std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = end_time - start_time;

  std::cout << " Elapsed time is : " << elapsed.count() << " " << std::endl;
  return 0;
}
