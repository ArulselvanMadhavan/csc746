Readme file for chrono_timer demo program

# Purpose:
chrono_timer is intended to demonstrate the use of the c++ std::chrono library
for doing time measurement.

# Usage:

chrono_timer [secs]

With no arguments, the program will sleep() for the default value of 2 
seconds. If the user gives an optional argument "secs" on the command line,
the program will use that value, interpreted as an integer, and sleep
that number of seconds.

The program will first measure the current time, then sleep() the indicated
number of seconds, then measure the current time again, the compute
and print the difference in time.

# Default build instructions:

% cd chrono_timer     # contains the source file and CMakeLists.txt file

% mkdir build

% cd build

% cmake ../           # cmake generates lots of output

% make                # to build the program

# Additional build options

Optimization level: as written in the CMakeLists.txt file, the gcc/g++ 
compiler is provided an optimization flag of -O0, which will effectively
disable most optimizations.

To change this behavior, uncomment the line in CMakeLists.txt:
  
   #set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

and comment out the line:
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O0")

and then repeat the Default build instructions above, but first delete
everything (\rm -rf) in the build directory (cmake retains a memory of
previous runs and that memory can sometimes be unpredictable..you're best
off just deleting everything and starting fresh if you change anything
in the CMakeLists.txt file)



# EOF
