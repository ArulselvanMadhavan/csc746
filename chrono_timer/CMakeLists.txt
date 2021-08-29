
cmake_minimum_required( VERSION 3.5)

project (chrono_timer LANGUAGES CXX)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-Wall")  # uncomment this line to turn on compiler warnings
set( CMAKE_VERBOSE_MAKEFILE on )
# for CSC 746 Fall 2021 Homework#1, the simplest way to specify either -O3 or -O0 compiler
# optimization is via the command line, thus:
#
# for -O3 optimization:
# % cmake ../   # no arguments to cmake will result in -O3 optimization
# for -\O0 optimization:
# % cmake ../ -DCMAKE_CXX_FLAGS_RELEASE="-O0"
#
# the stuff below here remains for reference
# uncomment the set(CMAKE_CXX_FLAGS ... ) line for either -O3 or -O0
#
# -O3 is full optimization in gcc/g++
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
# -O0 is no optimization in gcc/g++
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O0")
 
add_executable (no_mem no_mem.cpp)
add_executable (structured_mem structured_mem.cpp)

