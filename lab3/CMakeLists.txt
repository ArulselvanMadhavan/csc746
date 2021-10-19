cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
project(vec_add LANGUAGES CXX)

set(CMAKE_BUILD_TYPE "Release")

add_executable(vecadd_omp_offload vecadd-omp-offload.cpp)
target_compile_options(vecadd_omp_offload PRIVATE -Ofast -fopenmp -fopenmp-targets=nvptx64)
target_link_options(vecadd_omp_offload PRIVATE -Ofast -fopenmp -fopenmp-targets=nvptx64)
