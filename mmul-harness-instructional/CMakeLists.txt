cmake_minimum_required(VERSION 3.14)
project(hw2 LANGUAGES C CXX)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_VERBOSE_MAKEFILE on)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
find_package(BLAS REQUIRED)

# We require the GNU compiler for this assignment.
if (NOT "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set(Prg-Intel "PrgEnv-intel")
    set(Prg-Clang "PrgEnv-cray")
    message(WARNING
            "Must use GNU Compiler for submission. Make sure you ran:\n"
            "module swap ${Prg-${CMAKE_C_COMPILER_ID}} PrgEnv-gnu")
endif ()

# Common library target for benchmarking.
add_library(benchmark OBJECT benchmark.cpp)
target_compile_features(benchmark PRIVATE cxx_std_11)

add_executable(benchmark-basic dgemm-basic.cpp)
target_link_libraries(benchmark-basic PRIVATE benchmark ${BLAS_LIBRARIES})
target_compile_features(benchmark-basic PRIVATE c_std_11 c_restrict)
target_compile_options(benchmark-basic PRIVATE -Wall -pedantic -march=native)
set(BENCHMARK "basic")
configure_file(job.in job-${BENCHMARK})

add_executable(benchmark-basic-copy dgemm-basic-copy.cpp)
target_link_libraries(benchmark-basic-copy PRIVATE benchmark ${BLAS_LIBRARIES})
target_compile_features(benchmark-basic-copy PRIVATE c_std_11 c_restrict)
target_compile_options(benchmark-basic-copy PRIVATE -Wall -pedantic -march=native)
set(BENCHMARK "basic-copy")
configure_file(job.in job-${BENCHMARK})

add_executable(benchmark-blas dgemm-blas.cpp)
target_link_libraries(benchmark-blas PRIVATE benchmark ${BLAS_LIBRARIES})
target_compile_features(benchmark-blas PRIVATE c_std_11 c_restrict)
target_compile_options(benchmark-blas PRIVATE -Wall -pedantic -march=native)
set(BENCHMARK "blas")
configure_file(job.in job-${BENCHMARK})

add_executable(benchmark-blocked dgemm-blocked.cpp)
add_library(bench-blocked OBJECT benchmark.cpp)
target_compile_definitions(bench-blocked PRIVATE -DBLOCKED)
target_compile_features(bench-blocked PRIVATE cxx_std_11)
target_link_libraries(benchmark-blocked PRIVATE bench-blocked ${BLAS_LIBRARIES})
target_compile_features(benchmark-blocked PRIVATE c_std_11 c_restrict)
target_compile_options(benchmark-blocked PRIVATE -Wall -pedantic -march=native)
set(BENCHMARK "blocked")
configure_file(job.in job-${BENCHMARK})

