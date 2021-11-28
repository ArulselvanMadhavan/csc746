cmake_minimum_required(VERSION 3.10)
project(tsunami)

set(CMAKE_C_STANDARD 99)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# set(CMAKE_VERBOSE_MAKEFILE ON)

cmake_policy(SET CMP0060 NEW) ## use full path so it doesn't pick up hdf5 in /usr/lib
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O3 -pg")

# Find Packages
find_package(MPI REQUIRED)
set(HDF5_PREFER_PARALLEL true)
find_package(HDF5 1.10.1 REQUIRED)
# find_package(OpenMP)

if (NOT HDF5_IS_PARALLEL)
    message(FATAL_ERROR " -- HDF5 version is not parallel.")
endif (NOT HDF5_IS_PARALLEL)


add_executable(tsunami tsunami.c memory.c timer.c file_ops.c hdf5_file_ops.c)
set(IncludeDirs "${MPI_C_INCLUDE_PATH};${HDF5_C_INCLUDE_DIRS}")
set_target_properties(tsunami PROPERTIES INCLUDE_DIRECTORIES
"${IncludeDirs}")
set_target_properties(tsunami PROPERTIES COMPILE_FLAGS
"${MPI_C_COMPILE_FLAGS}")
set_target_properties(tsunami PROPERTIES LINK_FLAGS
  "${MPI_C_LINK_FLAGS}")

target_compile_options(tsunami PRIVATE -fopenmp -Wall -pedantic -march=native)
target_link_options(tsunami PRIVATE -fopenmp)
target_link_libraries(tsunami ${HDF5_C_LIBRARIES})
target_link_libraries(tsunami ${MPI_C_LIBRARIES})
target_link_libraries(tsunami -lm)

add_custom_command(
        TARGET tsunami POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
	${CMAKE_SOURCE_DIR}/compile_commands.json
	DEPENDS CMakeLists.txt)

configure_file(${CMAKE_SOURCE_DIR}/generate_image.py generate_image.py)
configure_file(${CMAKE_SOURCE_DIR}/generate_movie.sh generate_movie.sh)
