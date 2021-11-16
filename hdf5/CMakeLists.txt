cmake_minimum_required(VERSION 3.0)
project(hdf5block2d LANGUAGES C)

# Export compile commands
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_C_STANDARD 99)
cmake_policy(SET CMP0060 NEW) ## use full path so it doesn't pick up hdf5 in /usr/lib
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")

# Find Packages
find_package(MPI REQUIRED)
set(HDF5_PREFER_PARALLEL true)
find_package(HDF5 1.10.1 REQUIRED)

if (NOT HDF5_IS_PARALLEL)
    message(FATAL_ERROR " -- HDF5 version is not parallel.")
endif (NOT HDF5_IS_PARALLEL)

add_executable(hdf5block2d hdf5block2d.c)

set(IncludeDirs "${MPI_C_INCLUDE_PATH};${HDF5_C_INCLUDE_DIRS}")
set_target_properties(hdf5block2d PROPERTIES INCLUDE_DIRECTORIES
"${IncludeDirs}")
set_target_properties(hdf5block2d PROPERTIES COMPILE_FLAGS
"${MPI_C_COMPILE_FLAGS}")
set_target_properties(hdf5block2d PROPERTIES LINK_FLAGS
"${MPI_C_LINK_FLAGS}")
target_link_libraries(hdf5block2d ${HDF5_C_LIBRARIES})
target_link_libraries(hdf5block2d ${MPI_C_LIBRARIES})

add_custom_command(
        TARGET hdf5block2d POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
	${CMAKE_SOURCE_DIR}/compile_commands.json
	DEPENDS CMakeLists.txt)  
