cmake_minimum_required(VERSION 3.14)
project(mpi_2dmesh_harness LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_BUILD_TYPE "Release")
#set(CMAKE_BUILD_TYPE "Debug")

find_package(MPI REQUIRED)

add_executable(mpi_2dmesh mpi_2dmesh.cpp)

# include directories
target_include_directories(mpi_2dmesh PRIVATE ${MPI_C_INCLUDE_PATH} ${MPI_CXX_INCLUDE_PATH})

# link libraries
target_link_libraries(mpi_2dmesh PRIVATE ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})

add_custom_command(
        TARGET mpi_2dmesh POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
	${CMAKE_SOURCE_DIR}/compile_commands.json
	DEPENDS CMakeLists.txt)
      
add_custom_command(
	TARGET mpi_2dmesh PRE_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy
	${CMAKE_SOURCE_DIR}/hostfile
	${CMAKE_CURRENT_BINARY_DIR}/hostfile
	DEPENDS CMakeLists.txt
)
