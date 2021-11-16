cmake_minimum_required(VERSION 3.10)
project(tsunami)

set(CMAKE_C_STANDARD 99)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O3")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(tsunami tsunami.c memory.c)

add_custom_command(
        TARGET tsunami POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
	${CMAKE_SOURCE_DIR}/compile_commands.json
	DEPENDS CMakeLists.txt)  
