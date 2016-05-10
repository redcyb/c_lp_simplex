cmake_minimum_required(VERSION 3.3)
project(c_lp_simplex)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c99")

set(SOURCE_FILES main.c matrix.c matrix.h linear_program.c linear_program.h)
add_executable(c_lp_simplex ${SOURCE_FILES})
