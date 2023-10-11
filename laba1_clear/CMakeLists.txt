cmake_minimum_required(VERSION 3.26)
project(laba1_clear)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_executable(laba1_clear main.cpp)
find_package(Threads REQUIRED)
add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
message(STATUS "flags - ${CMAKE_CXX_FLAGS}")