cmake_minimum_required(VERSION 3.26)
project(lab1)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_executable(lab1 main.cpp)



find_package(Threads REQUIRED)
add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
message(STATUS "flags - '${CMAKE_CXX_FLAGS}'")