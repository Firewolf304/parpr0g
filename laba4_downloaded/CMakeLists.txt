cmake_minimum_required(VERSION 3.26)
project(laba1)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
add_executable(laba1 main.cpp)

target_link_libraries(laba1 omp)