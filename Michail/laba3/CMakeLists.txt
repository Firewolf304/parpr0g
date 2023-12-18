cmake_minimum_required(VERSION 3.27)
project(laba3)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lOpenCL")

add_executable(laba3 main.cpp)
