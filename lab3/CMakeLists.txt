cmake_minimum_required(VERSION 3.27)
project(lab3)

set(CMAKE_CXX_COMPILER "/usr/bin/mpic++")
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

add_executable(lab3 main.cpp)
