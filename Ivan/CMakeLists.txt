cmake_minimum_required(VERSION 3.27)
project(Ivan)
set(CMAKE_CXX_COMPILER "/usr/bin/mpic++")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

add_executable(Ivan main.cpp)
