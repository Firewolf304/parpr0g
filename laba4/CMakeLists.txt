cmake_minimum_required(VERSION 3.27)
project(laba4)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lOpenCL")
set(CMAKE_CXX_STANDARD 11)

add_executable(laba4 main.cpp)
#find_package(OpenCL REQUIRED)
#target_link_libraries(laba4 OpenCL)