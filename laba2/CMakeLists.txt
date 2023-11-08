cmake_minimum_required(VERSION 3.26)
project(laba2 CUDA)

set(CMAKE_CUDA_STANDARD 23)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -lgomp")
add_executable(laba2 main.cu
        #mat.h
        #matrixer.h
        #matrixer.cu
)
#find_package(Threads REQUIRED)
add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)

set_target_properties(laba2 PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
