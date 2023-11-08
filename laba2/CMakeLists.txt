cmake_minimum_required(VERSION 3.26)
project(laba2 CUDA)

set(CMAKE_CUDA_STANDARD 23)

add_executable(laba2 main.cu
        #mat.h
        #matrixer.h
        #matrixer.cu
)

set_target_properties(laba2 PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
