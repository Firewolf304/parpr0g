cmake_minimum_required(VERSION 3.27)
project(lab3 CUDA)

set(CMAKE_CUDA_STANDARD 11)

add_executable(lab3 src/main.cu
        include/includes.h)

set_target_properties(lab3 PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
