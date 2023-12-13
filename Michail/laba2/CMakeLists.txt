cmake_minimum_required(VERSION 3.27)
project(laba2 CUDA)

set(CMAKE_CUDA_STANDARD 11)

add_executable(laba2 main.cu)

set_target_properties(laba2 PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
