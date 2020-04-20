cmake_minimum_required(VERSION 3.15)
project(A2 C)

set(CMAKE_C_FLAGS -fopenmp)
set(CMAKE_C_STANDARD 99)

add_executable(A2 Exercise1.c)