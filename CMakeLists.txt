cmake_minimum_required(VERSION 3.15)
project(A2 C)

set(CMAKE_C_FLAGS -fopenmp)
set(CMAKE_C_STANDARD 99)

add_executable(Ex1 Exercise1.c)
add_executable(Ex2 stream.c)
add_executable(Ex3 Exercise3.c)