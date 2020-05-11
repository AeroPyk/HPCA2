cmake_minimum_required(VERSION 3.15)
project(A2 C)

set(CMAKE_C_FLAGS -fopenmp)
set(CMAKE_C_STANDARD 99)

add_executable(Ex1 Exercise1.c)
add_executable(Ex2 Ex2/stream.c)
add_executable(Ex3 Ex3/Exercise3.c)
add_executable(Ex4 Ex4/DFTW_parallel.c)
add_executable(Ex5.1 Ex5/serial.c)
add_executable(Ex5.2 Ex5/parallel.c)
add_executable(Ex5.schedule Ex5/parallel_schedule.c)
