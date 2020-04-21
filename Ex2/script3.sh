#!/bin/bash

module swap PrgEnv-cray PrgEnv-gnu
cc stream_dynamic.c -o stream_dynamic.out -fopenmp

export OMP_NUM_THREADS=32
srun -n 5 ./stream_dynamic.out > res32_dynamic.txt

cc stream_guided.c -o stream_guided.out -fopenmp
srun -n 5 ./stream_guided.out > res32_guided.txt

cc stream_static.c -o stream_static.out -fopenmp
srun -n 5 ./stream_static.out > res32_static.txt