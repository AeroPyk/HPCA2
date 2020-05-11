#!/bin/bash

module swap PrgEnv-cray PrgEnv-gnu

export OMP_NUM_THREADS=1
cc parallel_red.c -o parallel_red.out -fopenmp -O2 -lm
srun -n 8 ./parallel_red.out > res-parallel_red-T1.txt

export OMP_NUM_THREADS=4
srun -n 8 ./parallel_red.out > res-parallel_red-T4.txt

export OMP_NUM_THREADS=8
srun -n 8 ./parallel_red.out > res-parallel_red-T8.txt

export OMP_NUM_THREADS=12
srun -n 8 ./parallel_red.out > res-parallel_red-T12.txt

export OMP_NUM_THREADS=16
srun -n 8 ./parallel_red.out > res-parallel_red-T16.txt

export OMP_NUM_THREADS=20
srun -n 8 ./parallel_red.out > res-parallel_red-T20.txt

export OMP_NUM_THREADS=24
srun -n 8 ./parallel_red.out > res-parallel_red-T24.txt

export OMP_NUM_THREADS=28
srun -n 8 ./parallel_red.out > res-parallel_red-T28.txt

export OMP_NUM_THREADS=32
srun -n 8 ./parallel_red.out > res-parallel_red-T32.txt