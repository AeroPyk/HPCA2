#!/bin/bash

module swap PrgEnv-cray PrgEnv-gnu

export OMP_NUM_THREADS=1
cc parallel_500.c -o parallel_500.out -fopenmp -O2 -lm
srun -n 8 ./parallel_500.out > res-parallel_500-T1.txt

cc parallel_1000.c -o parallel_1000.out -fopenmp -O2 -lm
srun -n 8 ./parallel_1000.out > res-parallel_1000-T1.txt

export OMP_NUM_THREADS=4
srun -n 8 ./parallel_500.out > res-parallel_500-T4.txt
srun -n 8 ./parallel_1000.out > res-parallel_1000-T4.txt

export OMP_NUM_THREADS=8
srun -n 8 ./parallel_500.out > res-parallel_500-T8.txt
srun -n 8 ./parallel_1000.out > res-parallel_1000-T8.txt

export OMP_NUM_THREADS=12
srun -n 8 ./parallel_500.out > res-parallel_500-T12.txt
srun -n 8 ./parallel_1000.out > res-parallel_1000-T12.txt

export OMP_NUM_THREADS=16
srun -n 8 ./parallel_500.out > res-parallel_500-T16.txt
srun -n 8 ./parallel_1000.out > res-parallel_1000-T16.txt

export OMP_NUM_THREADS=24
srun -n 8 ./parallel_500.out > res-parallel_500-T24.txt
srun -n 8 ./parallel_1000.out > res-parallel_1000-T24.txt

export OMP_NUM_THREADS=32
srun -n 5 ./parallel_500.out > res-parallel_500-T32.txt
srun -n 5 ./parallel_1000.out > res-parallel_1000-T32.txt