#!/bin/bash

module swap PrgEnv-cray PrgEnv-gnu

cc DFTW_parallel.c -o DFTW_parallel.out -fopenmp -O2

export OMP_NUM_THREADS=32
srun -n 1 ./DFTW_parallel.out > res_parallel.txt