#!/bin/bash

module swap PrgEnv-cray PrgEnv-gnu

cc DFTW_serial.c -o DFTW_serial.out -fopenmp -O2

export OMP_NUM_THREADS=32
srun -n 1 ./DFTW_serial.out > res_serial.txt