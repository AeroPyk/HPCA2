#!/bin/bash

module swap PrgEnv-cray PrgEnv-gnu

cc Exercise3.c -o Exercise3.out -fopenmp

export OMP_NUM_THREADS=1
srun -n 1 ./Exercise3.out > res1.txt

export OMP_NUM_THREADS=2
srun -n 1 ./Exercise3.out > res2.txt

export OMP_NUM_THREADS=4
srun -n 1 ./Exercise3.out > res4.txt

export OMP_NUM_THREADS=8
srun -n 1 ./Exercise3.out > res8.txt

export OMP_NUM_THREADS=16
srun -n 1 ./Exercise3.out > res16.txt

export OMP_NUM_THREADS=20
srun -n 1 ./Exercise3.out > res20.txt

export OMP_NUM_THREADS=24
srun -n 1 ./Exercise3.out > res24.txt

export OMP_NUM_THREADS=28
srun -n 1 ./Exercise3.out > res28.txt

export OMP_NUM_THREADS=32
srun -n 1 ./Exercise3.out > res32.txt