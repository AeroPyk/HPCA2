#!/bin/bash

module swap PrgEnv-cray PrgEnv-gnu
cc stream.c -o stream.out -fopenmp

export OMP_NUM_THREADS=1
srun -n 5 ./stream.out > res1.txt

export OMP_NUM_THREADS=2
srun -n 5 ./stream.out > res2.txt

export OMP_NUM_THREADS=4
srun -n 5 ./stream.out > res4.txt

export OMP_NUM_THREADS=8
srun -n 5 ./stream.out > res8.txt

export OMP_NUM_THREADS=12
srun -n 5 ./stream.out > res12.txt

export OMP_NUM_THREADS=16
srun -n 5 ./stream.out > res16.txt

export OMP_NUM_THREADS=20
srun -n 5 ./stream.out > res20.txt

export OMP_NUM_THREADS=24
srun -n 5 ./stream.out > res24.txt

export OMP_NUM_THREADS=28
srun -n 5 ./stream.out > res28.txt

export OMP_NUM_THREADS=32
srun -n 5 ./stream.out > res32.txt