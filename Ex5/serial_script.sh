#!/bin/bash

module swap PrgEnv-cray PrgEnv-gnu

cc serial_500.c -o serial_500.out -fopenmp -O2 -lm
srun -n 4 ./serial_500.out > res-serial_500.txt

cc serial_1000.c -o serial_1000.out -fopenmp -O2 -lm
srun -n 4 ./serial_1000.out > res-serial_1000.txt

cc serial_2000.c -o serial_2000.out -fopenmp -O2 -lm
srun -n 4 ./serial_2000.out > res-serial_2000.txt