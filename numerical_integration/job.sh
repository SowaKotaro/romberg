#!/bin/sh

#$ -cwd
#$ -l node_f=4
#$ -l h_rt=0:20:00
#$ -N mpi_trapezoidal

module purge

module load gcc/14.2.0
module load openmpi/5.0.7-gcc

mpicc -o trapezoidal.out trapezoidal_mpi.c -lm

mpirun -npernode 192 -n 768 ./trapezoidal.out

