#!/bin/bash

#SBATCH --job-name=transport

#SBATCH --output=transport.out

#SBATCH --time=1:00:00

#SBATCH --exclusive

## receive an email when job starts, ends, and fails
##SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --ntasks=2

#SBATCH --partition=sandyb

# load the module(s)
module load openmpi
module load armadillo

# run
mpirun -n 2 ./transport

