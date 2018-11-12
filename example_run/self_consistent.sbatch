#!/bin/bash

#SBATCH --job-name=self_consistent

#SBATCH --output=self_consistent.out

#SBATCH --time=10:00:00

#SBATCH --exclusive

# receive an email when job starts, ends, and fails
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --ntasks=32

#SBATCH --partition=sandyb

# load the module(s)
module load openmpi
module load armadillo

# run
mpirun -n 32 ../self_consistent

