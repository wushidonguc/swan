#!/bin/sh

# SPECIFY WHERE EXECUTABLES ARE INSTALLED
INSTALL_DIR=../bin


# CREATE OUT/ DIRECTORY
if [ ! -d out ];then
  mkdir -p out
fi

# RUN

# Self-consistent simulation
mpirun -n 32 ${INSTALL_DIR}/self_consistent

# Calculate transport properties based on converged electrostatics
mpirun -n 32 ${INSTALL_DIR}/transport
