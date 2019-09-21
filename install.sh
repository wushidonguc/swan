#!/bin/sh

echo '\n'
echo "Swan (Self-consistent Wannier-function-based quantum transport solver)"
echo "by Wushi Dong, 2019"
echo 'Compiling from source ...'
echo '\n'

# SPECIFY WHERE TO INSTALL EXECUTABLES
mkdir -p bin

cd src
mkdir -p obj
make -f makefile.self_consistent
make -f makefile.transport

echo '\n'
echo "Done!"
echo '\n'

mv self_consistent transport ../bin
cd ..
