#!/bin/sh

# SPECIFY WHERE TO INSTALL EXECUTABLES
INSTALL_DIR=../bin

cd src
make -f makefile.self_consistent
make -f makefile.transport
mv self_consistent transport ${INSTALL_DIR}
cd ..
