#!/bin/bash

# This script will simply go into the CPU and GPU folder and compile all executables into the bin folder.

INSTALLDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

make all -C $INSTALLDIR/cpu
make all -C $INSTALLDIR/cuda
make all -C $INSTALLDIR/opencl

