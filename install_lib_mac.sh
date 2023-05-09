#!/bin/bash

#> monolis
git submodule update --init --recursive
cd submodule/monolis
./install_lib.sh METIS MUMPS
make FLAGS=MPI,METIS,MUMPS

