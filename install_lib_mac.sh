#!/bin/bash

#> monolis
git submodule update --init --recursive
cd submodule/monolis
./install_lib.sh MUMPS
make FLAGS=MUMPS

