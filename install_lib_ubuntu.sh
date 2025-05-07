#!/bin/bash

#> monolis
git submodule update --init --recursive

cd submodule/monolis

#./install_lib.sh

#> monolis_utils
cd submodule/monolis_utils/
make
cd ../..

#> gedatsu
cd submodule/gedatsu/
make FLAGS=SUBMODULE
cd ../..

make FLAGS=MUMPS
