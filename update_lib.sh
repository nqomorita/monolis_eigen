#!/bin/bash

#> monolis
cd submodule/monolis
git checkout .
git checkout master
git pull
make clean
make FLAGS=MPI,METIS
