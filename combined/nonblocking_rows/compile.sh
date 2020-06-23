#!/bin/bash

module load mpi/openmpi-x86_64
mpic++ -openmp main.cpp -O2 -o main.mpi
