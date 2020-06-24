#!/bin/bash

module load mpi/openmpi-x86_64
mpic++ -openmp combined.cpp -O2 -o combined.mpi
