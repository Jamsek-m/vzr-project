#!/bin/bash

module load mpi/openmpi-x86_64
mpic++ $1.cpp -O2 -o $1.mpi
