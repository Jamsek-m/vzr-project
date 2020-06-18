#!/bin/bash

mpic++ mpi_rows.cpp -O2 -o mpi_rows.mpi

#sbatch --reservation=fri batch_launcher.sh