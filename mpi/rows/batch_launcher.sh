#!/bin/bash

#SBATCH --job-name=steady_heat
#SBATCH --reservation=fri
#SBATCH --output=steady_heat.log

#SBATCH --constraint=AMD
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks-per-socket=8
#SBATCH --time=00:10:00             # Time limit hrs:min:sec

module load mpi/openmpi-x86_64

mpirun -np 8 --map-by ppr:8:node --mca btl openib,self,vader mpi_rows.mpi 8 8 1 1 0.1

