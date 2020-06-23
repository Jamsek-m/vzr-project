#!/bin/bash

#SBATCH --job-name=steady_heat
#SBATCH --reservation=fri
#SBATCH --output=steady_heat.log

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-socket=8
#SBATCH --time=00:10:00             # Time limit hrs:min:sec

module load mpi/openmpi-x86_64

mpirun -np 16 --map-by ppr:16:node --mca btl openib,self,vader blocks.mpi 32 32 8 8