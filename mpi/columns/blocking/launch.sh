#!/bin/bash

#SBATCH --ntasks=32 # cores
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32 # cores

#SBATCH --mem-per-cpu=100M
# --constraint=AMD
#SBATCH -J 'heat'
#SBATCH --reservation=fri
#SBATCH --output=out/log.txt

module load mpi/openmpi-x86_64
mpirun -np $SLURM_NTASKS --map-by ppr:32:node --mca coll ^tuned ./main.mpi 8096 8096
