#!/bin/bash

#SBATCH --ntasks=8 # cores
#SBATCH --nodes=1
#SBATCH --ntasks-per-socket=4
#SBATCH --ntasks-per-node=8 # cores

# --cpus-per-task=2
# --constraint=AMD

#SBATCH --mem-per-cpu=100M
#SBATCH -J 'heat'
#SBATCH --reservation=fri
#SBATCH --output=out/log.txt

module load mpi/openmpi-x86_64
mpirun --report-bindings -n $SLURM_JOB_NUM_NODES -np $SLURM_NTASKS --map-by ppr:8:node --mca btl openib ./main.mpi 6000 6000
