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

export OMP_PLACES=cores
export OMP_PROC_BIND=TRUE
export OMP_NUM_THREADS=8

module load mpi/openmpi-x86_64
mpirun --report-bindings -n $SLURM_JOB_NUM_NODES -np $SLURM_NTASKS --map-by ppr:8:node --mca btl openib ./main.mpi 512 512 1 1
