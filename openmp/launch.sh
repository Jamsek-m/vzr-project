#!/bin/sh
#SBATCH --job-name=heat
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --constraint=AMD
#SBATCH --output=out/run.log
#SBATCH --reservation=fri
#SBATCH --time=00:10:00             # Time limit hrs:min:sec

export OMP_PLACES=cores
export OMP_PROC_BIND=TRUE
export OMP_NUM_THREADS=64

srun prog 500 500 1 1 0.001

wait
