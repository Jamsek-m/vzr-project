#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-socket=8
#SBATCH -J 'heat'
#SBATCH --reservation=fri
#SBATCH --mem-per-cpu=1500M
#SBATCH --time=00:30:00
#SBATCH --output run_1_cores_1_nodes_8192_eth_rows_blocking.log
module load mpi/openmpi-x86_64
for i in {1..10}
do
mpirun -np $SLURM_NTASKS --map-by ppr:1:node --mca btl tcp,self,vader mpi_rows_blocking.mpi 8192 8192
done
