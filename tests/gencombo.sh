#!/bin/bash

EXEC_NAME=$1
for p in 1 2 8 32
do
    for b in 8192 1024
    do
        
        FILENAME="launchers/combination_${p}_cores_$b.sh"
        echo "#!/bin/bash" > $FILENAME
        echo "#SBATCH --ntasks=$p" >> $FILENAME
        echo "#SBATCH --ntasks-per-node=$p" >> $FILENAME

        echo "#SBATCH -J 'heat'" >> $FILENAME
        echo "#SBATCH --reservation=fri" >> $FILENAME
        echo "#SBATCH --mem-per-cpu=1500M" >> $FILENAME
        echo "#SBATCH --time=00:30:00" >> $FILENAME
        echo "#SBATCH --output results/run_combination_${p}_cores_$b.log" >> $FILENAME

        echo "export OMP_PLACES=cores" >> $FILENAME
        echo "export OMP_PROC_BIND=TRUE" >> $FILENAME
        echo "export OMP_NUM_THREADS=$p" >> $FILENAME

        echo "for i in {1..10}" >> $FILENAME
        echo "do" >> $FILENAME
        echo "module load mpi/openmpi-x86_64" >> $FILENAME
        echo "mpirun -np \$SLURM_NTASKS -n \$SLURM_JOB_NUM_NODES --map-by ppr:$p:node --mca btl openib $EXEC_NAME $b $b 1 1" >> $FILENAME
        echo "done" >> $FILENAME

        echo "Processes: $p, board size: $b"
        sbatch $FILENAME --reservation=fri

        # rm $FILENAME
        sleep 0.5
    done
done
