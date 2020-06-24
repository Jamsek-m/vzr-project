#!/bin/bash

EXEC_NAME=$1
for p in 1 2 8 32
do
    for b in 8192 1024
    do
        
        FILENAME="launchers/openmp_${p}_cores_$b.sh"
        echo "#!/bin/bash" > $FILENAME
        echo "#SBATCH --ntasks=1" >> $FILENAME
        echo "#SBATCH --cpus-per-task=$p" >> $FILENAME

        echo "#SBATCH -J 'heat'" >> $FILENAME
        echo "#SBATCH --reservation=fri" >> $FILENAME
        echo "#SBATCH --mem-per-cpu=1500M" >> $FILENAME
        echo "#SBATCH --time=00:30:00" >> $FILENAME
        echo "#SBATCH --output results/run_openmp_${p}_cores_$b.log" >> $FILENAME
        echo "#SBATCH --constraint=AMD" >> $FILENAME

        echo "export OMP_PLACES=cores" >> $FILENAME
        echo "export OMP_PROC_BIND=TRUE" >> $FILENAME
        echo "export OMP_NUM_THREADS=$p" >> $FILENAME

        echo "for i in {1..10}" >> $FILENAME
        echo "do" >> $FILENAME
        echo "srun $EXEC_NAME $b $b" >> $FILENAME
        echo "done" >> $FILENAME

        echo "Processes: $p, nodes: $n, board size: $b"
        sbatch $FILENAME --reservation=fri

        # rm $FILENAME
        sleep 0.5
    done
done
