#!/bin/bash

EXEC_NAME=$1
for p in 1 2 8 32
do
    for n in 1 2
    do
        for b in 8192 1024
        do
            if (( p == 1 && n == 2 )); then
                continue
            fi
            FILENAME="launchers/launch_${p}_cores_${n}_nodes_$b.sh"
            echo "#!/bin/bash" > $FILENAME
            echo "#SBATCH --ntasks=$p" >> $FILENAME
            echo "#SBATCH --nodes=$n" >> $FILENAME

            NTASKS_PER_NODE=$p
            if (( n == 2 )); then
                NTASKS_PER_NODE=$(( $p / 2 ))
            fi

            TILE_DIM="1024 2048"
            if (( b == 8192 && p == 1 )); then
                TILE_DIM="8192 8192"
            fi
            if (( b == 8192 && p == 2 )); then
                TILE_DIM="8192 4096"
            fi
            if (( b == 8192 && p == 8 )); then
                TILE_DIM="4096 2048"
            fi
            if (( b == 8192 && p == 32 )); then
                TILE_DIM="2048 1024"
            fi
            if (( b == 1024 && p == 1 )); then
                TILE_DIM="1024 1024"
            fi
            if (( b == 1024 && p == 2 )); then
                TILE_DIM="1024 512"
            fi
            if (( b == 1024 && p == 8 )); then
                TILE_DIM="512 256"
            fi
            if (( b == 1024 && p == 32 )); then
                TILE_DIM="256 128"
            fi

            echo "#SBATCH --ntasks-per-node=$NTASKS_PER_NODE" >> $FILENAME
            echo "#SBATCH --ntasks-per-socket=8" >> $FILENAME
            echo "#SBATCH -J 'heat'" >> $FILENAME
            echo "#SBATCH --reservation=fri" >> $FILENAME
            echo "#SBATCH --mem-per-cpu=1500M" >> $FILENAME
            echo "#SBATCH --time=00:30:00" >> $FILENAME
            echo "#SBATCH --output results/run_${p}_cores_${n}_nodes_$b.log" >> $FILENAME

            echo "module load mpi/openmpi-x86_64" >> $FILENAME

            echo "for i in {1..10}" >> $FILENAME
            echo "do" >> $FILENAME
            # coll ^tuned
            echo "mpirun -np \$SLURM_NTASKS --map-by ppr:$NTASKS_PER_NODE:node --mca btl openib,self,vader $EXEC_NAME $b $b $TILE_DIM" >> $FILENAME
            echo "done" >> $FILENAME

            echo "Processes: $p, nodes: $n, board size: $b"
            sbatch $FILENAME --reservation=fri

            # rm $FILENAME
            sleep 0.5
        done
    done
done
