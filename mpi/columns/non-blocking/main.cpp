#include <stdlib.h>
#include <stdio.h>
#include <chrono>
#include <math.h>
#include "/usr/include/openmpi-x86_64/mpi.h"

#define MAX_ITERS 1000


double ** initTable(int m, int n) {
    double * fields = (double *) malloc(sizeof(double) * m * n);
    double ** table = (double **) malloc(sizeof(double*) * m);
    for (int k = 0; k < m; k++) {
        table[k] = &fields[k * n];
    }
    // redundant
    for (int k = 0; k < m; k++) {
        for (int l = 0; l < n; l++) {
            table[k][l] = rand() < 0.25 * RAND_MAX;
        }
    }
    return table;
}

void initOuterTiles(double ** table, int m, int n) {
    for (int i = 1; i < m; i++) {
        table[i][0] = 100.0;
        table[i][n - 1] = 100.0;
    }
    for (int j = 0; j < n; j++) {
        table[m - 1][j] = 100.0;
        table[0][j] = 0.0;
    }

    double average = 0.0;
    for (int i = 0; i < m; i++) {
        average = average + table[i][0] + table[i][n- 1];
    }
    for (int j = 0; j < n; j++) {
        average = average + table[m - 1][j] + table[0][j];
    }
    average = average / (double) (2 * m + 2 * n - 4);
    for (int i = 1; i < m - 1; i++) {
        for (int j = 1; j < n - 1; j++) {
            table[i][j] = average;
        }
    }
}

void print2DArray(double ** table, int m, int n) {
    printf("[\n");
    for (int i = 0; i < m; i++) {
        printf("[");
        for (int j = 0; j < n; j++) {
            if (j != 0) {
                printf(", ");
            }
            printf("%6.2f", table[i][j]);
        }
        printf("]\n");
    }
    printf("]\n");
}

void print2DArrayParallel(double ** table, int m, int n, int procId) {
    printf("%d: [\n", procId);
    for (int i = 0; i < m; i++) {
        printf("%d: [", procId);
        for (int j = 0; j < n; j++) {
            if (j != 0) {
                printf(", ");
            }
            printf("%6.2f", table[i][j]);
            fflush(stdout);
        }
        printf("]: %d\n", procId);
        fflush(stdout);
    }
    printf("]: %d\n", procId);
    fflush(stdout);
}

void initArgs(int argc, char ** argv, int * table_w, int * table_h, int * tile_w, int * tile_h) {
    // Defaults
    *table_w = 20;
    *table_h = 20;
    *tile_h = 1;
    *tile_w = 1;

    if (argc < 3) {
        return;
    } else if (argc >= 3 && argc < 5) {
        *table_w = atoi(argv[1]);
        *table_h = atoi(argv[2]);
    } else if (argc == 5) {
        *table_w = atoi(argv[1]);
        *table_h = atoi(argv[2]);
        *tile_w = atoi(argv[3]);
        *tile_h = atoi(argv[4]);
    } else if (argc >= 6) {
        *table_w = atoi(argv[1]);
        *table_h = atoi(argv[2]);
        *tile_w = atoi(argv[3]);
        *tile_h = atoi(argv[4]);
    }
}

void updateTables(double *** b, double *** bn) {
    double** bt;
	bt = *b;
	*b = *bn;
	*bn = bt;
}

void freeTable(double ** table) {
    free(*table);
    free(table);
}

int main(int argc, char * argv[]) {
    int W, H, w, h;
    initArgs(argc, argv, &W, &H, &w, &h);
    int M = H / h;
    int N = W / w;

    double * pointerToTable = NULL;
    double ** table;

    int processId, processNum;
    int columnsNum;
    double ** localTable;
    double ** localTableNew;
    double * localLeftColumn, * localRightColumn;

    int i, j;
    int iters = 0;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &processId);
    MPI_Comm_size(MPI_COMM_WORLD, &processNum);

    MPI_Datatype columnType, column;
    MPI_Type_vector(M, 1, N, MPI_DOUBLE, &column);
    MPI_Type_commit(&column);
    MPI_Type_create_resized(column, 0, 1 * sizeof(double), &columnType);
    MPI_Type_commit(&columnType);

    // init table
    if (processId == 0) {
        table = initTable(M, N);
        pointerToTable = *table;
        initOuterTiles(table, M, N);
    }

    columnsNum = N / processNum;
    localTable = initTable(M, columnsNum);
    localTableNew = initTable(M, columnsNum);

    // print2DArrayParallel(localTableNew, M, columnsNum, processId);

    localLeftColumn = (double *) malloc(M * sizeof(double));
    localRightColumn = (double *) malloc(M * sizeof(double));

    // Set initial data in local table new
    if (processId == 0) {
        // first column - first element is 0.0, others are 100.0
        for (j = 0; j < M; j++) {
            localTableNew[j][0] = 100.0;
        }
    } else if (processId == (processNum - 1)) {
        // last column - first element is 0.0, others are 100.0
        for (j = 0; j < M; j++) {
            localTableNew[j][columnsNum - 1] = 100.0;
        }
    }
    for (i = 0; i < columnsNum; i++) {
        // set bottom row to 100.0
        localTableNew[M - 1][i] = 100.0;
        // set top row to 0.0
        localTableNew[0][i] = 0.0;
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    MPI_Scatter(pointerToTable, columnsNum, columnType,
                *localTable, M * columnsNum, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    // print2DArray(localTable, M, columnsNum);
    
    MPI_Request requestLeft, requestRight;
    double w_h = 1 + pow(w/h, 2.0);
    double h_w = 1 + pow(h/w, 2.0); 

    while (iters < MAX_ITERS) {
        /*
            Exchange data
        */
        if (processId == 0 && processId != processNum - 1) {
            // First column and if only one column, not last
            MPI_Isend(localTable[columnsNum - 1], N, MPI_DOUBLE, (processId + processNum + 1) % processNum, 0, MPI_COMM_WORLD, &requestLeft);
            MPI_Irecv(localRightColumn, N, MPI_DOUBLE, (processId + processNum + 1) % processNum, 1, MPI_COMM_WORLD, &requestRight);
            MPI_Wait(&requestLeft, MPI_STATUS_IGNORE);
        } else if (processId > 0 && processId == processNum - 1) {
            // Last column and if only one column, not first
            MPI_Isend(localTable[0], N, MPI_DOUBLE, (processId + processNum - 1) % processNum, 1, MPI_COMM_WORLD, &requestLeft);
            MPI_Irecv(localLeftColumn, N, MPI_DOUBLE, (processId + processNum - 1) % processNum, 0, MPI_COMM_WORLD, &requestRight);
            MPI_Wait(&requestRight, MPI_STATUS_IGNORE);
        } else if (processId > 0 && processId < processNum - 1) {
            // Middle columns
            MPI_Isend(localTable[0], N, MPI_DOUBLE, (processId + processNum - 1) % processNum, 1, MPI_COMM_WORLD, &requestLeft);
            MPI_Irecv(localRightColumn, N, MPI_DOUBLE, (processId + processNum + 1) % processNum, 1, MPI_COMM_WORLD, &requestRight);

            MPI_Isend(localTable[columnsNum - 1], N, MPI_DOUBLE, (processId + processNum + 1) % processNum, 0, MPI_COMM_WORLD, &requestLeft);
            MPI_Irecv(localLeftColumn, N, MPI_DOUBLE, (processId + processNum - 1) % processNum, 0, MPI_COMM_WORLD, &requestLeft);

            MPI_Wait(&requestLeft, MPI_STATUS_IGNORE);
            MPI_Wait(&requestRight, MPI_STATUS_IGNORE);
        }

        /*
            Compute new temperature
        */
        // left:
        // localLeftColumn[j]
        // right:
        // localRightColumn[j]
        // bottom:
        // localTable[i][j + 1]
        // top:
        // localTable[i][j - 1]
        for (i = 0; i < columnsNum; i++) {
            for (j = 1; j < M - 1; j++) {
                if (columnsNum == 1 && processId > 0 && processId < processNum - 1) {
                    // Only one column
                    localTableNew[i][j] = 0.5 * ((localTable[i][j + 1] + localTable[i][j - 1]) / w_h + (localRightColumn[j] + localLeftColumn[j]) / h_w);
                } else if (i == 0 && processId != 0 && columnsNum > 1) {
                    // Multiple rows
                    localTableNew[i][j] = 0.5 * ((localTable[i][j + 1] + localTable[i][j - 1]) / w_h + (localTable[i + 1][j] + localLeftColumn[j]) / h_w);
                } else if (columnsNum != 1 && i == columnsNum - 1 && processId != processNum - 1) {
                    // Multiple rows, last row in block, but not last row in table
                    localTableNew[i][j] = 0.5 * ((localTable[i][j + 1] + localTable[i][j - 1]) / w_h + (localRightColumn[j] + localTable[i - 1][j]) / h_w);
                } else if (columnsNum != 1 && i > 0 && i < columnsNum - 1) {
                    // All rows in between
                    localTableNew[i][j] = 0.5 * ((localTable[i][j + 1] + localTable[i][j - 1]) / w_h + (localTable[i + 1][j] + localTable[i - 1][j]) / h_w);
                }
            }
        }

        iters++;
        updateTables(&localTable, &localTableNew);

    }

    MPI_Gather(*localTable, columnsNum * M, MPI_DOUBLE,
               pointerToTable, columnsNum * M, MPI_DOUBLE,
               0, MPI_COMM_WORLD);
    
    if (processId == 0) {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto runningTime = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        std::cout << "Execution time: " << runningTime << " ms." << std::endl;

        print2DArray(table, M, N);
    }

    if (processId == 0) {
        freeTable(table);
    }
    freeTable(localTable);
    free(localLeftColumn);
    free(localRightColumn);

    MPI_Finalize();

    return 0;
}

