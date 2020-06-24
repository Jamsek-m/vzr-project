#include <stdlib.h>
#include <stdio.h>
#include <chrono>
#include <math.h>
#include "/usr/include/openmpi-x86_64/mpi.h"
#include <sstream>

#define MAX_ITERS 300


double ** initTable(int m, int n) {
    double * fields = (double *) malloc(sizeof(double) * m * n);
    double ** table = (double **) malloc(sizeof(double*) * m);
    for (int k = 0; k < m; k++) {
        table[k] = &fields[k * n];
    }
    // redundant
    for (int k = 0; k < m; k++) {
        for (int l = 0; l < n; l++) {
            table[k][l] = -1.0; // rand() < 0.25 * RAND_MAX;
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

void print2DArrayBufferred(double ** table, int m, int n, int procId) {
    std::stringstream ss;
    ss.precision(6);
    ss << procId << ": [\n";
    for (int i = 0; i < m; i++) {
        ss << procId << ": [";
        for (int j = 0; j < n; j++) {
            if (j != 0) {
                ss << ", ";
            }
            ss << table[i][j];
        }
        ss << "]\n";
    }
    ss << "]\n";
    std::cout << ss.str() << std::endl;
}

void printArrayBufferred(double * table, int len, int procId) {
    std::stringstream ss;
    ss.precision(6);

    ss << procId << ": [\n";
    for (int i = 0; i < len; i++) {
        if (i != 0) {
            ss << ",";
        }
        ss << table[i] << "\n";
    }
    ss << "]";
    std::cout << ss.str() << std::endl;
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
    const int M = H / h;
    const int N = W / w;

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

    columnsNum = N / processNum;

    MPI_Datatype columnType, column;
    MPI_Type_vector(M, columnsNum, N, MPI_DOUBLE, &column);
    MPI_Type_commit(&column);
    MPI_Type_create_resized(column, 0, columnsNum * sizeof(double), &columnType);
    MPI_Type_commit(&columnType);

    // init table
    if (processId == 0) {
        table = initTable(M, N);
        pointerToTable = *table;
        initOuterTiles(table, M, N);
    }

    localTable = initTable(M, columnsNum);
    localTableNew = initTable(M, columnsNum);

    localLeftColumn = (double *) malloc(M * sizeof(double));
    localRightColumn = (double *) malloc(M * sizeof(double));

    auto start_time = std::chrono::high_resolution_clock::now();

    MPI_Scatter(pointerToTable, 1, columnType,
                *localTable, M * columnsNum, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    for (int i = 0; i < columnsNum; i++) {
        for (int j = 0; j < M; j++) {
            localTableNew[j][i] = localTable[j][i];
        }
    }
    
    MPI_Request requestSendLeft, requestSendRight;
    double w_h = 1 + pow(w/h, 2.0);
    double h_w = 1 + pow(h/w, 2.0); 

    MPI_Datatype sendColumnType, sendColumn;
    MPI_Type_vector(M, 1, columnsNum, MPI_DOUBLE, &sendColumn);
    MPI_Type_commit(&sendColumn);

    while (iters < MAX_ITERS) {
        /*
            Exchange data
        */
        if (processId > 0) {
            // Send left
            MPI_Isend(&localTable[0][0], 1, sendColumn, (processId - 1), 0, MPI_COMM_WORLD, &requestSendLeft);

            // Receive left
            MPI_Irecv(localLeftColumn, M, MPI_DOUBLE, (processId - 1), 1, MPI_COMM_WORLD, &requestSendRight);
            MPI_Wait(&requestSendRight, MPI_STATUS_IGNORE);
        }
        if (processId < processNum - 1) {
            // Send right
            MPI_Isend(&localTable[0][columnsNum - 1], 1, sendColumn, (processId + 1), 1, MPI_COMM_WORLD, &requestSendRight);

            // Receive right
            MPI_Irecv(localRightColumn, M, MPI_DOUBLE, (processId + 1), 0, MPI_COMM_WORLD, &requestSendLeft);
            MPI_Wait(&requestSendLeft, MPI_STATUS_IGNORE);
        }
        /*
            Compute new temperature
        */

        for (i = 1; i < M - 1; i++) {
            for (j = 0; j < columnsNum; j++) {

                if (processId == 0) {
                    // Skip first column
                    if (j == columnsNum - 1) {
                        // if last column
                        localTableNew[i][j] = 0.5 * ((localTable[i + 1][j] + localTable[i - 1][j]) / w_h + (localRightColumn[i] + localTable[i][j - 1]) / h_w);
                    } else if (j != 0) {
                        // all other columns
                        localTableNew[i][j] = 0.5 * ((localTable[i + 1][j] + localTable[i - 1][j]) / w_h + (localTable[i][j + 1] + localTable[i][j - 1]) / h_w);
                    }
                } else if (processId == processNum - 1) {
                    if (j == 0) {
                        localTableNew[i][j] = 0.5 * ((localTable[i + 1][j] + localTable[i - 1][j]) / w_h + (localTable[i][j + 1] + localLeftColumn[i]) / h_w);
                    } else if (j != columnsNum - 1) {
                        localTableNew[i][j] = 0.5 * ((localTable[i + 1][j] + localTable[i - 1][j]) / w_h + (localTable[i][j + 1] + localTable[i][j - 1]) / h_w);
                    }
                } else {
                    if (j == 0) {
                        localTableNew[i][j] = 0.5 * ((localTable[i + 1][j] + localTable[i - 1][j]) / w_h + (localTable[i][j + 1] + localLeftColumn[i]) / h_w);
                    } else if (j == columnsNum - 1) {
                        localTableNew[i][j] = 0.5 * ((localTable[i + 1][j] + localTable[i - 1][j]) / w_h + (localRightColumn[i] + localTable[i][j - 1]) / h_w);
                    } else {
                        localTableNew[i][j] = 0.5 * ((localTable[i + 1][j] + localTable[i - 1][j]) / w_h + (localTable[i][j + 1] + localTable[i][j - 1]) / h_w);
                    }
                }
            }
        }

        iters++;
        updateTables(&localTable, &localTableNew);

    }

    MPI_Gather(*localTable, columnsNum * M, MPI_DOUBLE,
               pointerToTable, 1, columnType,
               0, MPI_COMM_WORLD);
    
    if (processId == 0) {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto runningTime = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        std::cout << "MPI Columns Nonblocking. Execution time: " << runningTime << " ms." << std::endl;

        // print2DArray(table, M, N);
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

