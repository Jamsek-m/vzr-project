#include <stdlib.h>
#include <stdio.h>
#include "/usr/include/openmpi-x86_64/mpi.h"

#define HEIGHT 3
#define WIDTH 4

void print2DArray(int ** table, int m, int n) {
    printf("[\n");
    for (int i = 0; i < m; i++) {
        printf("[");
        for (int j = 0; j < n; j++) {
            if (j != 0) {
                printf(", ");
            }
            printf("%d", table[i][j]);
        }
        printf("]\n");
    }
    printf("]\n");
}

int ** initTable(int m, int n) {
    int * fields = (int *) malloc(sizeof(int) * m * n);
    int ** table = (int **) malloc(sizeof(int) * m);
    for (int k = 0; k < m; k++) {
        table[k] = &fields[k * n];
    }
    return table;
}

void fillTableWithData(int ** table, int m, int n) {
    int count = 0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            table[i][j] = count;
            count++;
        }
    }
}

void logProcess(int procId, const char * msg) {
    std::cout << "PROC " << procId << ": " << msg << std::endl;
}

int main(int argc, char * argv[]) {

    int * pointerToTable = NULL;
    int ** table;
    int ** localTable;
    int processNum, processId;
    MPI_Status status;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &processNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &processId);

    if (processId == 0) {
        table = initTable(HEIGHT, WIDTH);
        pointerToTable = *table;
        fillTableWithData(table, HEIGHT, WIDTH);
        print2DArray(table, HEIGHT, WIDTH);
    }

    int columnsNum = WIDTH / processNum;

    // Each process initialize its own table
    localTable = initTable(HEIGHT, columnsNum);

    // Scatter data to processes
    MPI_Datatype columnType, column;
    MPI_Type_vector(HEIGHT, 1, WIDTH, MPI_INT, &column);
    MPI_Type_commit(&column);
    MPI_Type_create_resized(column, 0, 1 * sizeof(int), &columnType);
    MPI_Type_commit(&columnType);
    MPI_Scatter(pointerToTable, columnsNum, columnType,
                *localTable, HEIGHT * columnsNum, MPI_INT,
                0, MPI_COMM_WORLD);


    // Each process print its own table
    logProcess(processId, "PRINT ARRAY");
    print2DArray(localTable, HEIGHT, columnsNum);

    if (processId == 0) {
        free(*table);
        free(table);
    }
    free(*localTable);
    free(localTable);

    MPI_Finalize();

    return 0;
}