#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>


double ** initTable(int m, int n) {
    double ** table = (double **)malloc(m * sizeof(double *));
    for (int i = 0; i < m; i++) {
        table[i] = (double *)malloc(n * sizeof(double));
    }
    return table;
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

void initArgs(int argc, const char * argv[], int * table_w, int * table_h, int * tile_w, int * tile_h, double * eps) {
    // Defaults
    *table_w = 20;
    *table_h = 20;
    *tile_h = 1;
    *tile_w = 1;
    *eps = 0.001;

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
        *eps = atof(argv[5]);
    }
}

int main(int argc, const char * argv[]) {
    double average = 0.0;
    double localDiff;
    int i, j;
    int W, H, w, h;
    double EPSILON;
    initArgs(argc, argv, &W, &H, &w, &h, &EPSILON);
    int M = H / h;
    int N = W / w;

    double ** table = initTable(M, N);
    double ** tableOld = initTable(M, N);

    /*
      Initializing table borders
    */

    #pragma omp parallel shared ( table ) private ( i, j )
    {
        #pragma omp for
        for (i = 1; i < M - 1; i++) {
            table[i][0] = 100.0;
            table[i][N - 1] = 100.0;
        }

        #pragma omp for
        for (j = 0; j < N; j++) {
            table[M - 1][j] = 100.0;
            table[0][j] = 0.0;
        }

        // Setting average value in middle tiles
        #pragma omp for reduction ( + : average )
        for (i = 1; i < M - 1; i++) {
            average = average + table[i][0] + table[i][N - 1];
        }
        #pragma omp for reduction ( + : average )
        for (j = 0; j < N; j++) {
            average = average + table[M - 1][j] + table[0][j];
        }
    }

    // Values for average are summed up, need to divide
    average = average / (double)(2 * M + 2 * N - 4);

    // Set average value to all middle tiles
    #pragma omp parallel shared(average, table) private(i, j)
    {
        #pragma omp for
        for (i = 1; i < M - 1; i++) {
            for (j = 1; j < N - 1; j++) {
                table[i][j] = average;
            }
        }
    }

    /*
      Algorithm start
    */
    int iterations = 0;
    double diff = EPSILON;
    double startTime = omp_get_wtime();
    while (EPSILON <= diff) {

        #pragma omp parallel shared ( table, tableOld ) private ( i, j )
        {
            // Copy table to old table
            #pragma omp for
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < N; j++) {
                    tableOld[i][j] = table[i][j];
                }
            }

            // Calculate new values
            #pragma omp for
            for (int i = 1; i < M - 1; i++) {
                for (int j = 1; j < N - 1; j++) {
                    table[i][j] = 0.5 * ((tableOld[i + 1][j] + tableOld[i - 1][j]) / (1 + pow(w / h, 2.0)) + ((tableOld[i][j + 1] + tableOld[i][j - 1]) / (1 + pow(h / w, 2.0))));
                }
            }
        }

        diff = 0.0;
        #pragma omp parallel shared (diff, table, tableOld) private (i, j, localDiff)
        {
            localDiff = 0.0;
            #pragma omp for
            for (int i = 1; i < M - 1; i++) {
                for (int j = 1; j < N - 1; j++) {
                    double calculatedDiff = fabs(table[i][j] - tableOld[i][j]);
                    if (localDiff < calculatedDiff) {
                        localDiff = calculatedDiff;
                    }
                }
            }

            #pragma omp critical
            {
                if (diff < localDiff) {
                    diff = localDiff;
                }
            }
        }

        iterations++;
        // std::cout << "Diff: " << diff << ", Eps: " << EPSILON << std::endl;
    }

    double endTime = omp_get_wtime();

    printf("Running time: %f s\n", endTime - startTime);

    std::cout << "Number of iterations: " << iterations << std::endl;

    // print2DArray(table, M, N);

    return 0;
}
