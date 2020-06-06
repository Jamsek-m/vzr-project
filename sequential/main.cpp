#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


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
    for (int i = 1; i < M - 1; i++) {
        table[i][0] = 100.0;
    }
    for (int i = 1; i < M - 1; i++) {
        table[i][N - 1] = 100.0;
    }
    for (int j = 0; j < N; j++) {
        table[M - 1][j] = 100.0;
    }
    for (int j = 0; j < N; j++) {
        table[0][j] = 0.0;
    }

    /*
      Setting average value in middle tiles
    */
    double average = 0.0;
    for (int i = 1; i < M - 1; i++) {
        average = average + table[i][0] + table[i][N - 1];
    }
    for (int j = 0; j < N; j++) {
        average = average + table[M - 1][j] + table[0][j];
    }
    average = average / (double)(2 * M + 2 * N - 4);

    for (int i = 1; i < M - 1; i++) {
        for (int j = 1; j < N - 1; j++) {
            table[i][j] = average;
        }
    }

    // print2DArray(table);

    /*
      Algorithm start
    */
    int iterations = 0;
    double diff = EPSILON;

    while (EPSILON <= diff) {

        // Copy table to old table
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                tableOld[i][j] = table[i][j];
            }
        }

        // Calculate new values
        for (int i = 1; i < M - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                table[i][j] = (tableOld[i - 1][j] + tableOld[i + 1][j] + tableOld[i][j - 1] + tableOld[i][j + 1]) / 4.0;
            }
        }

        diff = 0.0;
        // 
        for (int i = 1; i < M - 1; i++) {
            for (int j = 1; j < N - 1; j++) {

                double localDiff = fabs(table[i][j] - tableOld[i][j]);
                if (diff < localDiff) {
                    // std::cout << "Changing diff from " << diff << ", to " << localDiff << std::endl;
                    diff = localDiff;
                }
            }
        }

        iterations++;
        // std::cout << "Diff: " << diff << ", Eps: " << EPSILON << std::endl;
    }

    std::cout << "Number of iterations: " << iterations << std::endl;

    print2DArray(table, M, N);

    return 0;
}
