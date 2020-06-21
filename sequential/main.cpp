#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <chrono>

#define MAXITERS 2000


double ** initTable(int m, int n) {
    double ** table = (double **)malloc(m * sizeof(double *));
    for (int i = 0; i < m; i++) {
        table[i] = (double *)malloc(n * sizeof(double));
    }
    return table;
}

void updateTables(double*** b, double*** bn) {
	double** bt;
	bt = *b;
	*b = *bn;
	*bn = bt;
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

void initArgs(int argc, const char * argv[], int * table_w, int * table_h, int * tile_w, int * tile_h) {
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

int main(int argc, const char * argv[]) {
    int W, H, w, h;
    int iters = 0;
    initArgs(argc, argv, &W, &H, &w, &h);
    int M = H / h;
    int N = W / w;

    double ** table = initTable(M, N);
    double ** tableNew = initTable(M, N);

    /*
      Initializing table borders
    */
    for (int i = 1; i < M - 1; i++) {
        table[i][0] = 100.0;
        table[i][N - 1] = 100.0;

        tableNew[i][0] = 100.0;
        tableNew[i][N - 1] = 100.0;
    }
    for (int j = 0; j < N; j++) {
        table[M - 1][j] = 100.0;
        table[0][j] = 0.0;

        tableNew[M - 1][j] = 100.0;
        tableNew[0][j] = 0.0;
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
            tableNew[i][j] = average;
        }
    }

    /*
      Algorithm start
    */
    auto start_time = std::chrono::high_resolution_clock::now();
    while (iters < MAXITERS) {
        for (int i = 1; i < M - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                tableNew[i][j] = 0.5 * ((table[i + 1][j] + table[i - 1][j]) / (1 + pow(w / h, 2.0)) + ((table[i][j + 1] + table[i][j - 1]) / (1 + pow(h / w, 2.0))));
            }
        }

        iters++;
        updateTables(&table, &tableNew);
        
        /*if (iters > 0 && iters % 10 == 0) {
            printf("Completed iteration %d/%d!\n", iters, MAXITERS);
        }*/
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::cout << "Algorithm finished!" << std::endl;

    auto runningTime = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "Execution time: " << runningTime << " ms." << std::endl;

    return 0;
}
