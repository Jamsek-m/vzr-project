#include <string.h>
#include "/usr/include/openmpi-x86_64/mpi.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <chrono>

#define MAXITERS 100

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
    }
}

void board_update(double*** b, double*** bn)
{
	double** bt;

	bt = *b;
	*b = *bn;
	*bn = bt;
}

void board_free(double** b)
{
	free(*b);
	free(b);
}

double** board_initialize(int n, int m)
{
	int k, l;

	double* bd = (double*)malloc(sizeof(double) * n * m);
	double** b = (double**)malloc(sizeof(double*) * n);
	for (k = 0; k < n; k++)
		b[k] = &bd[k * m];

	for (k = 0; k < n; k++)
		for (l = 0; l < m; l++)
			// b[k][l] = (double)((k * m) + l);
			b[k][l] = 0.0;
			// b[k][l] = 73.33;
			// b[k][l] = rand() < 0.25 * RAND_MAX;

	return b;
}

void board_print(double** b, int n, int m)
{
	int k, l;

	for (k = 0; k < n; k++)
	{
		for (l = 0; l < m; l++)
			printf("%6.2f ", b[k][l]);
		printf("\n");
	}
	printf("\n");
}

void borders_initialize(double** board, int m, int n)
{
	for (int i = 1; i < m - 1; i++) {
        board[i][0] = 100.0;
        board[i][n - 1] = 100.0;
    }
    for (int j = 0; j < n; j++) {
        board[m - 1][j] = 100.0;
        board[0][j] = 0.0;
    }

    // Setting average value in middle tiles
    double average = 0.0;
    for (int i = 1; i < m - 1; i++) {
        average = average + board[i][0] + board[i][n - 1];
    }
    for (int j = 0; j < n; j++) {
        average = average + board[m - 1][j] + board[0][j];
    }
    average = average / (double)(2 * m + 2 * n - 4);

    for (int i = 1; i < m - 1; i++) {
        for (int j = 1; j < n - 1; j++) {
            board[i][j] = average;
        }
    }
}

double calculate_temperature(int i, int j, int h, int w, double **myboard, double *myrow_top, double *myrow_bot, double *mycolumn_right, double *mycolumn_left)
{
    double wh = 1.0 + (double)((w/h)*(w/h));
    double hw = 1.0 + (double)((w/h)*(w/h));
    double temperature = 0;
	
    // gor desno
    if (j == w - 1 && i == 0){
        temperature = 0.5 * ((myboard[i + 1][j] + myrow_top[i])/wh + (myboard[i][j - 1] + mycolumn_right[i])/hw);
    // gor vmes
    } else if (j != 0 && j != w - 1 && i == 0){
        temperature = 0.5 * ((myboard[i + 1][j] + myrow_top[i])/wh + (myboard[i][j - 1] + myboard[i][j + 1])/hw);
    // gor levo
    } else if (j == 0 && i == 0){
        temperature = 0.5 * ((myboard[i + 1][j] + myrow_top[i])/wh + (myboard[i][j + 1] + mycolumn_left[i])/hw);
    // dol desno
    } else if (j == w - 1 && i == h - 1){
        temperature = 0.5 * ((myboard[i - 1][j] + myrow_bot[i])/wh + (myboard[i][j - 1] + mycolumn_right[i])/hw);
    // dol vmes
    } else if (j != 0 && j != w - 1 && i == h - 1){
        temperature = 0.5 * ((myboard[i - 1][j] + myrow_bot[i])/wh + (myboard[i][j - 1] + myboard[i][j + 1])/hw);
    // dol levo
    } else if (j == 0 && i == h - 1){
        temperature = 0.5 * ((myboard[i - 1][j] + myrow_bot[i])/wh + (myboard[i][j + 1] + mycolumn_left[i])/hw);
    // levo vmes
    } else if (j == 0 && i != h - 1 && i != 0){
        temperature = 0.5 * ((myboard[i + 1][j] + myboard[i - 1][j])/wh + (myboard[i][j + 1] + mycolumn_left[i])/hw);
    // desno vmes
    } else if (j == w - 1 && i != h - 1 && i != 0){
        temperature = 0.5 * ((myboard[i + 1][j] + myboard[i - 1][j])/wh + (myboard[i][j - 1] + mycolumn_right[i])/hw);
    // vmesne
    } else {
        temperature = 0.5 * ((myboard[i + 1][j] + myboard[i - 1][j])/wh + (myboard[i][j - 1] + myboard[i][j + 1])/hw);
    }

    return temperature;
}



int main(int argc, char** argv)
{
    int W, H, w, h;
    initArgs(argc, argv, &W, &H, &w, &h);
    int M = H / h;
    int N = W / w;

    double * boardptr = NULL;
	double ** board;

    int myid, procs;
    int mystart, myend, myrows;
    double ** myboard;
	double ** myboard_new;
    double * myrow_top, * myrow_bot, * mycolumn_left, * mycolumn_right;

    int i, j;
    double temperature;
    int iters = 0;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

	if (myid == 0){
        // Initialize table
        board = board_initialize(H, W);
        boardptr = *board;

        borders_initialize(board, H, W);
        board_print(board, H, W);
        printf("\nW: %d, H: %d, w: %d, h: %d, M: %d, N: %d, procs: %d\n", W,H,w,h,M,N,procs);
    }

    myboard = board_initialize(h, w);
	myboard_new = board_initialize(h, w);
    myrow_top = (double*)malloc(w * sizeof(double));
	myrow_bot = (double*)malloc(w * sizeof(double));
    mycolumn_left = (double*)malloc(h * sizeof(double));
	mycolumn_right = (double*)malloc(h * sizeof(double));

    // MPI_Datatype rowtype;
    MPI_Datatype columntype;
    MPI_Datatype blocktype;
    MPI_Datatype blocktype_resized;

    MPI_Type_vector(h, 1, w, MPI_DOUBLE, &columntype);
    MPI_Type_commit(&columntype);

    MPI_Type_vector(h, w, W, MPI_DOUBLE, &blocktype);
    MPI_Type_commit(&blocktype);
    MPI_Type_create_resized(blocktype, 0, w * sizeof(double), &blocktype_resized);
    MPI_Type_commit(&blocktype_resized);

    int counts[procs];
    int displacements[procs];
    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            counts[i*N + j] = 1;
            displacements[i*N + j] = j + (i * h * N) ;
        }
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();
    MPI_Scatterv(boardptr, counts, displacements, blocktype_resized, 
				*myboard, w*h, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // printf("Scatter is ok for proc %d\n", myid);
    // board_print(myboard, h, w);

    for (i = 0; i < h; i++){
        for (j = 0; j < w; j++){
            myboard_new[i][j] = myboard[i][j];
        }
    }

    // do the calculation
	while (iters < MAXITERS)
	{
        // ce nisi zadnja vrstica poslji dol in sprejmi od dol
        if (myid < (M - 1) * N){
            MPI_Sendrecv(myboard[h - 1], w, MPI_DOUBLE, myid + N, 1,
                        myrow_bot, w, MPI_DOUBLE, myid + N, 0, 
                        MPI_COMM_WORLD, &status);
        }
        // ce nisi 1.vrstica poslji gor in sprejmi od gor
        if (myid >= N){
            MPI_Sendrecv(myboard[0], w, MPI_DOUBLE, myid - N, 0,
                        myrow_top, w, MPI_DOUBLE, myid - N, 1, 
                        MPI_COMM_WORLD, &status);
        }    
        // ce nisi najbolj levi stolpec poslji levo in sprejmi od levega
        if (myid % N != 0){
            MPI_Sendrecv(&myboard[0][0], 1, columntype, myid - 1, 2,
                        mycolumn_left, w, MPI_DOUBLE, myid - 1, 3, 
                        MPI_COMM_WORLD, &status);
        }
        // ce nisi najbolj desni stolpec poslji desno in sprejmi od desnega
        if (myid % N != N - 1){
            MPI_Sendrecv(&myboard[0][w - 1], 1, columntype, myid + 1, 3,
                        mycolumn_right, w, MPI_DOUBLE, myid + 1, 2,
                        MPI_COMM_WORLD, &status);
        } 

		// do the computation of my part
        for (i = 0; i < h; i++){
            for (j = 0; j < w; j++){
                temperature = myboard[i][j];
                // 1.vrstica
                if (myid <  N){
                    // 1.stolpec
                    if (myid % N == 0 && i != 0 && j != 0){
                        temperature = calculate_temperature(i, j, h, w, myboard, myrow_top, myrow_bot, mycolumn_right, mycolumn_left);
                    // zadnji stolpec
                    } else if (myid % N == N - 1 && i != 0 && j != w - 1){
                        temperature = calculate_temperature(i, j, h, w, myboard, myrow_top, myrow_bot, mycolumn_right, mycolumn_left);
                    // vmesni stolpci
                    } else if (myid % N != 0 && myid % N != N - 1 && i != 0){
                        temperature = calculate_temperature(i, j, h, w, myboard, myrow_top, myrow_bot, mycolumn_right, mycolumn_left);
                    }

                // zadnja vrstica
                } else if (myid >= (M - 1) * N){
                    // 1.stolpec
                    if (myid % N == 0 && i != h-1 && j != 0){
                        temperature = calculate_temperature(i, j, h, w, myboard, myrow_top, myrow_bot, mycolumn_right, mycolumn_left);
                    // zadnji stolpec
                    } else if (myid % N == N - 1 && i != h-1 && j != w - 1){
                        temperature = calculate_temperature(i, j, h, w, myboard, myrow_top, myrow_bot, mycolumn_right, mycolumn_left);
                    // vmesni stolpci
                    } else if (myid % N != 0 && myid % N != N - 1 && i != h-1){
                        temperature = calculate_temperature(i, j, h, w, myboard, myrow_top, myrow_bot, mycolumn_right, mycolumn_left);
                    } 

                // vmesna vrstica
                } else {
                    // 1.stolpec
                    if (myid % N == 0 && j != 0){
                        temperature = calculate_temperature(i, j, h, w, myboard, myrow_top, myrow_bot, mycolumn_right, mycolumn_left);
                    // zadnji stolpec
                    } else if (myid % N == N - 1 && j != w - 1){
                        temperature = calculate_temperature(i, j, h, w, myboard, myrow_top, myrow_bot, mycolumn_right, mycolumn_left);
                    // vmesni stolpci
                    } else if (myid % N != 0 && myid % N != N - 1){
                        temperature = calculate_temperature(i, j, h, w, myboard, myrow_top, myrow_bot, mycolumn_right, mycolumn_left);
                    } 
                }
                myboard_new[i][j] = temperature;
            }
        }
        iters++;
        
		board_update(&myboard, &myboard_new);
	}

    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            counts[i*N + j] = w*h;
            displacements[i*N + j] = j + (i * h * N);
        }
    }

    MPI_Gatherv(*myboard, w*h, MPI_DOUBLE, boardptr,
                counts, displacements, blocktype_resized, 0, MPI_COMM_WORLD);
    
	// display
	if (myid == 0) {
        printf("print board\n");
        board_print(board, H, W);
        auto end_time = std::chrono::high_resolution_clock::now();
        auto runningTime = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        std::cout << "Execution time: " << runningTime << " ms." << std::endl;
    }

    // free memory
	if (myid == 0){
		board_free(board);
    }
    board_free(myboard);
    free(myrow_top);
    free(myrow_bot);

	MPI_Finalize();			// finalize MPI

	return 0;
}
