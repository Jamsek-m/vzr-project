#include <string.h>
#include "/usr/include/openmpi-x86_64/mpi.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define MAXITERS 100


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

void initArgs(int argc, char ** argv, int * table_w, int * table_h, int * tile_w, int * tile_h, double * eps) {
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

int main(int argc, char** argv)
{
    int W, H, w, h;
    double EPSILON;
    initArgs(argc, argv, &W, &H, &w, &h, &EPSILON);
    int M = H / h;
    int N = W / w;
    double localDiff = 0;
    double diff = 1.0;

    double * boardptr = NULL;
	double ** board;

    int myid, procs;
    int mystart, myend, myrows;
    double ** myboard;					
	double ** myboard_new;
    double * myrow_top, * myrow_bot;

    int i, j;
    double temperature;
    int iters = 0;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

	if (myid == 0){
        // Initialize table
        board = initTable(M, N);
        boardptr = *board;

        // Initializing table borders
        for (int i = 1; i < M - 1; i++) {
            board[i][0] = 100.0;
            board[i][N - 1] = 100.0;
        }
        for (int j = 0; j < N; j++) {
            board[M - 1][j] = 100.0;
            board[0][j] = 0.0;
        }
 
        // Setting average value in middle tiles
        double average = 0.0;
        for (int i = 1; i < M - 1; i++) {
            average = average + board[i][0] + board[i][N - 1];
        }
        for (int j = 0; j < N; j++) {
            average = average + board[M - 1][j] + board[0][j];
        }
        average = average / (double)(2 * M + 2 * N - 4);

        for (int i = 1; i < M - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                board[i][j] = average;
            }
        }
        print2DArray(board, M, N);
    }

    
    // divide work
	mystart = M / procs * myid;				// determine scope of work for each process; process 0 also works on its own part
	myend = M / procs * (myid + 1);
	myrows = M / procs;

    printf("\nprocs: %d, myrows: %d", procs, myrows);

    myboard = initTable(myrows, N);
	myboard_new = initTable(myrows, N);
    myrow_top = (double*)malloc(N * sizeof(double));
	myrow_bot = (double*)malloc(N * sizeof(double));
    

    // scatter initial matrix
	MPI_Scatter(boardptr, myrows * N, MPI_DOUBLE, 
				*myboard, myrows * N, MPI_DOUBLE, 
				0, MPI_COMM_WORLD);

    // do the calculation
	while (EPSILON <= diff && iters < MAXITERS)
	{
        if (myid > 0 && myid < (procs - 1)){
            // exchange borders with neigbouring processes
            MPI_Sendrecv(myboard[0], N, MPI_DOUBLE, (myid + procs - 1) % procs, 1,
                        myrow_bot, N, MPI_DOUBLE, (myid + 1) % procs, 1,
                        MPI_COMM_WORLD, &status);//MPI_STATUSES_IGNORE
            // ptr to send data, send data size, send data type, receiver, message tag,
            // ptr to received data, received data size, recevied data type, sender, message tag,
            // communicator, status
            MPI_Sendrecv(myboard[myrows - 1], N, MPI_DOUBLE, (myid + 1) % procs, 0,
                        myrow_top, N, MPI_DOUBLE, (myid + procs - 1) % procs, 0, 
                        MPI_COMM_WORLD, &status);
        }
		
        if (myid == 0 && myid < (procs - 1)){
            MPI_Sendrecv(myboard[myrows - 1], N, MPI_DOUBLE, (myid + 1) % procs, 0,
                        myrow_bot, N, MPI_DOUBLE, (myid + 1) % procs, 1, 
                        MPI_COMM_WORLD, &status);
        }    

        if (myid > 0 && myid == (procs - 1)){
            MPI_Sendrecv(myboard[0], N, MPI_DOUBLE, (myid - 1) % procs, 1,
                        myrow_bot, N, MPI_DOUBLE, (myid - 1) % procs, 0,
                        MPI_COMM_WORLD, &status);
        }   

		// do the computation of my part
        if (myid > 0 && myid < (procs - 1)){
            for (i = 0; i < myrows; i++){
                for (j = 1; j < (N-1); j++){
                    if(myrows == 1){
                        // temperature = 0.5 * ((myrow_bot[j] + myrow_top[j])/(1.0 + (double)((w/h)*(w/h))) + (myboard[i][j-1] + myboard[i][j+1])/(1.0 + (double)((h/w)*(h/w))));
                        temperature = (double) (0.25 * (myrow_bot[j] + myrow_top[j] + myboard[i][j-1] + myboard[i][j+1]));
                    }else if(i == 0){
                        temperature = 0.5 * ((myboard[i+1][j] + myrow_top[j])/(1.0 + (double)((w/h)*(w/h))) + (myboard[i][j-1] + myboard[i][j+1])/(1.0 + (double)((h/w)*(h/w))));
                    }else if(i == (myrows - 1)){
                        temperature = 0.5 * ((myrow_bot[j] + myboard[i-1][j])/(1.0 + (double)((w/h)*(w/h))) + (myboard[i][j-1] + myboard[i][j+1])/(1.0 + (double)((h/w)*(h/w))));
                    }
                    myboard_new[i][j] = temperature;
                }
            }
        }
        // printf("\nprint2DArray myboard\n");
        // print2DArray(myboard, 1, N);
        // printf("\nprint2DArray myboard_new\n");
        // print2DArray(myboard_new, 1, N);

        
        for (i = 0; i < myrows; i++){
			for (j = 1; j < (N-1); j++)
			{
				double calculatedDiff = fabs(myboard[i][j] - myboard_new[i][j]);
                if (localDiff < calculatedDiff) {
                    localDiff = calculatedDiff;
                }
			}
        }
        if (diff < localDiff) {
            diff = localDiff;
        }

		iters++;
		// swap boards (iter --> iter + 1)
        printf("\nmyid: %d, iteration: %d, diff: %2.0f\n", myid, iters, diff);
		// board_update(&myboard, &myboard_new);

        for (i = 0; i < myrows; i++){
			for (j = 0; j < N; j++){
                myboard[i][j] = myboard_new[i][j];
            }
        }

        printf("\nprint2DArray myboard_new start\n");
        print2DArray(myboard_new, myrows, N);
        printf("\nprint2DArray myboard_new success\n");
	}
	
	// gather results
	MPI_Gather(*myboard, myrows * N, MPI_DOUBLE, 
			   boardptr, myrows * N, MPI_DOUBLE, 
			   0, MPI_COMM_WORLD);
	// data to send, send data size, data type,
	// gathered data, recevied data size, data type,
	// gathering process, communicator

	// display
	if (myid == 0) {
        printf("\nprint2DArray start\n");
        print2DArray(board, M, N);
        printf("\nprint2DArray success\n");
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
