#include <string.h>
#include "/usr/include/openmpi-x86_64/mpi.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <chrono>

#define MAXITERS 20000

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
			b[k][l] = rand() < 0.25 * RAND_MAX;

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
        board = board_initialize(M, N);
        boardptr = *board;

        borders_initialize(board, M, N);
        // board_print(board, M, N);
    }

    // divide work
	mystart = M / procs * myid;				// determine scope of work for each process; process 0 also works on its own part
	myend = M / procs * (myid + 1);
    myrows = M / procs;

    myboard = board_initialize(myrows, N);
	myboard_new = board_initialize(myrows, N);
    myrow_top = (double*)malloc(N * sizeof(double));
	myrow_bot = (double*)malloc(N * sizeof(double));
    

    // myboard_new proces 0 prva vrstica 0.0
    if (myid == 0){
        for (j = 0; j < N; j++){
            myboard_new[0][j] = 0.0;
        }
    // myboard_new zadnji proces zadnja vrstica 100.0
    } else if (myid == (procs-1)){
        for (j = 0; j < N; j++){
            myboard_new[myrows - 1][j] = 100.0;
        }
    } 
    // myboard_new vsi procesi prvi in zadnji stolpec 100.0
    for (i = 0; i < myrows; i++){
        myboard_new[i][0] = 100.0;
        myboard_new[i][N-1] = 100.0;
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    // scatter initial matrix
	MPI_Scatter(boardptr, myrows * N, MPI_DOUBLE, 
				*myboard, myrows * N, MPI_DOUBLE, 
				0, MPI_COMM_WORLD);

    // do the calculation
	while (iters < MAXITERS)
	{
        // exchange borders with neigbouring processes
        if (myid > 0 && myid < (procs - 1)){
            MPI_Sendrecv(myboard[0], N, MPI_DOUBLE, (myid + procs - 1) % procs, 1,
                        myrow_bot, N, MPI_DOUBLE, (myid + 1 + procs) % procs, 1,
                        MPI_COMM_WORLD, &status);//MPI_STATUSES_IGNORE
            // ptr to send data, send data size, send data type, receiver, message tag,
            // ptr to received data, received data size, recevied data type, sender, message tag,
            // communicator, status
            MPI_Sendrecv(myboard[myrows - 1], N, MPI_DOUBLE, (myid + 1 + procs) % procs, 0,
                        myrow_top, N, MPI_DOUBLE, (myid + procs - 1) % procs, 0, 
                        MPI_COMM_WORLD, &status);
        }
		
        if (myid == 0 && myid < (procs - 1)){
            MPI_Sendrecv(myboard[myrows - 1], N, MPI_DOUBLE, (myid + 1) % procs, 0,
                        myrow_bot, N, MPI_DOUBLE, (myid + 1) % procs, 1, 
                        MPI_COMM_WORLD, &status);
        }    

        if (myid > 0 && myid == (procs - 1)){
            MPI_Sendrecv(myboard[0], N, MPI_DOUBLE, (myid - 1 + procs) % procs, 1,
                        myrow_top, N, MPI_DOUBLE, (myid - 1 + procs) % procs, 0,
                        MPI_COMM_WORLD, &status);
        }

		// do the computation of my part
        
        for (i = 0; i < myrows; i++){
            for (j = 1; j < (N-1); j++){
                //myrows == 1
                if(myrows == 1 && myid > 0 && myid < (procs - 1)){
                    temperature = 0.5 * ((myrow_bot[j] + myrow_top[j])/(1.0 + (double)((w/h)*(w/h))) + (myboard[i][j-1] + myboard[i][j+1])/(1.0 + (double)((h/w)*(h/w))));
                    myboard_new[i][j] = temperature;
                //myrows != 1 prva vrstica, ce je myid==0 prva vrstica ni pomembna
                }else if(myrows != 1 && i == 0 && myid != 0){
                    temperature = 0.5 * ((myboard[i+1][j] + myrow_top[j])/(1.0 + (double)((w/h)*(w/h))) + (myboard[i][j-1] + myboard[i][j+1])/(1.0 + (double)((h/w)*(h/w))));
                    myboard_new[i][j] = temperature;
                //myrows != 1 zadnja vrstica, ce je myid==(procs-1) tadnja vrstica ni pomembna
                }else if(myrows != 1 && i == (myrows - 1) && myid != (procs-1)){
                    temperature = 0.5 * ((myrow_bot[j] + myboard[i-1][j])/(1.0 + (double)((w/h)*(w/h))) + (myboard[i][j-1] + myboard[i][j+1])/(1.0 + (double)((h/w)*(h/w))));
                    myboard_new[i][j] = temperature;
                //myrows != 1 vmesne vrstice
                }else if(myrows != 1 && i > 0 && i < (myrows - 1)){
                    temperature = 0.5 * ((myboard[i+1][j] + myboard[i-1][j])/(1.0 + (double)((w/h)*(w/h))) + (myboard[i][j-1] + myboard[i][j+1])/(1.0 + (double)((h/w)*(h/w))));
                    myboard_new[i][j] = temperature;
                }
            }
        }
        iters++;
        
		board_update(&myboard, &myboard_new);
	}
    MPI_Barrier(MPI_COMM_WORLD);
	
	// gather results
	MPI_Gather(*myboard, myrows * N, MPI_DOUBLE, 
			   boardptr, myrows * N, MPI_DOUBLE, 
			   0, MPI_COMM_WORLD);
	// data to send, send data size, data type,
	// gathered data, recevied data size, data type,
	// gathering process, communicator

	// display
	if (myid == 0) {
        // board_print(board, M, N);
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
