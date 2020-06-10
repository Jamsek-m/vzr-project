#include <string.h>
#include "/usr/include/openmpi-x86_64/mpi.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define MAXITERS 100

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
    MPI_Request request;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    int * finished;
    bool finish = false;

	if (myid == 0){
        finished = (int*)malloc(procs * sizeof(int));
        // Initialize table
        board = board_initialize(M, N);
        boardptr = *board;

        borders_initialize(board, M, N);
        board_print(board, M, N);
    }

    
    // divide work
	mystart = M / procs * myid;				// determine scope of work for each process; process 0 also works on its own part
	myend = M / procs * (myid + 1);
	myrows = M / procs;

    myboard = board_initialize(myrows, N);
	myboard_new = board_initialize(myrows, N);
    myrow_top = (double*)malloc(N * sizeof(double));
	myrow_bot = (double*)malloc(N * sizeof(double));
    

    // scatter initial matrix
	MPI_Scatter(boardptr, myrows * N, MPI_DOUBLE, 
				*myboard, myrows * N, MPI_DOUBLE, 
				0, MPI_COMM_WORLD);

    if (myid == 0){
        for (i = 0; i < myrows; i++){
            for (j = 0; j < N; j++){
                myboard_new[i][j] = 0.0;
            }
        }
    } else if (myid == (procs-1)){
        for (i = 0; i < myrows; i++){
            for (j = 0; j < N; j++){
                myboard_new[i][j] = 100.0;
            }
        }
    }
    int message = myid;
    int messageReceived = -1;

    // MPI_Irecv(&message, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &request);

    // do the calculation
	while (iters < MAXITERS && finish == false)
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
            myboard_new[i][0] = 100.0;
            myboard_new[i][N-1] = 100.0;
            for (j = 1; j < (N-1); j++){
                if(myrows == 1){

                    if (myid > 0 && myid < (procs - 1)){
                        temperature = 0.5 * ((myrow_bot[j] + myrow_top[j])/(1.0 + (double)((w/h)*(w/h))) + (myboard[i][j-1] + myboard[i][j+1])/(1.0 + (double)((h/w)*(h/w))));
                            // temperature = (double) (0.25 * (myrow_bot[j] + myrow_top[j] + myboard[i][j-1] + myboard[i][j+1]));
                        // if(myid == (procs-2)){
                        //     printf("\nmyid: %d, temperature: %6.2f\n", myid, temperature);
                        //     printf("\nmyrow_bot[j]: %6.2f, myrow_top[j]: %6.2f, myboard[i][j-1]: %6.2f, myboard[i][j+1]: %6.2f\n", myrow_bot[j],myrow_top[j],myboard[i][j-1],myboard[i][j+1]);
                        // }
                        myboard_new[i][j] = temperature;
                    }
                }else if(i == 0 && myid != 0){
                    temperature = 0.5 * ((myboard[i+1][j] + myrow_top[j])/(1.0 + (double)((w/h)*(w/h))) + (myboard[i][j-1] + myboard[i][j+1])/(1.0 + (double)((h/w)*(h/w))));
                }else if(i == (myrows - 1) && myid != (procs-1)){
                    temperature = 0.5 * ((myrow_bot[j] + myboard[i-1][j])/(1.0 + (double)((w/h)*(w/h))) + (myboard[i][j-1] + myboard[i][j+1])/(1.0 + (double)((h/w)*(h/w))));
                }
                // myboard_new[i][j] = temperature;
            }
        }
        iters++;
        
        if(iters % 10 == 0){
            diff = 0.0;
            localDiff = 0.0;

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
                if(diff <= EPSILON){
                    printf("\nmyid: %d, iteration: %d, diff: %6.2f\n", myid, iters, (double)diff);
                    
                    // MPI_Test(&request, &messageReceived, MPI_STATUS_IGNORE);
                    // if (messageReceived) {
                    //     if(myid == 0){
                    //         printf("root found that process %d finished\n", message);
                    //         finished[message] = 1;
                    //         bool all = true;
                    //         for(int k = 1; k < procs; k++){
                    //             if(finished[k] != 1){
                    //                 all=false;
                    //                 break;
                    //             }
                    //         }
                    //         if (all == true){
                    //             MPI_Isend(message, 1, MPI_INT, myid, 0, MPI_COMM_WORLD, &request);
                    //             finish = true;
                    //         }
                    //     }else {
                    //         finish = true;
                    //     }
                    // }else {
                    //     MPI_Isend(message, 1, MPI_INT, myid, 0, MPI_COMM_WORLD, &request);
                    // }
                }
            }
        }
        
        // if(myid == 14){
        //     printf("\nmyid: %d, iteration: %d, diff: %6.2f\n", myid, iters, (double)diff);
        //     printf("\nprint2DArray myboard_new start\n");
        //     board_print(myboard_new, myrows, N);
        //     board_print(myboard, myrows, N);
        //     printf("\nprint2DArray myboard_new success\n");
        // }
		// swap boards (iter --> iter + 1)
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
        board_print(board, M, N);
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