#include <iostream>
#include <stdlib.h>
#include <math.h>

#define M_SIZE 20
#define N_SIZE 20
const double EPSILON = 0.001;

double table[M_SIZE][N_SIZE];
double tableOld[M_SIZE][N_SIZE];

void print2DArray(double table[M_SIZE][N_SIZE]) {
	std::cout << "[" << std::endl;
	for (int i = 0; i < M_SIZE; i++) {
		std::cout << "[";
		for (int j = 0; j < N_SIZE; j++) {
			if (j != 0) {
				std::cout << ", ";
			}
			std::cout << table[i][j];
		}
		std::cout << "]" << std::endl;
	}
	std::cout << "]" << std::endl;
}


int main() {

	/*
	  Initializing table borders
	*/
	for (int i = 1; i < M_SIZE - 1; i++) {
		table[i][0] = 100.0;
	}
	for (int i = 1; i < M_SIZE - 1; i++) {
		table[i][N_SIZE - 1] = 100.0;
	}
	for (int j = 0; j < N_SIZE; j++) {
		table[M_SIZE - 1][j] = 100.0;
	}
	for (int j = 0; j < N_SIZE; j++) {
		table[0][j] = 0.0;
	}

	/*
	  Setting average value in middle tiles
	*/
	double average = 0.0;
	for (int i = 1; i < M_SIZE - 1; i++) {
		average = average + table[i][0] + table[i][N_SIZE - 1];
	}
	for (int j = 0; j < N_SIZE; j++) {
		average = average + table[M_SIZE - 1][j] + table[0][j];
	}
	average = average / (double)(2 * M_SIZE + 2 * N_SIZE - 4);

	for (int i = 1; i < M_SIZE - 1; i++) {
		for (int j = 1; j < N_SIZE - 1; j++) {
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
		for (int i = 0; i < M_SIZE; i++) {
			for (int j = 0; j < N_SIZE; j++) {
				tableOld[i][j] = table[i][j];
			}
		}

		// Calculate new values
		for (int i = 1; i < M_SIZE - 1; i++) {
			for (int j = 1; j < N_SIZE - 1; j++) {
				table[i][j] = (tableOld[i - 1][j] + tableOld[i + 1][j] + tableOld[i][j - 1] + tableOld[i][j + 1]) / 4.0;
			}
		}

		diff = 0.0;
		for (int i = 1; i < M_SIZE - 1; i++) {
			for (int j = 1; j < N_SIZE - 1; j++) {
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

	print2DArray(table);

	return 0;
}
