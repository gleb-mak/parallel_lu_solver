#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrixfile.h"
#include "printmatrix.h"
#include "matrixformula.h"
#include "formb.h"
#include "solve.h"
#include "discrepancy.h"

int main(int argc, char *argv[]) {
	int n = atoi(argv[1]);
	int m = atoi(argv[2]);
	int k = atoi(argv[3]);
	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status Status;
	
	int error = 0;
	int numOfCorners = (rank >= n % size) ? n / size : 1 + n / size;
	double **cornerArray;
	double *buf, *b, *x;
	double r;
	cornerArray = (double**)malloc(numOfCorners * sizeof(double*));
	if (!cornerArray)
	{
		printf("Can not allocate memory\n");
        error = 4;
	}
	int i;
	for (i = 0; i < numOfCorners; i++)
	{
		cornerArray[i] = (double*)malloc((2 * rank + 2 * size * i + 1) * sizeof(double));
		if (!cornerArray[i])
		{
			printf("Can not allocate memory\n");
			error = 4;
		}
	}
	buf = (double*)malloc(2 * n * sizeof(double));
	if (!buf)
	{
		printf("Can not allocate memory\n");
        error = 4;
	}
	b = (double*)malloc(numOfCorners * sizeof(double));
	if (!b)
	{
		printf("Can not allocate memory\n");
        error = 4;	
	}
	x = (double*)malloc(numOfCorners * sizeof(double));
	if (!x)
	{
		printf("Can not allocate memory\n");
        error = 4;	
	}
	if (error)
	{
		MPI_Finalize();
		return error;
	}
	if (k == 0)
        error = matrix_from_file(cornerArray, n, argv[4], size, rank, buf);
	else
		form_matrix(cornerArray, n, k, size, rank);
	if (error)
	{
		MPI_Finalize();
		return error;
	}
	form_b(cornerArray, b, buf, n, size, rank);
	if (m > 0)
   	{
		if (rank == 0)  printf("\nMatrix a:\n");
        print_matrix(cornerArray, buf, n, n, m, size, rank);
		if (rank == 0)  printf("\nVector b:\n");
		print_vector(b, buf, n, m, size, rank);
   	}
	MPI_Barrier(MPI_COMM_WORLD);
	double start_time = MPI_Wtime();
	error = solve_gauss(cornerArray, buf, b, x, n, size, rank);
	MPI_Barrier(MPI_COMM_WORLD);
	double end_time = MPI_Wtime();
	double time = end_time - start_time;
	if (error)
	{
		MPI_Finalize();
		return error;
	}
	if (m > 0)
   	{
		if (rank == 0)  printf("\nMatrix LU:\n");
        print_matrix(cornerArray, buf, n, n, m, size, rank);
		if (rank == 0)  printf("\nVector x, LUx=b:\n");
		print_vector(x, buf, n, m, size, rank);
   	}
	if (m >= 0)
	{
		if (k != 0)
		{
			r = disc(cornerArray, b, x, buf, n, k, size, rank);
			if (rank == 0)
				printf("\nDiscrepancy: %10.3e\n", r);
		}
		if (rank == 0)	printf("\nTime: %lf sec.\n", time);
	}
	if (m == -1)
	{
		print_vector(x, buf, n, n, size, rank);
	}
	if ((m != -1))
    {
        r = solution_err(x, buf, n, size, rank);
		if (rank == 0)	printf("\nError of solution: %10.3e\n", r);
    }
	for (i = 0; i < numOfCorners; i++)
		free(cornerArray[i]);
	free(cornerArray);
	free(buf);
	free(b);
	free(x);
	MPI_Finalize();
	if (argc < 4)
		return 5;
	return 0;
}

