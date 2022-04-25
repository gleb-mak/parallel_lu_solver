#include <stdio.h>
#include <mpi.h>
#include "matrixfile.h"

int matrix_from_file(double **cornerArray, int n, char *file_name, int size, int rank, double *buf) //input from file
{
	int i = 0;
	int j = 0;
	FILE* ff;
	if (rank == 0)
	{
		ff = fopen(file_name, "r");
		if (ff == NULL)
		{
			printf("\nFile %s was not opened\n", file_name);
			return 1;
		}
	}
	for (i = 0; i < n; i++) 
	{
		if (rank == 0)
		{
			j = 0;
			while (j < n && !feof(ff) && fscanf(ff, "%lf", &buf[j]) == 1) // read current line
				j++;
			if (j < n)
			{
				printf("\nNot enough data in the file, j = %d\n", j);
				return 2;
			}
			MPI_Send(buf, i + 1, MPI_DOUBLE, i % size, 1, MPI_COMM_WORLD); //send left part of line to (i % size)th process
		}
		if (rank == i % size)
			MPI_Recv(cornerArray[i / size], i + 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for (j = i + 1; j < n; j++) // send the rest line 
		{
			if (rank == 0)
				MPI_Send(buf + j, 1, MPI_DOUBLE, j % size, 1, MPI_COMM_WORLD);
			if (rank == j % size)
				MPI_Recv(cornerArray[j / size] + 2 * j - i, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	if (rank == 0)
		fclose(ff);
    return 0;
}


