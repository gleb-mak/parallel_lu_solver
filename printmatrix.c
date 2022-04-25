#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void print_matrix(double **cornerArray, double *buf, int n, int l, int m, int size, int rank)
{
	int i, j;
	int ll = l, nn = n;
   	if (ll > m) ll = m;
   	if (nn > m) nn = m;
	for (i = 0; i < nn; i++)
	{
		if (i % size == rank) // send to root process left part of line
			MPI_Send(cornerArray[i / size], i + 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		if (rank == 0)
			MPI_Recv(buf, i + 1, MPI_DOUBLE, i % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for (j = i + 1; j < nn; j++)
		{
			if (rank == j % size) //send the rest part
				MPI_Send(cornerArray[j / size] + 2 * j - i, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
			if (rank == 0)
				MPI_Recv(buf + j, 1, MPI_DOUBLE, j % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		if (rank == 0)
		{ 
			for (j = 0; j < ll; j++) // root process print current line
			{
				printf("%10.3e ", buf[j]);
			}
			printf("\n");
		}
	}
}

void print_vector(double *x, double *buf, int n, int m, int size, int rank)
{
	int i, nn;
	if (m < n)	n = m;
	if (rank == 0)
		printf("\n");
	for (i = 0; i < n; i++)
	{
		if (rank == i % size)
			MPI_Send(x + i / size, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		if (rank == 0)
		{
			MPI_Recv(buf, 1, MPI_DOUBLE, i % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("%10.3e\n", buf[0]);
		}
	}
}
