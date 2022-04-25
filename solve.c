#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "solve.h"
#include "printmatrix.h"

int solve_gauss(double **cornerArray, double *buf, double *b, double *x, int n, int size, int rank)
{
	int i, j, k;

	int numOfCorners = (rank >= n % size) ? n / size : 1 + n / size; //check applicability
	for (j = 0; j < numOfCorners; j++)	
	{
		buf[0] = 0; 
		if(fabs(cornerArray[j][j * size + rank]) < 1e-30)
		{
			buf[0] = 1;
		}
		MPI_Send(buf, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}
	if (rank == 0)
		for (j = 0; j < n; j++)
		{
			MPI_Recv(buf, 1, MPI_DOUBLE, j % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (buf[0] == 1)
			{
				printf("\nCan not use the Gauss method\n");
				MPI_Abort(MPI_COMM_WORLD, 6);
				return 6;
			}
		}
	
	for (i = 0; i < n; i++) //LU
	{
			int cornerIndex = i / size;
			if (rank == i % size)
			{
				for (k = 0; k < i; k++)
					cornerArray[cornerIndex][i] -= cornerArray[cornerIndex][k] * cornerArray[cornerIndex][2 * i - k];
				for (j = 0; j < 2 * i + 1; j++)
					buf[j] = cornerArray[cornerIndex][j];
			}
			MPI_Bcast(buf, 2 * i + 1, MPI_DOUBLE, i % size, MPI_COMM_WORLD);
			for (j = i + 1; j < n; j++)
			{
				if (rank == j % size)
				{
					cornerIndex = j / size;
					for (k = 0; k < i; k++)
					{
						cornerArray[cornerIndex][i] -= cornerArray[cornerIndex][k] * buf[2 * i - k];	
						cornerArray[cornerIndex][2 * j - i] -= buf[k] * cornerArray[cornerIndex][2 * j - k];
					}
					cornerArray[cornerIndex][2 * j - i] /= buf[i];
				}
			}
	}

	for (i = 0; i < n; i++) //Ly = b, y in x
	{
		if (rank == i % size)
		{
			x[i / size] = b[i / size] / cornerArray[i / size][i];
			buf[0] = x[i / size];
		}
		MPI_Bcast(buf, 1, MPI_DOUBLE, i % size, MPI_COMM_WORLD);
		for (j = i + 1; j < n; j++)
		{
			if (rank == j % size)
				b[j / size] -= buf[0] * cornerArray[j / size][i];
		}
	}

	for (i = n - 1; i >= 0; i--) //Ux = y, y in x
	{
			for (j = i - 1; j >= 0; j--)
			{
				if (rank == i % size)
				{
					buf[0] = cornerArray[i / size][2 * i - j] * x[i / size];
					MPI_Send(buf, 1, MPI_DOUBLE, j % size, 1, MPI_COMM_WORLD);
				}
				if (rank == j % size)
				{
					MPI_Recv(buf, 1, MPI_DOUBLE, i % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					x[j / size] -= buf[0];
				}
			}
	}
	return 0;
}
