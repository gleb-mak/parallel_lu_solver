#include <stdio.h>
#include <mpi.h>
#include "formb.h"

void form_b(double **cornerArray, double *b, double *buf, int n, int size, int rank) // system right part
{
    int i, j;
    for (i = 0; i < n; i++)
    {
		if (i % size == rank) // compute left part of row (from 0 to i)
		{
     		b[i / size] = 0;
        	for (j = 0; j <= i; j += 2)
			{
            	b[i / size] += cornerArray[i / size][j];
			}
		}
		j = (i % 2 == 0) ? i + 2 : i + 1;
		for (j; j < n; j += 2) // compute the rest part
		{
			if (j % size == rank)
				MPI_Send(cornerArray[j / size] + 2 * j - i, 1, MPI_DOUBLE, i % size, 1, MPI_COMM_WORLD);
			if (i % size == rank)
			{
				MPI_Recv(buf, 1, MPI_DOUBLE, j % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				b[i / size] += buf[0];
			}
		}
    }
}
