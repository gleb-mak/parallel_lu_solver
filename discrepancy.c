#include <math.h>
#include <mpi.h>
#include "discrepancy.h"

double disc(double **cornerArray, double *b, double *x, double *buf, int n, int k, int size, int rank)
{
	double s = 0;
	double t = 0;
	double r = 0;
	int i, j;
	double bi;
	for (i = 0; i < n; i++) // send to root process vector x
	{
		if (i % size == rank)
			MPI_Send(x + i / size, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		if (rank == 0)
			MPI_Recv(buf + n + i, 1, MPI_DOUBLE, i % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	}
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++) //compute row i (write in buf)
		{
			if (k == 1)
				buf[j] = n - fmax(i, j);
			if (k == 2)
				buf[j] = fmax(i, j);
			if (k == 3)
				buf[j] = fabs(i - j);
			if (k == 4)
				buf[j] = (double)(1. / (i+j+1));
		}
		bi = 0;
		for (j = 0; j < n; j += 2) // compute b[i]
			bi += buf[j];
		if (rank == 0) // root process compute discrepancy
		{
			t += bi * bi;
        	s = bi;
        	for (j = 0; j < n; j++)
			{
           		s -= buf[j] * buf[n + j];
			}
       		r += s * s;
		}		
	}
	t = sqrt(t);
    r = sqrt(r);
	if (rank == 0)    return(r/t);
	else	return 1;
}

double solution_err(double* x, double *buf, int n, int size, int rank)
{
    double s = 0;
	int i;
    for (i = 0; i < n; i++)
	{
		if (rank == i % size) // send x to root process
			MPI_Send(x + i / size, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		if (rank == 0) // root process compute error
		{
			MPI_Recv(buf, 1, MPI_DOUBLE, i % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        	s += ((1-i%2) - buf[0])*((1-i%2)-buf[0]);
		}
	}
    return(sqrt(s));
}
