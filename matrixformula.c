#include <stdio.h>
#include <math.h>
#include "matrixformula.h"

void form_matrix(double** cornerArray, int n, int k, int size, int rank) //input from formula
{
    	int i, j;
		if (k == 1)
		{
			for (i = rank; i < n; i+= size)
			{
				for (j = 0; j <= 2 * i; j++)
					cornerArray[i / size][j] = n - i;
			}
		}

		if (k == 2)
		{
			for (i = rank; i < n; i+= size)
			{
				for (j = 0; j <= 2 * i; j++)
					cornerArray[i / size][j] = i + 1;
			}
		}

		if (k == 3)
		{
			for (i = rank; i < n; i+= size)
				for (j = 0; j <= i; j++)
				{
					cornerArray[i / size][j] = i - j;
					cornerArray[i / size][2 * i - j] = i - j;
				}
		}

		if (k == 4)
		{
			for (i = rank; i < n; i+= size)
				for (j = 0; j <= i; j++)
				{
					cornerArray[i / size][j] = (double)(1. / (i+j+1));
					cornerArray[i / size][2 * i - j] = (double)(1. / (i+j+1));
				}
		}
}
