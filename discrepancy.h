#pragma once

double disc(double **cornerArray, double *b, double *x, double *buf, int n, int k, int size, int rank);

double solution_err(double* x, double *buf, int n, int size, int rank);
