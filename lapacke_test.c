#include "mkl.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>

int main(int argc, char** argv)
{
    int n = 8192, i;
    if (argc == 2) n = atoi(argv[1]);
    double *d = (double *)malloc(n * sizeof(double));
    double *e = (double *)malloc((n - 1) * sizeof(double));
    for (i = 0; i < n; ++i)
    {
        d[i] = 2.0 * cos(2.0 * M_PI * (double)i * 0.34);
        if (i < n - 1) e[i] = 1.0;
    }
    struct timeval start, end;
    gettimeofday(&start, NULL);
    LAPACKE_dsterf (n, d, e);
    gettimeofday(&end, NULL);
    double used_seconds = (double)((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    printf("\nUsing MKL lapacke to solve the eigenvalues of a n=%d tridiagonal matrix consumes time %fs", n, used_seconds);
}
