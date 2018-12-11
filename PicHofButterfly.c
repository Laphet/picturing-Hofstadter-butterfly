#include <mkl.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "omp.h"

void save_data(FILE* f, int k, int l)
{
    double* e = (double *)malloc((2 * k + 1) * l * sizeof(double));
    double* off_diag = (double *)malloc(2 * k * sizeof(double));
    double delta_alpha = 1.0 / (double)(l + 1);
    int i, j;

    for (j = 0; j < l; ++j)
        for (i = -k; i <= k; ++i)
            e[j * (2 * k + 1) + i + k] = 2.0 * cos(2.0 * M_PI * (double)i * (double)(j + 1) * delta_alpha);

    for (j = 0; j < l; ++j)
    {
        for (i = 0; i < 2 * k; ++i) off_diag[i] = 1.0;
        LAPACKE_dsterf(2 * k + 1, &e[j * (2 * k + 1)], off_diag);
    }

    double alpha;
    for (j = 0; j < l; ++j)
    {
        alpha = delta_alpha * (double) (j + 1);
        for (i = 0; i < 2 * k + 1; ++i)
            fprintf(f, "%.8f,%.8f\n", e[j * (2 * k + 1) + i], alpha);
    }
    free(e);
    free(off_diag);
}

double get_consumed_time(int k, int l)
{
    clock_t start, end;
    double* e = (double *)malloc((2 * k + 1) * l * sizeof(double));
    double* off_diag = (double *)malloc(2 * k * l * sizeof(double));
    double delta_alpha = 1.0 / (double) (l + 1);
    int i, j;
    for (j = 0; j < l; ++j)
    {
        for (i = -k; i <= k; ++i)
            e[j * (2 * k + 1) + i + k] = 2.0 * cos(2.0 * M_PI * (double)i * (double)(j + 1) * delta_alpha);
        for (i = 0; i < 2 * k; ++i)
            off_diag[j * 2 * k + i] = 1.0;
    }
    start = clock();
    #pragma omp parallel for
    for (j = 0; j < l; ++j)
    {
        LAPACKE_dsterf(2 * k + 1, &e[j * (2 * k + 1)], &off_diag[j * 2 * k]);

    }
    end = clock();
    free(e);
    free(off_diag);
    return (double)(end - start) / CLOCKS_PER_SEC;
}

int main(int argc, char** argv)
{
    int k = 1024, l = 1024;
    if (argc == 3)
    {
        k = atoi(argv[1]);
        l = atoi(argv[2]);
    }
    //FILE *f = fopen("data.csv", "w+");
    //save_data(f, k, l);
    //fclose(f);
    printf("%.6f\n", get_consumed_time(k, l));
    int num_threads;
    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
    }
    printf("%d\n", num_threads);
}
