#include "TriMatEigen_omp.h"
#include <time.h>

int main()
{
    double *alpha;
    double *beta;
    int k = 1600, i, groups = 3;
    alpha = (double *)malloc(k * sizeof(double));
    beta = (double *)malloc((k - 1) * sizeof(double));
    for (i = 0; i < k; ++i)
    {
        alpha[i] = 2.0 * cos(2.0 * M_PI * (double)i * 0.34);
        if (i < k - 1) beta[i] = 1.0;
    }
    Context ctx = {.order = k, .low_bound = -4.0, .up_bound = 4.0, .tol = -1.0, .alpha = &alpha[0], .beta = &beta[0]};
    clock_t start = clock();
    EigenArray eigens = solve_trimateigen_omp(groups, &ctx);
    clock_t end = clock();
    //FILE* f = fopen("test.dat", "w+");
    //write_EigenArray(f, &eigens);
    //fclose(f);
    //print_EigenArray(&eigens);
    printf("Using divide-conquer algorithm to solve n=%d tridiagonal matrix consumes time %fs", k, (double)(end - start) / CLOCKS_PER_SEC);
    free_EigenArray(&eigens);
    free(alpha);
    free(beta);
    /* printf("omp_get_num_threads, return %d\n", omp_get_num_threads());
      #pragma omp parallel
      {
          printf("omp_get_num_threads, return %d\n", omp_get_num_threads());
      }

      double *share_var;
      double private_var;
      #pragma omp parallel
      {
          int k = omp_get_thread_num();
          private_var = (double)k + M_PI;
          share_var = &private_var;
      }
      printf("share is %f\n", *share_var);
      */
}
