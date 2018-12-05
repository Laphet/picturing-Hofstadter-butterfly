#include "TriMatEigen_omp.h"
#include <time.h>

int main(int argc, char** argv)
{
    double *alpha;
    double *beta;
    int is_use_omp = 0;
    int k = 8, i, groups = 64;
    if (argc == 3)
    {
        k = atoi(argv[1]);
        groups = atoi(argv[2]);
        is_use_omp = 1;
    }
    if (argc == 2)
        k = atoi(argv[1]);

    alpha = (double *)malloc(k * sizeof(double));
    beta = (double *)malloc((k - 1) * sizeof(double));
    for (i = 0; i < k; ++i)
    {
        alpha[i] = 2.0 * cos(2.0 * M_PI * (double)i * 0.34);
        if (i < k - 1) beta[i] = 1.0;
    }
    Context ctx = {.order = k, .low_bound = -4.0, .up_bound = 4.0, .tol = -1.0, .alpha = &alpha[0], .beta = &beta[0]};
    EigenArray eigens;
    clock_t start = clock();
    //EigenArray eigens = solve_trimateigen_omp(groups, &ctx);
    if (is_use_omp) eigens = solve_trimateigen_omp(groups, &ctx);
    else eigens = solve_trimateigen(&ctx);
    clock_t end = clock();
    //FILE* f = fopen("test.dat", "w+");
    //write_EigenArray(f, &eigens);
    //fclose(f);
    //print_EigenArray(&eigens);
    if (is_use_omp)
        printf("\nUsing divide-conquer algorithm openMP version to solve the eigenvalues of a n=%d tridiagonal matrix consumes time %fs", k, (double)(end - start) / CLOCKS_PER_SEC);
    else
        printf("\nUsing divide-conquer algorithm version to solve the eigenvalues of a n=%d tridiagonal matrix consumes time %fs", k, (double)(end - start) / CLOCKS_PER_SEC);
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
