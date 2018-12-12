#include "TriMatEigen_groups.h"
#include <sys/time.h>

int main(int argc, char** argv)
{
    double *alpha;
    double *beta;
    char option = 'a';
    int k = 4096, i, groups = 64;
    if (argc > 1)
    {
        option = *argv[1];
        if (argc > 2) groups = atoi(argv[2]);
        if (argc > 3) k = atoi(argv[3]);
    }

    alpha = (double *)malloc((2 * k + 1) * sizeof(double));
    beta = (double *)malloc(2 * k * sizeof(double));
    for (i = 0; i < 2 * k + 1; ++i)
    {
        alpha[i] = 2.0 * cos(2.0 * M_PI * (double)(i - k) * 0.34);
        if (i < 2 * k) beta[i] = 1.0;
    }
    Context ctx = {.order = 2 * k + 1, .low_bound = -4.0, .up_bound = 4.0, .tol = -1.0, .alpha = &alpha[0], .beta = &beta[0]};
    EigenArray eigens;
    struct timeval start, end;
    gettimeofday(&start, NULL);
    switch(option)
    {
    case 'a':
        eigens = solve_trimateigen(&ctx);
        break;
    case 'b':
        eigens = solve_trimateigen_omp(&ctx);
        break;
    case 'c':
        eigens = solve_trimateigen_groups(groups, &ctx);
        break;
    default :
        eigens = solve_trimateigen(&ctx);
    }
    gettimeofday(&end, NULL);
    double used_seconds = (double)((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    //FILE* f = fopen("test.dat", "w+");
    //write_EigenArray(f, &eigens);
    //fclose(f);
    //print_EigenArray(&eigens);
    //if (is_use_omp)
    //    printf("\nUsing divide-conquer algorithm openMP version to solve the eigenvalues of a n=%d tridiagonal matrix consumes time %fs\n", k, used_seconds);
    //else
    //    printf("\nUsing divide-conquer algorithm version to solve the eigenvalues of a n=%d tridiagonal matrix consumes time %fs\n", k, used_seconds);
    printf("%f\n", used_seconds);
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
