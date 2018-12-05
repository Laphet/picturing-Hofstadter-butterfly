#include "TriMatEigen_omp.h"
#include <time.h>

EigenArray solve_trimateigen_omp(int groups, Context *ctx)
{
    int k = ctx->order / groups, r = ctx->order % groups, i, j;
    if (k > 0)
    {
        Context *ctx_list = (Context *)malloc(groups * sizeof(Context));
        EigenArray *eigens_list = (EigenArray *)malloc(groups * sizeof(EigenArray));
        #pragma omp parallel for
        for (i = 0; i < groups; ++i)
        {
            //int ID = omp_get_thread_num();
            Context ctx_local;
            if (i == 0)
            {
                ctx_local = (Context)
                {
                    .order = k + r, .low_bound = ctx->low_bound, .up_bound = ctx->up_bound, .tol = -1.0,
                    .alpha = &(ctx->alpha[0]),
                    .beta = &(ctx->beta[0])
                };
            }
            else
            {
                ctx_local = (Context)
                {
                    .order = k, .low_bound = ctx->low_bound, .up_bound = ctx->up_bound, .tol = -1.0,
                    .alpha = &(ctx->alpha[i * k + r - 1]),
                    .beta = (i * k + r - 1 < ctx->order - 1) ? &(ctx->beta[i * k + r - 1]) : NULL
                };
            }
            //clock_t start = clock();
            eigens_list[i] = solve_trimateigen(&ctx_local);
            //clock_t end = clock();
            //printf("thread %d at leaf %d uses time %f\n", ID, i, (double)(end - start) / CLOCKS_PER_SEC);
            ctx_list[i] = ctx_local;
        }

        for (i = 1; i < groups; i = 2 * i)
            #pragma omp parallel for
            for (j = i; j < groups; j = j + (2 * i))
            {
                eigens_list[j - i] = get_merged_eigenvalues(&eigens_list[j - i], &eigens_list[j]);
                ctx_list[j - i].order = ctx_list[j - i].order + ctx_list[j].order;
                ctx_list[j - i].tol = -1.0;
                get_eigenvalues(&eigens_list[j - i], &ctx_list[j - i]);
            }

        #pragma omp parallel for
        for (i = 1; i < groups; ++i)
            free_EigenArray(&eigens_list[i]);
        EigenArray output = eigens_list[0];
        free(eigens_list);
        free(ctx_list);
        return output;
    }
    else return solve_trimateigen(ctx);
}
