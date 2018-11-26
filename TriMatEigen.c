#include "TriMatEigen.h"

typedef struct OpenInterval
{
    double a;
    double b;
} OpenInterval;


typedef struct DetEvl
{
    double df;
    double d2f;
    int count;
} DetEvl;


void free_EigenArray(EigenArray *e)
{
    free(e->data);
    free(e->eigenvalue_index);
}

void print_EigenArray(EigenArray *e)
{
    int i;
    printf("\n");
    for (i = 0; i < e->used_size; ++i) printf("%f\t", e->data[i]);
    printf("\n");
}

EigenArray get_merged_eigenvalues(EigenArray *a, EigenArray *b)
{
    double *c = (double *)malloc((a->used_size + b->used_size) * sizeof(double));
    int i = 0, j = 0, k = 0;
    while (i < a->used_size && j < b->used_size)
    {
        if (a->data[i] <= b->data[j])
        {
            if (k >= 1 && fabs(a->data[i] - c[k - 1]) < TOL)
                c[k - 1] = 0.5 * (c[k - 1] + a->data[i++]);
            else
                c[k++] = a->data[i++];
        }
        else
        {
            if (k >= 1 && fabs(b->data[j] - c[k - 1]) < TOL)
                c[k - 1] = 0.5 * (c[k - 1] + b->data[j++]);
            else
                c[k++] = b->data[j++];
        }
    }
    while (i < a->used_size) c[k++] = a->data[i++];
    while (j < b->used_size) c[k++] = b->data[j++];
    EigenArray output = {.total_size = a->used_size + b->used_size, .used_size = k, .data = c, .eigenvalue_index = NULL};
    free_EigenArray(a);
    free_EigenArray(b);
    return output;
}

int get_low_eigenvalues_count(double lambda, Context *ctx)
{
    double xi;
    int count = 0, i, n = ctx->order;
    double *alpha = ctx->alpha, *beta = ctx->beta, temp;
    xi = alpha[0] - lambda;
    if (xi == 0.0) xi = beta[0] * beta[0] * EPSILON * EPSILON;
    if (xi < 0.0) count += 1;
    for (i = 1; i < n; ++i)
    {
        temp = beta[i - 1] * beta[i - 1] / xi;
        xi = alpha[i] - lambda - temp;
        if (xi == 0.0) xi = temp * EPSILON * EPSILON;
        else if(xi < 0.0) count += 1;
    }
    return count;

}

DetEvl get_sturm_sequence(double lambda, Context *ctx)
{
    double xi, eta, zeta, eta_ = 0.0, zeta_ = 0.0, eta_temp = 0.0, zeta_temp = 0.0, temp = 0.0;
    int count = 0, i, n = ctx->order;
    double *alpha = ctx->alpha, *beta = ctx->beta;
    xi = alpha[0] - lambda;
    if (xi == 0.0) xi = beta[0] * beta[0] * EPSILON * EPSILON;
    eta = 1.0 / xi;
    zeta = 0.0;
    if (xi < 0.0) count++;
    for (i = 1; i < n; ++i)
    {
        temp = beta[i - 1] * beta[i - 1] / xi;
        xi = alpha[i] - lambda - temp;
        if (xi == 0.0)
            xi = temp * EPSILON * EPSILON;
        else if (xi < 0.0)
            count += 1;

        eta_temp = 1.0 / xi * ((alpha[i] - lambda) * eta + 1.0 - temp * eta_);
        zeta_temp = 1.0 / xi * ((alpha[i] - lambda) * zeta + 2.0 * eta - temp * zeta_);
        eta_ = eta;
        eta = eta_temp;
        zeta_ = zeta;
        zeta = zeta_temp;
    }
    return (DetEvl)
    {
        .df = eta, .d2f = zeta, .count = count
    };
}

double get_laguerre_iteration(double x, int r, int is_positive, Context *ctx)
{
    double n = (double)ctx->order;
    DetEvl det_evl_output = get_sturm_sequence(x, ctx);
    double df = det_evl_output.df, d2f = det_evl_output.d2f;
    if (is_positive)
        x += n / (df + sqrt(fabs((n - r) / r * ((n - 1) * df * df - n * d2f))));
    else
        x += n / (df - sqrt(fabs((n - r) / r * ((n - 1) * df * df - n * d2f))));
    return x;
}

double get_ith_eigenvalue(int i, double x, Context *ctx)
{
    double x_ = x, delta, delta_ = EPSILON, *beta = ctx->beta;
    int kappa, kappa_ = get_low_eigenvalues_count(x_, ctx), n = ctx->order, j, r = 1;
    double temp;
    if (ctx->tol <= 0.0)
    {
        double  max_beta_jj1 = fabs(beta[0]) + fabs(beta[1]);
        for (j = 1; j < n - 3; ++j)
            max_beta_jj1 = (fabs(beta[j]) + fabs(beta[j + 1]) > max_beta_jj1) ? fabs(beta[j] + fabs(beta[j + 1])) : max_beta_jj1;
        ctx->tol = (2.5 * max_beta_jj1 + fabs(x)) * EPSILON;
    }

    while (1)
    {
        if (kappa_ == i) x = get_laguerre_iteration(x_, r, 1, ctx);
        else if (kappa_ == i + 1) x = get_laguerre_iteration(x_, r, 0, ctx);
        else
        {
            printf("Warning: the start point is badly chosen, may result in obtaining incorrect eigenvalue, try jumping");
            if (kappa_ > i + 1) x -= (ctx->up_bound - ctx->low_bound) / ctx->order * (kappa_ - i - 1);
            else x += (ctx->up_bound - ctx->low_bound) / ctx->order * (i - kappa_);
            continue;
        }
        // printf("current x at %f\n", x);
        delta = fabs(x - x_);
        kappa = get_low_eigenvalues_count(x, ctx);
        if (delta < ctx->tol || delta * delta / delta_ < ctx->tol ) break;
        while(1)
        {
            if (abs(kappa - kappa_) > 1)
            {
                r = abs(kappa - kappa_);
                temp = 0.5 * (x + x_);
                x_ = x;
                x = temp;
                // printf("--inner jump to x at %f\n", x);
            }
            else break;
        }
        kappa_ = kappa;
        delta_ = delta;
        x_ = x;
    }
    return x;
}

OpenInterval get_interval(int i, EigenArray *lambdas_, Context *ctx)
{
    double a, b;
    int k;
    if (lambdas_->eigenvalue_index == NULL)
    {
        lambdas_->eigenvalue_index = (int *) malloc(lambdas_->used_size * sizeof(double));
        for (k = 0; k < lambdas_->used_size; ++k)
            lambdas_->eigenvalue_index[k] = get_low_eigenvalues_count(lambdas_->data[k], ctx);
    }
    if (i < lambdas_->eigenvalue_index[0])
    {
        a = ctx->low_bound - TOL;
        b = lambdas_->data[0];
        return (OpenInterval)
        {
            .a = a, .b = b
        };
    }
    else if (lambdas_->eigenvalue_index [lambdas_->used_size - 1] <= i)
    {
        a = lambdas_->data[lambdas_->used_size - 1];
        b = ctx->up_bound + TOL;
        return (OpenInterval)
        {
            .a = a, .b = b
        };
    }
    else
    {
        for (k = 0; k < lambdas_->used_size - 1; ++k)
        {
            if (lambdas_->eigenvalue_index[k] <= i && i < lambdas_->eigenvalue_index[k + 1])
            {
                a = lambdas_->data[k];
                b = lambdas_->data[k + 1];
                return (OpenInterval)
                {
                    .a = a, .b = b
                };
            }
        }
    }
    printf("Error: cannot find an interval contains i-th eigenvalue, use default instead");
    return (OpenInterval)
    {
        .a = ctx->low_bound - TOL, .b = ctx->up_bound + TOL
    };

}

EigenArray get_eigenvalues(EigenArray *lambdas_, Context *ctx)
{
    int index_low_bound = get_low_eigenvalues_count(ctx->low_bound - TOL, ctx), index_up_bound = get_low_eigenvalues_count(ctx->up_bound + TOL, ctx);
    int eigenvalues_count = index_up_bound - index_low_bound;
    int i;
    double *eigens = (double *)malloc(eigenvalues_count * sizeof(double));
    double x;
    DetEvl det_evl_output;
    OpenInterval interval;
    for (i = index_low_bound; i < index_up_bound; ++i)
    {
        interval = get_interval(i, lambdas_, ctx);
        while (1)
        {

            if  (fabs(interval.b - interval.a) > TOL)
            {

                x = 0.5 * (interval.a + interval.b);
                det_evl_output = get_sturm_sequence(x, ctx);
                if (det_evl_output.count <= i)
                    interval.a = x;
                else
                    interval.b = x;
                if ((det_evl_output.count == i && det_evl_output.df > 0.0) || ( det_evl_output.count == i + 1 && det_evl_output.df < 0.0))
                {
                    eigens[i - index_low_bound] = get_ith_eigenvalue(i, x, ctx);
                    break;
                }
            }
            else
            {
                eigens[i - index_low_bound] = 0.5 * (interval.a + interval.b);
                break;
            }
        }
    }
    free_EigenArray(lambdas_);
    return (EigenArray)
    {
        .total_size = eigenvalues_count, .used_size = i, .data = eigens, .eigenvalue_index = NULL
    };
}

EigenArray solve_trimateigen(Context *ctx)
{
    double *eigens;
    if (ctx->order == 1)
    {
        eigens = (double *)malloc(sizeof(double));
        eigens[0] = ctx->alpha[0];
        return (EigenArray)
        {
            .total_size = 1, .used_size = 1, .data = eigens, .eigenvalue_index = NULL
        };
    }
    if (ctx->order == 2)
    {
        double a0 = ctx->alpha[0], a1 = ctx->alpha[1], b0 = ctx->beta[0];
        eigens = (double *)malloc(2 * sizeof(double));
        eigens[0] = 0.5 * (a0 + a1) - 0.5 * sqrt(4.0 * b0 * b0 + (a0 - a1) * (a0 - a1));
        eigens[1] = 0.5 * (a0 + a1) + 0.5 * sqrt(4.0 * b0 * b0 + (a0 - a1) * (a0 - a1));
        return (EigenArray)
        {
            .total_size = 2, .used_size = 2, .data = eigens, .eigenvalue_index = NULL
        };
    }
    if (ctx->order > 2)
    {
        int n = ctx->order;
        Context ctx_child0 = {.order = n / 2, .low_bound = ctx->low_bound, .up_bound = ctx->up_bound, .tol = -1.0, .alpha = &ctx->alpha[0], .beta = n / 2 > 1 ? &ctx->beta[0] : NULL};
        Context ctx_child1 = {.order = n - n / 2, .low_bound = ctx->low_bound, .up_bound = ctx->up_bound, .tol = -1.0, .alpha = &ctx->alpha[n / 2], .beta = &ctx->beta[n / 2 - 1]};
        EigenArray lambda0 = solve_trimateigen(&ctx_child0);
        EigenArray lambda1 = solve_trimateigen(&ctx_child1);
        EigenArray lambdas_ = get_merged_eigenvalues(&lambda0, &lambda1);
        EigenArray lambdas = get_eigenvalues(&lambdas_, ctx);
        return lambdas;
    }
    printf("Error: ~");
    return (EigenArray)
    {
        .total_size = 0, .used_size = 1, .data = NULL, .eigenvalue_index = NULL
    };
}
