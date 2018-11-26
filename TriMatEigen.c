#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define EPSILON 1.0e-7
#define TOL 1.0e-12

typedef struct Context
{
    int order;
    double* alpha;
    double* beta;
    double low_bound;
    double up_bound;
} Context;

typedef struct DetEvl
{
    double df;
    double d2f;
    int count;
} DetEvl;

typedef struct EigenArray
{
    int total_size;
    int used_size;
    double *data;
    int *eigenvalue_index;
} EigenArray;

typedef struct OpenInterval
{
    double a;
    double b;
} OpenInterval;

EigenArray get_merged_eigenvalues(EigenArray a, EigenArray b)
{
    double *c = (double *)malloc((a.used_size + b.used_size) * sizeof(double));
    int i = 0, j = 0, k = 0;
    while (i < a.used_size && j < b.used_size)
    {
        if (a.data[i] <= b.data[j])
        {
            if (k >= 1 && fabs(a.data[i] - c[k - 1]) < TOL)
                c[k - 1] = 0.5 * (c[k - 1] + a.data[i++]);
            else
                c[k++] = a.data[i++];
        }
        else
        {
            if (k >= 1 && fabs(b.data[j] - c[k - 1]) < TOL)
                c[k - 1] = 0.5 * (c[k - 1] + b.data[j++]);
        }
    }
    while (i < a.used_size) c[k++] = a.data[i++];
    while (j < b.used_size) c[k++] = b.data[j++];
    return (EigenArray)
    {
        .total_size = a.used_size + b.used_size, .used_size = k, .data = c
    };
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
    int n = ctx->order;
    DetEvl det_evl_output = get_sturm_sequence(x, ctx);
    double df = det_evl_output.df, d2f = det_evl_output.d2f;
    if (is_positive)
        x += n / (df + sqrt(fabs((n - r) / r * ((n - 1) * df * df - n * d2f))));
    else
        x += n / (df - sqrt(fabs((n - r) / r * ((n - 1) * df * df - n * d2f))));
    return x;
}

int get_closest_upindex(int i, EigenArray lambdas_, Context *ctx)
{
    int j;
    for (j = 1; j < lambdas_.used_size; ++i)
    {
        if (get_low_eigenvalues_count(lambdas_.data[j], ctx) >= i) break;
    }
    return i;
}

double get_ith_eigenvalue(int i, double x, Context *ctx)
{
    double x_ = x, delta, delta_ = EPSILON, tol, *beta = ctx->beta;
    int kappa, kappa_, n = ctx->order, j, r = 1;
    double max_beta_jj1 = fabs(beta[0]) + fabs(beta[1]), temp;
    for (j = 1; j < n - 3; ++j)
        max_beta_jj1 = (fabs(beta[j]) + fabs(beta[j + 1]) > max_beta_jj1) ? fabs(beta[j] + fabs(beta[j + 1])) : max_beta_jj1;

    tol = (2.5 * max_beta_jj1 + fabs(x)) * EPSILON;
    while (1)
    {
        kappa_ = get_low_eigenvalues_count(x_, ctx);
        if (kappa_ == i) x = get_laguerre_iteration(x_, r, 1, ctx);
        else if (kappa_ == i + 1) x = get_laguerre_iteration(x_, r, 0, ctx);
        else
        {
            printf("Warning: the start point is badly chosen, may result in obtaining incorrect eigenvalue, try jumping");
            if (kappa_ > i + 1) x -= (ctx->up_bound - ctx->low_bound) / ctx->order * (kappa_ - i - 1);
            else x += (ctx->up_bound - ctx->low_bound) / ctx->order * (i - kappa_);
            continue;
        }
        printf("current x at %f\n", x);
        delta = fabs(x - x_);
        kappa = get_low_eigenvalues_count(x, ctx);
        if (delta < tol || delta * delta / delta_ < tol ) break;
        while(1)
        {
            if (abs(kappa - kappa_) > 1)
            {
                r = abs(kappa - kappa_);
                temp = (x + x_) / 2.0;
                x_ = x;
                x = temp;
                printf("--inner jump to x at %f\n", x);
            }
            else break;
        }
        delta_ = delta;
        x_ = x;
    }
    return x;
}

OpenInterval get_interval(int i, EigenArray lambdas_, Context *ctx)
{
    double a, b;
    int k;
    if (lambdas_.eigenvalue_index == NULL)
    {
        for (k = 0; k < lambdas_.used_size; ++k)
            lambdas_.eigenvalue_index[k] = get_low_eigenvalues_count(lambdas_.data[k], ctx);
    }
    if (i < lambdas_.eigenvalue_index[0])
    {
        a = ctx->low_bound - TOL;
        b = lambdas_.data[0];
        return (OpenInterval)
        {
            .a = a, .b = b
        };
    }
    if (lambdas_.eigenvalue_index [lambdas_.used_size - 1] <= i)
    {
        a = lambdas_.data[lambdas_.used_size - 1];
        b = ctx->up_bound + TOL;
        return (OpenInterval)
        {
            .a = a, .b = b
        };
    }
    for (k = 0; k < lambdas_.used_size - 1; ++k)
    {
        if (lambdas_.eigenvalue_index[k] <= i && i < lambdas_.eigenvalue_index[k + 1])
        {
            a = lambdas_.data[k];
            b = lambdas_.data[k + 1];
            return (OpenInterval)
            {
                .a = a, .b = b
            };
        }
    }
    printf{"Error: cannot find an interval contains i-th eigenvalue, use default instead"};
    return (OpenInterval)
    {
        .a = ctx->low_bound - TOL, .b = ctx->up_bound + TOL
    };

}

EigenArray get_eigenvalues(EigenArray lambdas_, Context *ctx)
{
    int index_low_bound = get_low_eigenvalues_count(ctx->up_bound + TOL, ctx), index_up_bound = get_low_eigenvalues_count(ctx->up_bound + TOL, ctx);
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
                if (det_evl_output.count < i)
                    interval.a = x;
                else
                    interval.b = x;
                if (det_evl_output.count == i && det_evl_output.df > 0.0)
                {
                    eigens[i - index_low_bound] = get_laguerre_iteration(x, 1, 0, ctx);
                    break;
                }
                else if ( det_evl_output.count == i + 1 && det_evl_output.df < 0.0)
                {
                    eigens[i - index_low_bound] = get_laguerre_iteration(x, 1, 0, ctx);
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
    return (EigenArray)
    {
        .total_size = eigenvalues_count, .used_size = i, .data = eigens, .eigenvalue_index = NULL
    };
}

int main (int argc, char **argv)
{
    int i = 5, j;
    double alpha[5] = {1.0, 1.1, 1.2, 1.3, 1.4};
    double beta[4] = {1.0, 1.0, 1.0, 1.0};
    double low_bound = -4.0;
    double up_bound = 4.0;
    Context ctx = {.order = i, .alpha = &alpha[0], .beta = &beta[0], .low_bound = low_bound, .up_bound = up_bound};
    EigenArray a = {.total_size = 5, .used_size = 5, .data = &alpha[0]};
    EigenArray b = {.total_size = 4, .used_size = 4, .data = &beta[0]};
    EigenArray c = get_merged_eigenvalues(a, b);
    for (j = 0; j < c.used_size; ++j) printf("\t%f", c.data[j]);
    printf("\n%0*d\n", 20, 0);
    printf("the first eigenvalue is %f\n", get_ith_eigenvalue(0, -4.0, 1, &ctx));
    printf("the last  eigenvalue is %f\n", get_ith_eigenvalue(4, 4.0, 1, &ctx));
}
