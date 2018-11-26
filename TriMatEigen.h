#ifndef _TRI_MAT_EIGEN_H
#define _TRI_MAT_EIGEN_H
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define EPSILON 1.0e-7
#define TOL 1.0e-12

typedef struct Context
{
    int order;
    double low_bound;
    double up_bound;
    double tol; //default -1.0
    double* alpha;
    double* beta;
} Context;

typedef struct EigenArray
{
    int total_size;
    int used_size;
    double *data;
    int *eigenvalue_index;
} EigenArray;

void free_EigenArray(EigenArray*);

void print_EigenArray(EigenArray*);

EigenArray get_merged_eigenvalues(EigenArray*, EigenArray*);

EigenArray get_eigenvalues(EigenArray*, Context*);

EigenArray solve_trimateigen(Context*);

#endif