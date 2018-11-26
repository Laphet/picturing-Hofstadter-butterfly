/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2018, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Picturing Hofstadter's butterfly"
                     "The command line options are:\n"
                     "  -n <int>,  take center 2n+1 truncation of the infinite matrix, \n"
                     "the matrix needed to compute eigenvalues is 2n + 1.\n\n"
                     "  -m <int>,  slice [0, 1) into m pieces, where delta_alpha = 1 / m.\n\n";
#include <math.h>
#include <slepceps.h>
#define EIGENVALUE_LOWBOUND -4.0
#define EIGENVALUE_UPBOUND   4.0

int main(int argc, char **argv)
{
    Mat            A;           /* problem matrix */
    EPS            eps;         /* eigenproblem solver context */
    PetscReal      delta_alpha, alpha = 0.0;
    PetscInt       n = 30, m = 100, i, Istart, Iend, nev;
    PetscErrorCode ierr;
    ST             st;
    KSP            ksp;
    PC             pc;
    ierr = SlepcInitialize(&argc, &argv, (char*)0, help);
    if (ierr) return ierr;

    ierr = PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL);
    CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-m", &m, NULL);
    CHKERRQ(ierr);
    delta_alpha = 1.0 / m;

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Compute the operator matrix that defines the eigensystem, Ax=kx
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    ierr = MatCreate(PETSC_COMM_WORLD, &A);
    CHKERRQ(ierr);
    ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, 2 * n + 1, 2 * n + 1);
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);
    CHKERRQ(ierr);
    ierr = MatSetUp(A);
    CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(A, &Istart, &Iend);
    CHKERRQ(ierr);
    alpha += delta_alpha;
    for (i = Istart; i < Iend; i++)
    {
        if (i > 0)
        {
            ierr = MatSetValue(A, i, i - 1, 1.0, INSERT_VALUES);
            CHKERRQ(ierr);
        }
        if (i < n - 1)
        {
            ierr = MatSetValue(A, i, i + 1, 1.0, INSERT_VALUES);
            CHKERRQ(ierr);
        }
        ierr = MatSetValue(A, i, i, 2.0 * cos(2.0 * M_PI * alpha * (i - n)), INSERT_VALUES);
        CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                  Create the eigensolver and set various options
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
       Create eigensolver context
    */
    ierr = EPSCreate(PETSC_COMM_WORLD, &eps);
    CHKERRQ(ierr);

    /*
       Set operators. In this case, it is a standard eigenvalue problem
    */
    ierr = EPSSetOperators(eps, A, NULL);
    CHKERRQ(ierr);
    ierr = EPSSetProblemType(eps, EPS_HEP);
    CHKERRQ(ierr);
    /*
       Set solver parameters at runtime
    */
    ierr = EPSSetWhichEigenpairs(eps, EPS_ALL);
    CHKERRQ(ierr);
    ierr = EPSSetInterval(eps, EIGENVALUE_LOWBOUND, EIGENVALUE_UPBOUND);
    CHKERRQ(ierr);
    ierr = EPSGetST(eps, &st);
    CHKERRQ(ierr);
    ierr = STSetType(st, STSINVERT);
    CHKERRQ(ierr);
    ierr = STGetKSP(st, &ksp);
    CHKERRQ(ierr);
    ierr = KSPGetPC(ksp, &pc);
    CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPPREONLY);
    CHKERRQ(ierr);
    ierr = PCSetType(pc, PCLU);
    CHKERRQ(ierr);
    ierr = EPSSetFromOptions(eps);
    CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Solve the eigensystem
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = EPSSolve(eps);
    CHKERRQ(ierr);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Display solution and clean up
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = EPSGetDimensions(eps, &nev, NULL, NULL);
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, " Found %D eigenvalues converged", nev);
    CHKERRQ(ierr);
    /*
       Free work space
    */
    ierr = EPSDestroy(&eps);
    CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    CHKERRQ(ierr);
    ierr = SlepcFinalize();
    return ierr;
}
