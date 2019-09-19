#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "implicit_step.h"

struct tridiag_row_t {
    double a;
    double b;
    double c;
    double d;
};

struct tridiag_t {
    struct tridiag_row_t *A;
    double *x;
    size_t n;
};

/**
 * Solve a general tridiagonal system using Thomas algorithm
 */
void inplace_thomas(struct tridiag_t system) {
    struct tridiag_row_t *A = system.A;
    double *x = system.x;
    size_t n = system.n;

    // Forward reduction. Not vectorizable -- read after write data dependency
    for (int i = 1; i < n; i++) {
        double w = A[i].a/A[i-1].b;
        A[i].b -= w*A[i-1].c;
        A[i].d -= w*A[i-1].d;
    }

    // Backward substitution. Again, not vectorizable
    x[n-1] = A[n-1].d/A[n-1].b;
    for (int i = n - 2; i >= 0; i--) {
        x[i] = (A[i].d - A[i].c*x[i+1])/A[i].b;
    }
}

/**
 * Solve Ax_2 = x_1 in serial using Thomas algorithm where A is the tridiagonal matrix
 *     |    1                                    |
 *     |   -a   (1+2a)  -a                       |
 * A = |          -a   (1+2a)  -a                |
 *     |                  .....                  |
 *     |                       -a   (1+2a)   -a  |
 *     |                                      1  |
 * that arises from applying the backward Euler stencil.
 * 
 * Knowing the form of the matrix a priori saves the memory needed to store
 * the matrix.
 * 
 * @param a   The parameter a that defines the tridiagonal matrix A
 * @param x1  The known right-hand side
 * @param x2  The unknown to solve for
 * @param N   The dimension of the the system
 */
void implicit_step_serial(double a, double *x1, double *x2, int N) {
    x1[1] += a*x1[0];
    x1[N-2] += a*x1[N-1];
    x2[0] = x1[0];
    x2[N-1] = x1[N-1];

    x2[1] = 1+2*a; // Reuse memory to store the diagonal of A
    // Forward reduction. Not vectorizable -- read after write data dependency
    for (int i = 2; i < N-1; i++) {
        double w = -a/x2[i-1];
        x2[i] = (1+2*a) + w*a;
        x1[i] -= w*x1[i-1];
    }

    // Backward substitution. Again, not vectorizable
    x2[N-2] = x1[N-2]/x2[N-2];
    for (int i = N - 3; i > 0; i--) {
        x2[i] = (x1[i] + a*x2[i+1])/x2[i];
    }
}

/**
 * Solve Ax_2 = x_1 in parallel where A is the tridiagonal matrix
 *     |    1                                    |
 *     |   -a   (1+2a)  -a                       |
 * A = |          -a   (1+2a)  -a                |
 *     |                  .....                  |
 *     |                       -a   (1+2a)   -a  |
 *     |                                      1  |
 * that arises from applying the backward Euler stencil.
 * 
 * The algorithm is based off of an adaptation of the Thomas algorithm presented
 * in "A Memory Efficient Parallel Tridiagonal Solver" by T. Austin, et al.
 * 
 * It works by preprocessing the matrix A into a block-diagonal matrix of
 * tridiagonal submatrices. This reduces the system into k indepedent tri-
 * diagonal systems that can be solved in parallel using the traditional
 * Thomas algorithm.
 * 
 * @param a   The parameter a that defines the tridiagonal matrix A
 * @param x1  The known right-hand side
 * @param x2  The unknown to solve for
 * @param N   The dimension of the the system
 * @param n   The number of processors to parallelize the computation over
 */
void implicit_step_parallel(double a, double *x1, double *x2, int N, int n) {
    struct tridiag_t interface_sys;
    interface_sys.n = 2*n;
    interface_sys.A = malloc(interface_sys.n * sizeof(struct tridiag_row_t));
    interface_sys.x = malloc(interface_sys.n * sizeof(double));

    #pragma omp parallel shared(interface_sys, x1, x2)
    {
        // Compute the interface equations on each processor
        #pragma omp for schedule(static)
        for (int partition = 0; partition < n; partition++) {
            // Divide the matrix as evenly as possible among the processors
            int k = N/n;
            int remainder = N - k*n;
            int start = k*partition;
            if (partition < remainder) {
                start += partition;
                k++;
            } else {
                start += remainder;
            }

            struct tridiag_row_t *upper = &interface_sys.A[partition*2];
            struct tridiag_row_t *lower = &interface_sys.A[partition*2+1];

            upper->a = -a;
            upper->b = 1+2*a;
            upper->c = -a;
            upper->d = x1[start + k - 2];
            for (int i = start + k - 3; i >= start; i--) {
                double beta = -a/upper->b;
                upper->d = x1[i] - beta*upper->d;
                upper->c = - beta*upper->c;
                upper->b = (1+2*a) - beta*upper->a;
            }
            if (partition == 0) {
                upper->a = 0;
                upper->b = 1;
                upper->c = 0;
                upper->d = x1[0];
            }

            lower->a = -a;
            lower->b = 1+2*a;
            lower->c = -a;
            lower->d = x1[start + 1];
            for (int i = start + 2; i < start + k; i++) {
                double alpha = -a/lower->b;
                lower->a = -alpha*lower->a;
                lower->b = (1+2*a) - alpha*lower->c;
                lower->d = x1[i] - alpha*lower->d;
            }
            if (partition == n-1) {
                lower->a = 0;
                lower->b = 1;
                lower->c = 0;
                lower->d = x1[N-1];
            }
        }

        // Collect and solve the interface system on the master thread
        #pragma omp barrier
        if (omp_get_thread_num() == 0) {
            inplace_thomas(interface_sys);
        }
        #pragma omp barrier

        // Distribute the solution of the interface system back among the processors
        #pragma omp for schedule(static)
        for (int partition = 0; partition < n; partition++) {
            // Divide the matrix as evenly as possible among the processors
            int k = N/n;
            int remainder = N - k*n;
            int start = k*partition;
            if (partition < remainder) {
                start += partition;
                k++;
            } else {
                start += remainder;
            }

            x1[start]        = interface_sys.x[2*partition];
            x1[start + k -1] = interface_sys.x[2*partition + 1];

            // Solve the n independent tridiagonal k x k systems which now have
            // the same form as the original system
            double *x1_part = &x1[start];
            double *x2_part = &x2[start];
            implicit_step_serial(a, x1_part, x2_part, k);
        }
    }
    free(interface_sys.A);
    free(interface_sys.x);
}