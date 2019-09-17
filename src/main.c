#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <assert.h>

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

void inplace_thomas(struct tridiag_t system) {
    struct tridiag_row_t *A = system.A;
    double *x = system.x;
    size_t n = system.n;

    // Forward reduction
    // Not vectorizable. Read after write data dependency
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

void print_system(struct tridiag_t system) {
    printf("| % 0.2f % 0.2f", system.A[0].a, system.A[0].b);
    printf("%*s|\n", (system.n-2)*6 + 1, " ");
    for (int i = 1; i < system.n-1; i++) {
        printf("|%*s", (i-1)*6+1, " ");
        printf("% 0.2f % 0.2f % 0.2f", system.A[i].a, system.A[i].b, system.A[i].c);
        printf("%*s|\n", (system.n-i-2)*6 + 1, " ");
    }
    printf("|%*s", (system.n-2)*6 + 1, " ");
    printf("% 0.2f % 0.2f |\n", system.A[system.n-1].b, system.A[system.n-1].c);
}

void parallel_thomas(struct tridiag_t system, double *x, int n_part) {
    int part_size = system.n/n_part;
    struct tridiag_row_t *A = system.A;

    struct tridiag_t interface_sys;
    interface_sys.A = malloc(sizeof(struct tridiag_row_t)*n_part*2);
    interface_sys.n = 2*n_part;
    interface_sys.x = malloc(interface_sys.n*sizeof(double));

    #pragma omp parallel shared(interface_sys) num_threads(n_part)
    {
        // Compute the interface equations on each processor
        #pragma omp for schedule(static)
        for (int partition = 0; partition < n_part; partition++) {
            int start = partition*part_size;

            struct tridiag_row_t *upper = &interface_sys.A[partition*2];
            struct tridiag_row_t *lower = &interface_sys.A[partition*2+1];

            upper->a = A[start + part_size - 2].a;
            upper->b = A[start + part_size - 2].b;
            upper->c = A[start + part_size - 2].c;
            upper->d = A[start + part_size - 2].d;
            for (int i = start + part_size - 3; i >= start; i--) {
                double beta = A[i].c/upper->b;
                upper->d = A[i].d - beta*upper->d;
                upper->c = - beta*upper->c;
                upper->b = A[i].b - beta*upper->a;
                upper->a = A[i].a;
            }

            lower->a = A[start + 1].a;
            lower->b = A[start + 1].b;
            lower->c = A[start + 1].c;
            lower->d = A[start + 1].d;
            for (int i = start + 2; i < start + part_size; i++) {
                double alpha = A[i].a/lower->b;
                lower->a = -alpha*lower->a;
                lower->b = A[i].b - alpha*lower->c;
                lower->c = A[i].c;
                lower->d = A[i].d - alpha*lower->d;
            }
        }

        // Collect and solve the interface equation on the master thread
        #pragma omp barrier
        if (omp_get_thread_num() == 0) {
            inplace_thomas(interface_sys);
        }
        #pragma omp barrier

        // Scatter the solved interface system among the processors
        #pragma omp for schedule(static)
        for (int partition = 0; partition < n_part; partition++) {
            system.A[partition*part_size].a = 0;
            system.A[partition*part_size].b = 1;
            system.A[partition*part_size].c = 0;
            system.A[partition*part_size].d = interface_sys.x[2*partition];

            system.A[(partition + 1)*part_size - 1].a = 0;
            system.A[(partition + 1)*part_size - 1].b = 1;
            system.A[(partition + 1)*part_size - 1].c = 0;
            system.A[(partition + 1)*part_size - 1].d = interface_sys.x[2*partition + 1];

            // Solve the independent tridiagonal systems
            struct tridiag_t sub_sys;
            sub_sys.n = part_size;
            sub_sys.A = &system.A[partition*part_size];
            sub_sys.x = &x[partition*part_size];
            inplace_thomas(sub_sys);
        }
    }
}


/**
 * Solve Ax_2 = x_1 in serial using Thomas algorithm where A is the tridiagonal matrix
 *     | (1-2a)   -a                        |
 *     |   -a   (1-2a)  -a                  |
 * A = |          -a   (1-2a)  -a           |
 *     |                  .....             |
 *     |                        -a   (1-2a) |
 * that arises from applying the backward Euler stencil.
 * 
 * @param a   The parameter a that defines the tridiagonal matrix A
 * @param x1  The known right-hand side
 * @param x2  The unknown to solve for
 * @param N   The dimension of the the system
 */
void implicit_stencil_serial(double a, double *x1, double *x2, int N) {
    x2[0] = 1-2*a; // Reuse memory to store the diagonal of A
    // Forward reduction. Not vectorizable -- read after write data dependency
    for (int i = 1; i < N; i++) {
        double w = -a/x2[i-1];
        x2[i] = (1-2*a) + w*a;
        x1[i] -= w*x1[i-1];
    }

    // Backward substitution. Again, not vectorizable
    x2[N-1] = x1[N-1]/x2[N-1];
    for (int i = N - 2; i >= 0; i--) {
        x2[i] = (x1[i] + a*x2[i+1])/x2[i];
    }
}

/**
 * Solve Ax_2 = x_1 in parallel where A is the tridiagonal matrix
 *     | (1-2a)   -a                        |
 *     |   -a   (1-2a)  -a                  |
 * A = |          -a   (1-2a)  -a           |
 *     |                  .....             |
 *     |                        -a   (1-2a) |
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
void implicit_stencil_parallel(double a, double *x1, double *x2, int N, int n) {
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
            upper->b = 1-2*a;
            upper->c = -a;
            upper->d = x1[start + k - 2];
            for (int i = start + k - 3; i >= start; i--) {
                double beta = -a/upper->b;
                upper->d = x1[i] - beta*upper->d;
                upper->c = - beta*upper->c;
                upper->b = (1-2*a) - beta*upper->a;
            }
            if (partition == 0) upper->a = 0;

            lower->a = -a;
            lower->b = 1-2*a;
            lower->c = -a;
            lower->d = x1[start + 1];
            for (int i = start + 2; i < start + k; i++) {
                double alpha = -a/lower->b;
                lower->a = -alpha*lower->a;
                lower->b = (1-2*a) - alpha*lower->c;
                lower->d = x1[i] - alpha*lower->d;
            }
            if (partition == n-1) lower->c = 0;
        }

        // Collect and solve the interface system on the master thread
        #pragma omp barrier
        if (omp_get_thread_num() == 0) {
            inplace_thomas(interface_sys);
        }
        #pragma omp barrier

        // Distribution the solution of the interface system back among the processors
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

            // Solve the n independent tridiagonal k x k systems of the form
            // |    1                                  | |  x2_1  |   |  x1_1  |
            // |   -a   (1-2a)  -a                     | |  x2_2  |   |  x1_2  |
            // |          -a   (1-2a)  -a              | |  x2_3  | = |  x1_3  |
            // |                  .....                | |  ...   |   |  ...   |
            // |                       -a   (1-2a)  -a | |x2_(k-1)|   |x1_(k-1)|
            // |                                     1 | |  x2_k  |   |  x1_k  |
            x1[start+1] += a*x1[start];
            x1[start+k-2] += a*x1[start+k-1];
            x2[start] = x1[start];
            x2[start+k-1] = x1[start+k-1];

            double *x1_part = &x1[start + 1];
            double *x2_part = &x2[start + 1];
            implicit_stencil_serial(a, x1_part, x2_part, k-2);
        }
    }
    free(interface_sys.A);
    free(interface_sys.x);
}


int main() {
    int n = 1<<5;
    double a = 0.2;
    double *x1 = malloc(sizeof(double)*n);
    double *x2 = malloc(sizeof(double)*n);

    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        x1[i] = (i < n/2) ? ((double) (i+1))/(n+1) : 1.0 - ((double) (i+1))/(n+1);
    }

    printf("Running parallel solver\n");
    double start = omp_get_wtime();
    implicit_stencil_parallel(0.2, x1, x2, n, 5);
    double total = omp_get_wtime() - start;

    printf("%0.4f: Total time = %f\n", x2[n/2], total);

    for (int i = 0; i < n; i++) {
        printf("%0.5f\n", x2[i]);
    }
}
