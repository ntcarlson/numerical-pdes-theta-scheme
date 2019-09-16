#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <time.h>

struct tridiag_row_t {
    double a;
    double b;
    double c;
    double d;
};

struct tridiag_t {
    struct tridiag_row_t *A;
    size_t n;
};

void inplace_thomas(struct tridiag_t system, double *x) {
    struct tridiag_row_t *A = system.A;
    size_t n = system.n;

    // Not vectorizable. Read after write data dependency
    for (int i = 1; i < n; i++) {
        double w = A[i].a/A[i-1].b;
        A[i].b -= w*A[i-1].c;
        A[i].d -= w*A[i-1].d;
    }

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

void phase1(struct tridiag_t system, struct tridiag_t interface_sys, int n_part) {
    int part_size = system.n/n_part;
    int n = system.n;
    struct tridiag_row_t *A = system.A;

    for (int partition = 0; partition < n_part; partition++) {
        int start = partition*part_size;
        struct tridiag_row_t *lower = &interface_sys.A[partition*2+1];
        struct tridiag_row_t *upper = &interface_sys.A[partition*2];
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
    }
}

void parallel_thomas(struct tridiag_t system, double *x, int n_part) {
    int part_size = system.n/n_part;
    struct tridiag_row_t *A = system.A;

    struct tridiag_t interface_sys;
    interface_sys.A = malloc(sizeof(struct tridiag_row_t)*n_part*2);
    interface_sys.n = 2*n_part;

    #pragma omp parallel
    {
        // Compute the interface equations on each processor
        #pragma omp for
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
         inplace_thomas(interface_sys, x);
        }
        #pragma omp barrier

        // Scatter the solved interface system among the processors
        #pragma omp for
        for (int partition = 0; partition < n_part; partition++) {
            system.A[partition*part_size].a = 0;
            system.A[partition*part_size].b = 1;
            system.A[partition*part_size].c = 0;
            system.A[partition*part_size].d = x[2*partition];

            system.A[(partition + 1)*part_size - 1].a = 0;
            system.A[(partition + 1)*part_size - 1].b = 1;
            system.A[(partition + 1)*part_size - 1].c = 0;
            system.A[(partition + 1)*part_size - 1].d = x[2*partition + 1];
        }

        // Solve the independent tridiagonal systems
        #pragma omp for
        for (int partition = 0; partition < n_part; partition++) {
            struct tridiag_t sub_sys;
            sub_sys.n = part_size;
            sub_sys.A = &system.A[partition*part_size];
            double *sub_x = &x[partition*part_size];
            inplace_thomas(sub_sys, sub_x);
        }
    }
}

int main() {
    int n = 1<<27;
    //int n = 16;
    double a = 0.2;
    double *x = malloc(sizeof(double)*n);

    struct tridiag_t system;
    system.A = malloc(sizeof(struct tridiag_row_t)*n);
    system.n = n;

    for (int i = 0; i < n; i++) {
        system.A[i].a = -a;
        system.A[i].b = 1-2*a;
        system.A[i].c = -a;
        system.A[i].d = (i < n/2) ? ((double) (i+1))/(n+1) : 1.0 - ((double) (i+1))/(n+1);
    }
    system.A[0].a = 0;
    system.A[n-1].c = 0;


    printf("Running parallel solver\n");
    double start = omp_get_wtime();
    parallel_thomas(system, x, 8);
    double total = omp_get_wtime() - start;
    //inplace_thomas(system, x);

    printf("%0.4f: Total time = %f\n", x[n-1], total);

    // for (int i = 0; i < n; i++) {
    //     printf("%0.5f\n", x[i]);
    // }
}
