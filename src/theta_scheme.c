#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "theta_scheme.h"
#include "implicit_step.h"

/**
 * Apply the forward Euler step which involves computing U_2 = AU_1 where
 *     | (1-2a)    a                        |
 *     |    a   (1-2a)   a                  |
 * A = |           a   (1-2a)    a          |
 *     |                  .....             |
 *     |                         a   (1-2a) |
 */
void explicit_step(double a, double *U1, double *U2, int N) {
    // Boundary conditions
    U2[0]   = 0;
    U2[N-1] = 0;

    #pragma omp parallel for schedule(static)
    for (int i = 1; i < N-1; i++) {
        U2[i] = a*U1[i-1] + (1-2*a)*U1[i] + a*U1[i+1];
    }
}

double max_error(struct bvp_t bvp, double *U, double (*u)(double, double), double t) {
    double error = 0;
    #pragma omp parallel shared(error)
    {
        #pragma omp for reduction(max : error)
        for (int i = 0; i <= bvp.J; i++) {
            double x = ((double) i)/bvp.J;
            error = fabs(U[i] - u(x,t));
        }
    }
    return error;
}

double theta_scheme(struct bvp_t bvp, double (*u)(double, double)) {
    int N = bvp.J+1;
    double *U1 = malloc(sizeof(double)*N);
    double *U2 = malloc(sizeof(double)*N);

    double dx = 1.0/bvp.J;
    double mu = bvp.dt/(dx*dx);

    // Apply initial conditions
    U1[0]   = 0;
    U1[N-1] = 0;
    #pragma omp parallel for schedule(static)
    for (int i = 1; i < N-1; i++) {
        U1[i] = bvp.IC(i*dx);
    }
    
    double t = 0;
    double start_time = omp_get_wtime();
    double max_err = 0;
    if (bvp.theta == 0) {
        printf("Using forward Euler scheme\n");
        // Forward Euler method
        while (t < bvp.stop_time) {
            explicit_step(mu, U1, U2, N);
            double * tmp = U2;
            U2 = U1;
            U1 = tmp;
            if (t > 0.1) {
                double err = max_error(bvp, U1, u, t);
                if (err > max_err) {
                    max_err = err;
                }
            }

            t +=  bvp.dt;
        }
    } else if (bvp.theta == 1) {
        printf("Using backward Euler scheme\n");
        // Backward Euler method
        while (t < bvp.stop_time) {
            implicit_step_parallel(mu, U1+1, U2+1, N-2, 8);
            double * tmp = U2;
            U2 = U1;
            U1 = tmp;
            if (t > 0.1) {
                double err = max_error(bvp, U1, u, t);
                if (err > max_err) {
                    max_err = err;
                }
            }

            t +=  bvp.dt;
        }
    } else {
        // Theta scheme
        printf("Using theta scheme (theta = %0.2f)\n", bvp.theta);
        while (t < bvp.stop_time) {
            explicit_step((1-bvp.theta)*mu, U1, U2, N);
            implicit_step_parallel(bvp.theta*mu, U2+1, U1+1, N-2, 8);
            if (t > 0.1) {
                double err = max_error(bvp, U1, u, t);
                if (err > max_err) {
                    max_err = err;
                }
            }

            t +=  bvp.dt;
        }
    }
    double total_time = omp_get_wtime() - start_time;
    printf("Computed %d timesteps in %f seconds\n", (int) (bvp.stop_time/bvp.dt), total_time);

    free(U1);
    free(U2);
    return max_err;
}