#include <stdlib.h>
#include <stdio.h>
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

    #pragma omp for
    for (int i = 1; i < N-1; i++) {
        U2[i] = a*U1[i-1] + (1-2*a)*U1[i] + a*U1[i+1];
    }
}

void theta_scheme(struct bvp_t bvp) {
    int N = bvp.J+1;
    double *U1 = malloc(sizeof(double)*N);
    double *U2 = malloc(sizeof(double)*N);

    double dx = 1.0/bvp.J;
    double mu = bvp.dt/(dx*dx);
    // Apply initial conditions
    U1[0]   = 0;
    U1[N-1] = 0;
    for (int i = 1; i < N-1; i++) {
        U1[i] = bvp.IC(i*dx);
    }
    
    double t = 0;
    if (bvp.theta == 0) {
        printf("Using forward Euler scheme\n");
        // Forward Euler method
        while (t < bvp.stop_time) {
            explicit_step(mu, U1, U2, N);
            double * tmp = U2;
            U2 = U1;
            U1 = tmp;
            t +=  bvp.dt;
        }
    } else if (bvp.theta == 1) {
        printf("Using backward Euler scheme\n");
        // Backward Euler method
        while (t < bvp.stop_time) {
            implicit_step_parallel(mu, U1+1, U2+1, N-2, 4);
            double * tmp = U2;
            U2 = U1;
            U1 = tmp;
            t +=  bvp.dt;
        }
    } else {
        // Theta scheme
        printf("Using theta scheme\n");
        while (t < bvp.stop_time) {
            explicit_step((1-bvp.theta)*mu, U1, U2, N);
            implicit_step_parallel(bvp.theta*mu, U2+1, U1+1, N-2, 4);
            t +=  bvp.dt;
        }
    }

    for (int i = 0; i < N; i++) {
        printf("%0.5f\n", U1[i]);
    }
}