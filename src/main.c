#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include "theta_scheme.h"

#define NUM_TERMS 5

double initial_condition(double x) {
    return x*(1-x);
}

double exact_solution(double x, double t) {
    double u = 0;
    for (int n = 1; n < 2*NUM_TERMS; n += 2) {
        double an = 8.0/pow(n*M_PI,3);
        u += an*sin(n*M_PI*x)*exp(-pow(M_PI*n,2)*t);
    }
    return u;
}

int main() {
    struct bvp_t bvp;
    bvp.J = 20;
    bvp.dt = 0.0012;
    bvp.theta = 0;
    bvp.IC = initial_condition;
    bvp.stop_time = 1;

    double err = theta_scheme(bvp, exact_solution);
    printf("max err = %0.9f\n", err);
}