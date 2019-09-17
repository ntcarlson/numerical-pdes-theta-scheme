#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "theta_scheme.h"

double initial_condition(double x) {
    return x*(1-x);
}

int main() {
    struct bvp_t bvp;
    bvp.J = 20;
    bvp.dt = 0.0012;
    bvp.theta = 0.5;
    bvp.IC = initial_condition;
    bvp.stop_time = 0.05;

    theta_scheme(bvp);
}