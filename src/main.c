#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include "theta_scheme.h"

#define NUM_TERMS 4

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
    int max_J = 200;
    struct bvp_t bvp;
    bvp.IC = initial_condition;
    bvp.stop_time = 1;
    { // theta = 0, mu = 1/2
        FILE *fp = fopen("analysis/1.txt", "w"); 
        bvp.theta = 0;
        double mu = 0.5;
        for (bvp.J = 10; bvp.J <= max_J; bvp.J += 10) {
            double dx = 1.0/((double) bvp.J);
            bvp.dt = mu*dx*dx;
            double err = theta_scheme(bvp, exact_solution);
            fprintf(fp, "%d\t%0.9f\n", bvp.J, log10(err));
        }
        fclose(fp);
    }

    { // theta = 1, mu = 5
        FILE *fp = fopen("analysis/2.txt", "w"); 
        bvp.theta = 1;
        double mu = 5;
        for (bvp.J = 10; bvp.J <= max_J; bvp.J += 10) {
            double dx = 1.0/((double) bvp.J);
            bvp.dt = mu*dx*dx;
            double err = theta_scheme(bvp, exact_solution);
            fprintf(fp, "%d\t%0.9f\n", bvp.J, log10(err));
        }
        fclose(fp);
    }

    { // theta = 0.5, nu = 1/20
        FILE *fp = fopen("analysis/3.txt", "w"); 
        bvp.theta = 0.5;
        double nu = 0.05;
        for (bvp.J = 10; bvp.J <= max_J; bvp.J += 10) {
            double dx = 1.0/((double) bvp.J);
            bvp.dt = nu*dx;
            double err = theta_scheme(bvp, exact_solution);
            fprintf(fp, "%d\t%0.9f\n", bvp.J, log10(err));
        }
        fclose(fp);
    }

}