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

void error_analysis(struct bvp_t bvp, char *output, double mu, double nu) {
    int max_J = 80;
    FILE *fp = fopen(output, "w"); 
    for (bvp.J = 10; bvp.J <= max_J; bvp.J += 2) {
        double dx = 1.0/bvp.J;
        bvp.dt = fmax(mu*dx*dx, nu*dx);
        double err = theta_scheme(bvp, exact_solution);
        fprintf(fp, "%ld\t%f\t%0.9f\n", bvp.J, bvp.J/bvp.dt, log10(err));
    }
    fclose(fp);
}

int main() {
    int max_J = 80;
    { // theta = 0, mu = 1/2
        double mu = 0.5;
        struct bvp_t bvp;
        bvp.IC = initial_condition;
        bvp.stop_time = 1;
        bvp.theta = 0;
        error_analysis(bvp, "analysis/1_err.txt", mu, 0);

        bvp.J = max_J;
        bvp.dt = mu/(bvp.J*bvp.J);
        bvp.stop_time = 0.5;
        plot_solution(bvp, exact_solution, "analysis/1_sol.txt");
    }

    printf("\n");
    { // theta = 1, mu = 5
        double mu = 5;
        struct bvp_t bvp;
        bvp.IC = initial_condition;
        bvp.stop_time = 1;
        bvp.theta = 1;
        error_analysis(bvp, "analysis/2_err.txt", mu, 0);

        bvp.J = max_J;
        bvp.dt = mu/(bvp.J*bvp.J);
        bvp.stop_time = 0.5;
        plot_solution(bvp, exact_solution, "analysis/2_sol.txt");
    }

    printf("\n");
    { // theta = 0.5, nu = 1/20
        double nu = 0.05;
        struct bvp_t bvp;
        bvp.IC = initial_condition;
        bvp.stop_time = 1;
        bvp.theta = 0.5;
        error_analysis(bvp, "analysis/3_err.txt", 0, nu);

        bvp.J = max_J;
        bvp.dt = nu/bvp.J;
        bvp.stop_time = 0.5;
        plot_solution(bvp, exact_solution, "analysis/3_sol.txt");
    }

}