#pragma once

// Numerical parameters for the BVP u_t = u_xx
struct bvp_t {
    size_t J;             // Number of spatial grid divisions
    double (*IC)(double); // Initial condition
    double dt;            // Timestep size
    double stop_time;
    double theta;
};

void explicit_step(double a, double *U1, double *U2, int N);
double theta_scheme(struct bvp_t bvp, double (*u)(double,double));