#pragma once

/**
 * Algorithms to solve the tridiagonal system Ax_2 = x_1 where
 *     | (1+2a)   -a                        |
 *     |   -a   (1+2a)  -a                  |
 * A = |          -a   (1+2a)  -a           |
 *     |                  .....             |
 *     |                        -a   (1+2a) |
 *  
 * This specific form of matrix arises from the backward Euler method.
 * 
 * @param a   The parameter a that defines the tridiagonal matrix A
 * @param x1  The known right-hand side
 * @param x2  The unknown to solve for
 * @param N   The dimension of the the system
 * @param n   The number of processors to parallelize the computation over
 */
void implicit_step_serial(double a, double *x1, double *x2, int N);
void implicit_step_parallel(double a, double *x1, double *x2, int N, int n);