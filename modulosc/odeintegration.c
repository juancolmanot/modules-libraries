#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "odeintegration.h"

// RKF45 step function
void rkf45_step(
    RKF45State *state,
    deriv_func f,
    jacobian_func jf,
    long double p[]
) {
    int i;
    long double t = state->t;
    long double h = state->h;
    int n = state->n;
    long double x[n];
    for (int j = 0; j < n; j++) {
        x[j] = state->x[j];
    }

    long double k1[n], k2[n], k3[n], k4[n], k5[n], k6[n];
    long double x_temp[n];
    long double J[n][n];

    // Coefficients for RKF45
    const long double a2 = 1.0 / 4.0;
    const long double a3 = 3.0 / 8.0;
    const long double a4 = 12.0 / 13.0;
    const long double a5 = 1.0;
    const long double a6 = 1.0 / 2.0;

    const long double b21 = 1.0 / 4.0;
    const long double b31 = 3.0 / 32.0, b32 = 9.0 / 32.0;
    const long double b41 = 1932.0 / 2197.0, b42 = -7200.0 / 2197.0, b43 = 7296.0 / 2197.0;
    const long double b51 = 439.0 / 216.0, b52 = -8.0, b53 = 3680.0 / 513.0, b54 = -845.0 / 4104.0;
    const long double b61 = -8.0 / 27.0, b62 = 2.0, b63 = -3544.0 / 2565.0, b64 = 1859.0 / 4104.0, b65 = -11.0 / 40.0;

    const long double c1 = 16.0 / 135.0, c3 = 6656.0 / 12825.0, c4 = 28561.0 / 56430.0, c5 = -9.0 / 50.0, c6 = 2.0 / 55.0;

    
    // Compute the stages
    f(t, x, k1, p);
    for (i = 0; i < n; i++) x_temp[i] = x[i] + h * b21 * k1[i];
    f(t + a2 * h, x_temp, k2, p);
    for (i = 0; i < n; i++) x_temp[i] = x[i] + h * (b31 * k1[i] + b32 * k2[i]);
    f(t + a3 * h, x_temp, k3, p);
    for (i = 0; i < n; i++) x_temp[i] = x[i] + h * (b41 * k1[i] + b42 * k2[i] + b43 * k3[i]);
    f(t + a4 * h, x_temp, k4, p);
    for (i = 0; i < n; i++) x_temp[i] = x[i] + h * (b51 * k1[i] + b52 * k2[i] + b53 * k3[i] + b54 * k4[i]);
    f(t + a5 * h, x_temp, k5, p);
    for (i = 0; i < n; i++) x_temp[i] = x[i] + h * (b61 * k1[i] + b62 * k2[i] + b63 * k3[i] + b64 * k4[i] + b65 * k5[i]);
    f(t + a6 * h, x_temp, k6, p);

    // Compute the jacobian Matrix
    jf(t, x, J, p);

    // Compute the next value of y and the error estimate
    for (i = 0; i < n; i++) {
        x_temp[i] = x[i] + h * (c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c5 * k5[i] + c6 * k6[i]);
    }

    state->t += h;
    for (i = 0; i < n; i++) {
        state->x[i] = x_temp[i];
    }
}

// RKF45 step function
void rkf45_step_adap(
    RKF45State *state,
    deriv_func f,
    jacobian_func jf,
    long double p[],
    long double abs_err,
    long double rel_err
) {
    int i;
    long double t, h;
    int n = state->n;
    long double x[n];
    long double k1[3], k2[3], k3[3], k4[3], k5[3], k6[3];
    long double x_temp[3], x_err[3];

    // Coefficients for RKF45
    const long double a2 = 1.0 / 4.0;
    const long double a3 = 3.0 / 8.0;
    const long double a4 = 12.0 / 13.0;
    const long double a5 = 1.0;
    const long double a6 = 1.0 / 2.0;

    const long double b21 = 1.0 / 4.0;
    const long double b31 = 3.0 / 32.0, b32 = 9.0 / 32.0;
    const long double b41 = 1932.0 / 2197.0, b42 = -7200.0 / 2197.0, b43 = 7296.0 / 2197.0;
    const long double b51 = 439.0 / 216.0, b52 = -8.0, b53 = 3680.0 / 513.0, b54 = -845.0 / 4104.0;
    const long double b61 = -8.0 / 27.0, b62 = 2.0, b63 = -3544.0 / 2565.0, b64 = 1859.0 / 4104.0, b65 = -11.0 / 40.0;

    const long double c1 = 16.0 / 135.0, c3 = 6656.0 / 12825.0, c4 = 28561.0 / 56430.0, c5 = -9.0 / 50.0, c6 = 2.0 / 55.0;
    const long double dc1 = c1 - 25.0 / 216.0, dc3 = c3 - 1408.0 / 2565.0, dc4 = c4 - 2197.0 / 4104.0, dc5 = c5 + 1.0 / 5.0;

    
    while (1) {
        t = state->t;
        h = state->h;

        for (i = 0; i < n; i++) {
            x[i] = state->x[i];
        }

        // Compute the stages
        f(t, x, k1, p);
        for (i = 0; i < n; i++) x_temp[i] = x[i] + h * b21 * k1[i];
        f(t + a2 * h, x_temp, k2, p);
        for (i = 0; i < n; i++) x_temp[i] = x[i] + h * (b31 * k1[i] + b32 * k2[i]);
        f(t + a3 * h, x_temp, k3, p);
        for (i = 0; i < n; i++) x_temp[i] = x[i] + h * (b41 * k1[i] + b42 * k2[i] + b43 * k3[i]);
        f(t + a4 * h, x_temp, k4, p);
        for (i = 0; i < n; i++) x_temp[i] = x[i] + h * (b51 * k1[i] + b52 * k2[i] + b53 * k3[i] + b54 * k4[i]);
        f(t + a5 * h, x_temp, k5, p);
        for (i = 0; i < n; i++) x_temp[i] = x[i] + h * (b61 * k1[i] + b62 * k2[i] + b63 * k3[i] + b64 * k4[i] + b65 * k5[i]);
        f(t + a6 * h, x_temp, k6, p);

        // Compute the next value of y and the error estimate
        for (i = 0; i < n; i++) {
            x_err[i] = h * (dc1 * k1[i] + dc3 * k3[i] + dc4 * k4[i] + dc5 * k5[i] + c6 * k6[i]);
            x_temp[i] = x[i] + h * (c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c5 * k5[i] + c6 * k6[i]);
        }

        // Estimate the error
        long double err = 0.0;
        for (i = 0; i < n; i++) {
            long double sk = abs_err + rel_err * fmaxl(fabsl(x[i]), fabsl(x_temp[i]));
            err += (x_err[i] / sk) * (x_err[i] / sk);
        }

        err = sqrtl(err / (long double)n);

        // Compute the optimal step size
        long double safety = 0.9;
        long double beta = 0.5;
        long double min_factor = 0.2;
        long double max_factor = 10.0;
        long double h_new = h * fmaxl(min_factor, fminl(max_factor, safety * powl(1.0 / err, beta)));
        
        // Accept the step if the error is within the tolerance
        if (err <= 1.0) {
            state->t += h;
            for (i = 0; i < n; i++) {
                state->x[i] = x_temp[i];
            }
            state->h = h_new;
            break;
        } else {
            // If the step is rejected, reduce the step size
            state->h = h_new;
        }
    }
}

void rkf45_integrate(
    RKF45State *state,
    deriv_func f,
    jacobian_func jf,
    long double p[],
    long double t1,
    long double abs_err,
    long double rel_err
) {
    while (state->t < t1) {
        if (state->t + state->h > t1) {
            state->h = t1 - state->t;
        }
        rkf45_step_adap(state, f, p, abs_err, rel_err);
    }
}