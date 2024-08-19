#ifndef ODEINTEGRATION_H
#define ODEINTEGRATION_H

typedef struct {
    long double t;    // current time
    long double h;    // current step size
    long double *x;   // current state vector
    int n;       // dimension of the state vector
} RKF45State;

typedef void (*deriv_func)(
    long double t, 
    long double x[],
    long double dxdt[],
    long double p[]
);

typedef void (*jacobian_func)(
    long double t, 
    long double x[],
    long double dxdx[],
    long double p[]
);

void rkf45_step(
    RKF45State *state,
    deriv_func f,
    jacobian_func jf,
    long double p[]
);

void rkf45_step_adap(
    RKF45State *state,
    deriv_func f,
    jacobian_func jf,
    long double p[],
    long double abs_err,
    long double rel_err
);

void rkf45_integrate(
    RKF45State *state,
    deriv_func f,
    jacobian_func jf,
    long double p[],
    long double t1,
    long double abs_err,
    long double rel_err
);

#endif