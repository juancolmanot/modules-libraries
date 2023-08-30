#ifndef LORENZ_H
#define LORENZ_H

typedef struct {
    double a;
    double b;
    double c;
} Parameters;


int lorenz(double t, double x[], double f[], void *params);

int jaclorenz(double t, double x[], double *dfdx, double dfdt[], void *params);

void lorenzH(double x[], double f[], double *U, double *K);

#endif