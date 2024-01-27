#ifndef LORENZ_H
#define LORENZ_H
#include <gsl/gsl_matrix.h>

typedef struct {
    double a;
    double b;
    double c;
} Parameters;


int lorenz(double t, double x[], double f[], void *params);

gsl_matrix *lorenz_J(double x[], void *params);

int jaclorenz(double t, double x[], double *dfdx, double dfdt[], void *params);

void lorenzH(double x[], double f[], double *U, double *K);

gsl_matrix *fundamental(double x[], void *params);

#endif