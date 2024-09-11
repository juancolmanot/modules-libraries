#ifndef GSL_UTILITIES_H
#define GSL_UTILITIES_H

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

typedef struct {
    gsl_interp *interp;
    gsl_interp_accel *acc;
    double *x;
    double *y;
    unsigned int n;
} InterpData;

double interp_funct(
    double xi,
    void *params
);

double gsl_integrate_from_arrays(
    double xfit[],
    double yfit[],
    unsigned int n
);

void gsl_linear_regression(
    double *x,
    double *y,
    unsigned int size,
    double *m,
    double *b
);

void gsl_regression_quadratic(
    double xdat[],
    double ydat[],
    unsigned int n,
    double coeff[]
);

#endif