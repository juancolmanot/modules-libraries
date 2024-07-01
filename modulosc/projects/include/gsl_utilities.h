#ifndef GSL_UTILITIES_H
#define GSL_UTILITIES_H

void gsl_regression_quadratic(
    double xdat[],
    double ydat[],
    unsigned int n,
    double coeff[]
);

#endif