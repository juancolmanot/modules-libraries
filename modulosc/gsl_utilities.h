#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H
#include <gsl/gsl_vector.h>

/* Linspace generator */
gsl_vector *gsl_linspace(double x0, double xN, double N);

/* Quadratic regression */
void gsl_regression_quadratic(double *x_data, double *y_data, unsigned int n, double *coefficients);

#endif