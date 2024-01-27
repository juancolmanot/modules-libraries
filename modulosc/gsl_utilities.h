#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

/* Linspace generator */
gsl_vector *gsl_linspace(double x0, double xN, double N);

/* Quadratic regression */
void gsl_regression_quadratic(double *x_data, double *y_data, unsigned int n, double *coefficients);

/* Eigenvalues  of nxn square matrix A*/
gsl_vector *gsl_eigenvalues(gsl_matrix *A, int n);

/* Eigenvectors of nxn square matrix A*/
gsl_matrix *gsl_eigenvectors(gsl_matrix *A, int n);


#endif