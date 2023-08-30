#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>

gsl_vector *gsl_linspace(double x0, double xN, double N){
    if (N <= 1) {
        fprintf(stderr, "Error: Number of points (N) should be greater than 1.\n");
        return NULL;
    }
    if (xN <= x0) {
        fprintf(stderr, "Error: xN should be greater than x0.\n");
        return NULL;
    }
    
    gsl_vector *x_space = gsl_vector_alloc(N);
    double dx = (double)(xN - x0) / (N - 1);
    for (unsigned int i = 0; i < N; i++){
        gsl_vector_set(x_space, i, x0 + (double) i * dx);
    }
    return x_space; 
}

void gsl_regression_quadratic(double *x_data, double *y_data, unsigned int n, double *coefficients) {
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, 3);
    gsl_matrix *X = gsl_matrix_alloc(n, 3);
    gsl_vector *y = gsl_vector_alloc(n);
    gsl_vector *coeff = gsl_vector_alloc(3);
    gsl_matrix *cov = gsl_matrix_alloc(3, 3);
    double chisq;

    for (size_t i = 0; i < n; ++i) {
        double xi = x_data[i];
        double xi2 = xi * xi;
        gsl_matrix_set(X, i, 0, 1.0);
        gsl_matrix_set(X, i, 1, xi);
        gsl_matrix_set(X, i, 2, xi2);
        gsl_vector_set(y, i, y_data[i]);
    }

    gsl_multifit_linear(X, y, coeff, cov, &chisq, work);
    coefficients[0] = gsl_vector_get(coeff, 0);
    coefficients[1] = gsl_vector_get(coeff, 1);
    coefficients[2] = gsl_vector_get(coeff, 2);

    gsl_multifit_linear_free(work);
    gsl_matrix_free(X);
    gsl_vector_free(y);
    gsl_matrix_free(cov);
    gsl_vector_free(coeff);
}