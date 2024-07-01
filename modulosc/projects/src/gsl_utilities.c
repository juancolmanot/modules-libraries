#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_blas.h>
#include "gsl_utilities.h"

void gsl_regression_quadratic(
    double xdat[],
    double ydat[],
    unsigned int n,
    double ci[]
) {
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, 3);
    gsl_matrix *X = gsl_matrix_alloc(n, 3);
    gsl_vector *y = gsl_vector_alloc(n);
    gsl_vector *coeff = gsl_vector_alloc(3);
    gsl_matrix *cov = gsl_matrix_alloc(3, 3);
    double chisq;

    for (size_t i = 0; i < n; ++i) {
        double xi = xdat[i];
        double xi2 = xi * xi;
        gsl_matrix_set(X, i, 0, 1.0);
        gsl_matrix_set(X, i, 1, xi);
        gsl_matrix_set(X, i, 2, xi2);
        gsl_vector_set(y, i, ydat[i]);
    }

    gsl_multifit_linear(X, y, coeff, cov, &chisq, work);
    ci[0] = gsl_vector_get(coeff, 0);
    ci[1] = gsl_vector_get(coeff, 1);
    ci[2] = gsl_vector_get(coeff, 2);

    gsl_multifit_linear_free(work);
    gsl_matrix_free(X);
    gsl_vector_free(y);
    gsl_matrix_free(cov);
    gsl_vector_free(coeff);
}