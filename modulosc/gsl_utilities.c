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
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include "gsl_utilities.h"

void gsl_linear_regression(
    double *x,
    double *y,
    unsigned int size,
    double *m,
    double *b
)
{
    double c0, c1, cov00, cov01, cov11, sumsq;

    gsl_fit_linear(x, 1, y, 1, size, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

    *m = c1;
    *b = c0;
}

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

double interp_funct(
    double xi,
    void *params
) {
    InterpData *data = (InterpData *)params;
    return gsl_interp_eval(data->interp, data->x, data->y, xi, data->acc);
}

double gsl_integrate_from_arrays(
    double xfit[],
    double yfit[],
    unsigned int n
) {
    // Create and initialize the interpolation object
    gsl_interp *interp = gsl_interp_alloc(gsl_interp_linear, n);
    gsl_interp_init(interp, xfit, yfit, n);

    // Create and initialize the accelerator
    gsl_interp_accel *acc = gsl_interp_accel_alloc();

    // Create InterpData struct
    InterpData data = {interp, acc, xfit, yfit, n};

    // Define the integration workspace
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(10000);

    // Define the integration function
    gsl_function F;
    F.function = &interp_funct;
    F.params = &data;

    // Perform the integration
    double result, error;
    int status = gsl_integration_qags(&F, xfit[0], xfit[n - 1], 0, 1e-4, 10000, workspace, &result, &error);

    if (status != GSL_SUCCESS) {
        fprintf(stderr, "GSL integration error: %s\n", gsl_strerror(status));
        // Handle error appropriately, e.g., by returning a special value or exiting
    }

    gsl_interp_free(interp);
    gsl_interp_accel_free(acc);
    gsl_integration_workspace_free(workspace);

    return result;
}