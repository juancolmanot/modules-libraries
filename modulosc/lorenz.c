#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

typedef struct {
    double a;
    double b;
    double c;
} Parameters;


int lorenz(double t, double x[], double f[], void *params) {
    Parameters *p = (Parameters*)params;
    double sigma = p->a;
    double rho = p->b;
    double beta = p->c;
    
    f[0] = sigma * (x[1] - x[0]);
    f[1] = ((rho - x[2]) * x[0]) - x[1];
    f[2] = (x[0] * x[1]) - beta * x[2];

    return GSL_SUCCESS; 
}

int jaclorenz(double t, double x[], double *dfdx, double dfdt[], void *params) {
    Parameters *p = (Parameters*)params;
    double sigma = p->a;
    double rho = p->b;
    double beta = p->c;

    gsl_matrix_view dfdx_mat = gsl_matrix_view_array(dfdx, 3, 3);
    gsl_matrix *J = &dfdx_mat.matrix;
    gsl_matrix_set(J, 0, 0, - sigma);
    gsl_matrix_set(J, 0, 1, sigma);
    gsl_matrix_set(J, 0, 2, 0.0);
    gsl_matrix_set(J, 1, 0, rho - x[2]);
    gsl_matrix_set(J, 1, 1, - 1.0);
    gsl_matrix_set(J, 1, 2, - x[0]);
    gsl_matrix_set(J, 2, 0, x[1]);
    gsl_matrix_set(J, 2, 1, x[0]);
    gsl_matrix_set(J, 2, 2, - beta);

    dfdt[0] = dfdt[1] = dfdt[2] = 0.0;
    return GSL_SUCCESS;
}

void lorenzH(double x[], double f[], double *U, double *K) {

    double Ki, Ui;

    Ki = 0.5 * (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
    Ui = 0.5 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);

    *U = Ui;
    *K = Ki;
}