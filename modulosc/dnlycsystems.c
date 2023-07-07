#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

int sys1(double t, double x[], double f[], void *params) {
    double mu = *(double *)params;
    f[0] = mu * x[0] + x[0] * x[0];
    return GSL_SUCCESS;
}

int jac1(double t, double x[], double *dfdx, double dfdt[], void *params) {
    double mu = *(double *)params;
    gsl_matrix_view dfdx_mat = gsl_matrix_view_array(dfdx, 1, 1);
    gsl_matrix * J = &dfdx_mat.matrix;
    gsl_matrix_set(J, 0, 0, mu + 2.0 * x[0]);
    dfdt[0] = 0.0;
    return GSL_SUCCESS;
}

int sys2(double t, double x[], double f[], void *params) {
    double mu = *(double *)params;
    f[0] = - mu + x[0] * x[0];
    return GSL_SUCCESS;
}

int jac2(double t, double x[], double *dfdx, double dfdt[], void *params) {
    double mu = *(double *)params;
    gsl_matrix_view dfdx_mat = gsl_matrix_view_array(dfdx, 1, 1);
    gsl_matrix * J = &dfdx_mat.matrix;
    gsl_matrix_set(J, 0, 0, 2.0 * x[0]);
    dfdt[0] = 0.0;
    return GSL_SUCCESS;
}

int sys3(double t, double x[], double f[], void *params) {
    double mu = *(double *)params;
    f[0] = - mu * x[0] + x[0] * x[0] * x[0];
    return GSL_SUCCESS;
}

int jac3(double t, double x[], double *dfdx, double dfdt[], void *params) {
    double mu = *(double *)params;
    gsl_matrix_view dfdx_mat = gsl_matrix_view_array(dfdx, 1, 1);
    gsl_matrix *J = &dfdx_mat.matrix;
    gsl_matrix_set(J, 0, 0, -mu + 3.0 * x[0] * x[0]);
    dfdt[0] = 0.0;
    return GSL_SUCCESS;
}

int sys4(double t, double x[], double f[], void *params) {
    double mu = *(double *)params;
    f[0] = - mu * x[0] - x[0] * x[0] * x[0];
    return GSL_SUCCESS;
}

int jac4(double t, double x[], double *dfdx, double dfdt[], void *params) {
    double mu = *(double *)params;
    gsl_matrix_view dfdx_mat = gsl_matrix_view_array(dfdx, 1, 1);
    gsl_matrix *J = &dfdx_mat.matrix;
    gsl_matrix_set(J, 0, 0, -mu - 3.0 * x[0] * x[0]);
    dfdt[0] = 0.0;
    return GSL_SUCCESS;
}

int sys5(double t, double x[], double f[], void *params) {
    double mu = *(double *)params;
    f[0] = mu - x[0] * x[0] * x[0];
    return GSL_SUCCESS;
}

int jac5(double t, double x[], double *dfdx, double dfdt[], void *params) {
    double mu = *(double *)params;
    gsl_matrix_view dfdx_mat = gsl_matrix_view_array(dfdx, 1, 1);
    gsl_matrix *J = &dfdx_mat.matrix;
    gsl_matrix_set(J, 0, 0, - 3.0 * x[0] * x[0]);
    return GSL_SUCCESS;
}