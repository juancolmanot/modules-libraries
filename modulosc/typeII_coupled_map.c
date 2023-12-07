#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

typedef struct {
    double mu;
    double gamma;
    double b;
    double eps;
} Parameters;

double *LMBK(double *x, void *params) {
    Parameters *p = (Parameters*)params;
    double mu = p->mu;
    double gamma = p->gamma;
    double b = p->b;
    double eps = p->eps;
    double *x1;
    x1 = calloc(2, sizeof(double));
    x1[0] = -((1 + mu + gamma) * x[0] + x[0] * x[0] * x[0]) * exp(-b * x[0] * x[0]) + gamma * x[0] + eps * x[1];
    x1[1] = -((1 + mu + gamma) * x[1] + x[1] * x[1] * x[1]) * exp(-b * x[1] * x[1]) + gamma * x[1] - eps * x[0];
    return x1;
}