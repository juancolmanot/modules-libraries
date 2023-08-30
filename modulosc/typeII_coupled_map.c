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
} Parameters;

double LMBKx(double x, void *params) {
    Parameters *p = (Parameters*)params;
    double mu = p->mu;
    double gamma = p->gamma;
    double b = p->b;
    double x1;

    x1 = -((1 + mu + gamma) * x + x * x * x) * exp(-b * x * x) + gamma * x;
    
    return x1;
}

double LMBKy(double y, void *params) {
    Parameters *p = (Parameters*)params;
    double mu = p->mu;
    double gamma = p->gamma;
    double b = p->b;
    double y1;

    y1 = -((1 + mu + gamma) * y + y * y * y) * exp(-b * y * y) + gamma * y;
    
    return y1;
}

