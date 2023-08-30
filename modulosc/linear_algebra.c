#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void distance_line(double *x, double *x1, double *xerr) {

    xerr[0] = fabs(x[0] - x1[0]) / sqrt(2);
    xerr[1] = fabs(x[1] - x1[1]) / sqrt(2);
}

double *linspace(double x0, double xn, unsigned int n) {

    double dx;
    double* x_space = malloc(n * sizeof(double));
    
    dx = (xn - x0) / (double)(n - 1);

    for (unsigned int i = 0; i < n; i++) {
        
        x_space[i] = x0 + dx * (double)i;
    }

    return x_space;
    free(x_space);
}

double* logspace(double exp0, double expn, unsigned int n) {

    double delta_exp = (expn - exp0) / (double) (n - 1);
    double *x_logspace = calloc(n, sizeof(double));

    double exp_i = exp0;
    
    for (unsigned int i = 0; i < n; i++) {
        x_logspace[i] = pow(10, exp_i);
        exp_i += delta_exp;
    }

    return x_logspace;
    free(x_logspace);
}

double* linspace_discrete(double x0, double xn, unsigned int n) {

    double* x_space = calloc(n, sizeof(double));

    for (unsigned int i = 0; i < n; i++){
        x_space[i] = x0 + i;
    }

    return x_space;
    free(x_space);
}

void rotation_operator(double *xrot, double *x, double theta) {
    
    xrot[0] = cos(theta) * x[0] - sin(theta) * x[1];
    xrot[1] = sin(theta) * x[0] + cos(theta) * x[1];
}

void antirotation_operator(double *xrot, double *x, double theta) {
    
    xrot[0] = cos(theta) * x[0] + sin(theta) * x[1];
    xrot[1] = - sin(theta) * x[0] + cos(theta) * x[1];
}

double min(double* x, unsigned int length) {

    double minimum = x[0];

    for (unsigned int i = 0; i < length; i++) {
        if (x[i] < minimum) {
            minimum = x[i];
        }
    }

    return minimum;

}

double max(double* x, unsigned int length) {

    double maximum = 0;

    for (unsigned int i = 0; i < length; i++) {
        if (x[i] > maximum) {
            maximum = x[i];
        }
    }
    
    return maximum;
}

double distance(double* x) {

    double d;

    d = sqrt(x[0] * x[0] + x[1] * x[1]);

    return d;
}

void rel_err(double* x, double* x1, double* xerr) {

    xerr[0] = fabs((x1[0] - x[0]) / x[0]);
    xerr[1] = fabs((x1[1] - x[1]) / x[1]);
}

double rel_err_scalar(double x, double x1) {

    double err;

    err = fabs((x1 - x) / x);
    
    return err;
}