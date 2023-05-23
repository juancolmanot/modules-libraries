#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double *linspace(double x0, double xn, unsigned int n) {

    double dx;
    double* x_space = malloc(n * sizeof(double));
    
    dx = (xn - x0) / (double)(n - 1);

    for (unsigned int i = 0; i < n; i++) {
        
        x_space[i] = x0 + dx * (double)i;
    }

    return x_space;
}

double* linspace_discrete(double x0, double xn, unsigned int n) {

    double* x_space = malloc(n * sizeof(double));

    for (unsigned int i = 0; i < n; i++){
        x_space[i] = x0 + i;
    }

    return x_space;
}

double min(double* x) {

    double minimum = x[0];

    for (unsigned int i = 0; i < sizeof x; i++) {
        if (x[i] < minimum) {
            minimum = x[i];
        }
    }

    return minimum;

}

double max(double* x) {

    double maximum = 0;

    for (unsigned int i = 0; i < sizeof x; i++) {
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