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

double min(double* x) {

    double minimum = 0;

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