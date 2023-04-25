#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

bool reinjection(double xn1, double xn, double ubr, double lbr) {

    if (xn1 < ubr && xn1 > lbr && xn > ubr && xn < lbr){
        return true;
    }
    else {
        return false;
    }
}

bool ejection(double xn1, double xn, double ubr, double lbr) {

    if (xn1 > ubr && xn1 < lbr && xn < ubr && xn > lbr) {
        return true;
    }
    else {
        return false;
    }
}

/* Two dimensional map from Kim & Kye */

void map_2d(double* x1, double* x, double alpha, double beta) {
    
    x1[0] = 4.f * alpha * x[0] * (1.f - x[0]) + beta * x[1] * (1.f - x[0]);
    x1[1] = 4.f * alpha * x[1] * (1.f - x[1]) + beta * x[0] * (1.f - x[1]);
}

void map_2d_n(double* x1, double* x, unsigned int n, double alpha, double beta) {

    double* x_aux = malloc(2 * sizeof(double));
    
    x_aux[0] = x[0];
    x_aux[1] = x[1];

    for (unsigned int i = 0; i < n; i++){
        map_2d(x1, x_aux, alpha, beta);
        x_aux[0] = x1[0];
        x_aux[1] = x1[1];
    }
}