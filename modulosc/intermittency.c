#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

bool reinjection(double xn1, double xn, double ubr) {

    if (xn1 < ubr && xn > ubr){
        return true;
    }
    else {
        return false;
    }
}

bool ejection(double xn1, double xn, double ubr) {

    if (xn1 > ubr && xn < ubr) {
        return true;
    }
    else {
        return false;
    }
}

/* Two dimensional map from Kim & Kye */

void map_2d(double* x1, double* x, double alpha, double beta) {
    
    x1[0] = 4 * alpha * x[0] * (1 - x[0]) + beta * x[1] * (1 - x[0]);
    x1[1] = 4 * alpha * x[1] * (1 - x[1]) + beta * x[0] * (1 - x[1]);
}

void map_2d_n(double* x1, double* x, unsigned int n, double alpha, double beta) {

    double* x_aux = calloc(2, sizeof(double));
    
    x_aux[0] = x[0];
    x_aux[1] = x[1];

    for (unsigned int i = 0; i < n; i++){
        map_2d(x1, x_aux, alpha, beta);
        x_aux[0] = x1[0];
        x_aux[1] = x1[1];
    }
    free(x_aux);
}

void map_2d_deacoplated(double* x1, double* x, double alpha, double beta) {
    
    x1[0] = 4 * alpha * x[0] * (1 - x[0]); // + beta * x[1] * (1 - x[0]);
    x1[1] = 4 * alpha * x[1] * (1 - x[1]); // + beta * x[0] * (1 - x[1]);
}

void map_2d_n_deacoplated(double* x1, double* x, unsigned int n, double alpha, double beta) {

    double* x_aux = calloc(2, sizeof(double));
    
    x_aux[0] = x[0];
    x_aux[1] = x[1];

    for (unsigned int i = 0; i < n; i++){
        map_2d_deacoplated(x1, x_aux, alpha, beta);
        x_aux[0] = x1[0];
        x_aux[1] = x1[1];
    }
    free(x_aux);
}


/* Mapa con terminos cruzados anulados para variable u1 */
void map_2d_rotated_u(double *xn1, double *xn, double alpha, double beta) {

    double cu, cu2, cv2;
    double dv, duv;

    cu = 4 * alpha + beta;
    cu2 = - 2.8284271247461903 * alpha - 0.707107 * beta;
    cv2 = 0.707107 * beta - 2.8284271247461903 * alpha;

    dv = 4 * alpha - beta;
    duv = - 5.65685 * alpha;

    xn1[0] = cu * xn[0] + cu2 * xn[0] * xn[0] + cv2 * xn[1] * xn[1];
    xn1[1] = dv * xn[1] + duv * xn[0] * xn[1];
}

void map_2d_n_rotated_u(double *xn1, double *xn, unsigned int n, double alpha, double beta) {
    
    double *xn_aux = calloc(2, sizeof(double));

    xn_aux[0] = xn[0];
    xn_aux[1] = xn[1];

    for (unsigned int i = 0; i < n; i++) {
        map_2d_rotated_u(xn1, xn_aux, alpha, beta);
        xn_aux[0] = xn1[0];
        xn_aux[1] = xn1[1];
    }
    free(xn_aux);
}