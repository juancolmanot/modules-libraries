#ifndef INTERMITTENCY_H
#define INTERMITTENCY_H

#include <stdbool.h>

bool reinjection(double xn1, double xn, double ubr);

bool ejection(double xn1, double xn, double ubr);

/* Two dimensional map from Kim & Kye */
void map_2d(double* x1, double* x, double alpha, double beta);

void map_2d_n(double* xn1, double* xn, unsigned int n, double alpha, double beta);

void map_2d_deacoplated(double* x1, double* x, double alpha, double beta);

void map_2d_n_dacoplated(double* xn1, double* xn, unsigned int n, double alpha, double beta);

void map_2d_rotated_u(double *xn1, double *xn, double alpha, double beta);

void map_2d_n_rotated_u(double *xn1, double *xn, unsigned int n, double alpha, double beta);

#endif