#ifndef INTERMITTENCY_H
#define INTERMITTENCY_H

#include <stdbool.h>

bool reinjection(double xn1, double xn, double ubr, double lbr);

bool ejection(double xn1, double xn, double ubr, double lbr);

/* Two dimensional map from Kim & Kye */
void map_2d(double* x1, double* x, double alpha, double beta);

void map_2d_n(double* xn1, double* xn, unsigned int n, double alpha, double beta);

#endif