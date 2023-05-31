#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "linear_algebra.h"

void iz_map(double* x1, double* x, double a, double b, double c, double d, double I) {

    double u = x[0];
    double v = x[1];

    double v1, u1;

    if (v < 30) {
        v1 = minof(0.04 * v * v + 6 * v + 140 - I - u, 30);
        u1 = u + a * (b * v - u);
    }

    else if (v >= 30) {
        v1 = c;
        u1 = u + d;
    }

    x1[0] = u1;
    x1[1] = v1;

}
