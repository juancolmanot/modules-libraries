#ifndef STATS_H
#define STATS_H

void histogram_3d(
    long double **z,
    long double *x,
    long double *y,
    long double *datax,
    long double *datay,
    unsigned int data_size,
    unsigned int bin_size
);

void stats_histogram(
    long double *y,
    long double *x,
    long double *data,
    unsigned int data_size,
    unsigned int xsize
);

void stats_histogram_double(
    double *y,
    double *x,
    double *data,
    unsigned int data_size,
    unsigned int xsize
);

void linear_regression(
    long double *x,
    long double *y,
    unsigned int size,
    long double *m,
    long double *b
);

void quadratic_regression(
    long double xdat[],
    long double ydat[],
    unsigned int n,
    long double ci[]
);

#endif