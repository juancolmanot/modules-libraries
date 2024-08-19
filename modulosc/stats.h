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

void linear_regression_double(
    double *x,
    double *y,
    unsigned int size,
    double *m,
    double *b
);

void quadratic_regression(
    long double xdat[],
    long double ydat[],
    unsigned int n,
    long double ci[]
);

void bubble_sort(
    long double *arr,
    int n
);

void bubble_sort_double_unsigned(
    double *arr,
    unsigned int n
);

//===================================================
// Swap function for quick sort algorithm.
//===================================================
void swap(
    double *a,
    double *b
);

void swap_long(
    long double *a,
    long double *b
);

//===================================================
// Partition function for quick sort algorithm.
//===================================================
int partition(
    double arr[],
    int low,
    int high
);

int partition_long(
    long double arr[],
    int low,
    int high
);

//===================================================
// Quicksort algorithm function.
//===================================================
void quicksort(
    double arr[],
    int low,
    int high
);

void quicksort_long(
    long double arr[],
    int low,
    int high
);

//===================================================
// Wrapper function to call quicksort.
//===================================================
void quicksort_double_unsigned(
    double *arr,
    unsigned int n
);

void quicksort_long_unsigned(
    long double *arr,
    unsigned int n
);

//===================================================
// MONTECARLO integrator.
//===================================================
double montecarlo_integration(
    double x[],
    double fx[],
    unsigned int n,
    unsigned int n_samples
);

long double montecarlo_integration_long(
    long double x[],
    long double fx[],
    unsigned int n,
    unsigned int n_samples
);

#endif