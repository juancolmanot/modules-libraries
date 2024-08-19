#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "stats.h"
#include "linear_algebra.h"

void histogram_3d(
    long double **z,
    long double *x,
    long double *y,
    long double *datax,
    long double *datay,
    unsigned int data_size,
    unsigned int bin_size
){
    long double xmin, xmax, ymin, ymax;
    xmin = la_min(datax, data_size);
    xmax = la_max(datax, data_size);
    ymin = la_min(datay, data_size);
    ymax = la_max(datay, data_size);

    long double dx, dy;
    dx = (xmax - xmin) / (long double)(bin_size - 1);
    dy = (ymax - ymin) / (long double)(bin_size - 1);
    
    for (unsigned int i = 0; i < bin_size; i++){
        x[i] = xmin + dx * (long double)i;
        y[i] = ymin + dy * (long double)i;
    }

    for (unsigned int k = 0; k < data_size; k++){
        // Find the bin index for datax[k]
        unsigned int i = (unsigned int)((datax[k] - xmin) / dx);
        if (i >= bin_size) i = bin_size - 1;

        // Find the bin index for datay[k]
        unsigned int j = (unsigned int)((datay[k] - ymin) / dy);
        if (j >= bin_size) j = bin_size - 1;

        // Increment the bin count
        z[i][j]++;
    }
}

void stats_histogram(
    long double *y,
    long double *x,
    long double *data,
    unsigned int data_size,
    unsigned int xsize
){
    long double xmax = 0, xmin = 0;
    xmax = la_max(data, data_size);
    xmin = la_min(data, data_size);
    
    long double dx = (xmax - xmin) / (long double) (xsize - 1);

    for (unsigned int i = 0; i < xsize - 1; i++) {
        for (unsigned int j = 0; j < data_size; j++) {
            if (data[j] > (xmin + dx * i) && data[j] < (xmin + dx * (i + 1))) {
                y[i]++;
            }
        }
        x[i] = xmin + dx * i;
    }
}

void stats_histogram_double(
    double *y,
    double *x,
    double *data,
    unsigned int data_size,
    unsigned int xsize
){
    double xmax = 0, xmin = 0;
    xmax = la_max_d(data, data_size);
    xmin = la_min_d(data, data_size);
    
    double dx = (xmax - xmin) / (double) (xsize - 1);

    for (unsigned int i = 0; i < xsize - 1; i++) {
        for (unsigned int j = 0; j < data_size; j++) {
            if (data[j] > (xmin + dx * i) && data[j] < (xmin + dx * (i + 1))) {
                y[i]++;
            }
        }
        x[i] = xmin + dx * i;
    }
}

void linear_regression(
    long double *x,
    long double *y,
    unsigned int size,
    long double *m,
    long double *b
){
    long double sumx, sumy, sumxy, sumx2;
    sumx = sumy = sumxy = sumx2 = 0;

    for (unsigned int i = 0; i < size; i++){
        sumx += x[i];
        sumy += y[i];
        sumxy += x[i] * y[i];
        sumx2 += x[i] * x[i];
    }
    
    *m = ((long double)size * sumxy - sumx * sumy) / ((long double)size * sumx2 - sumx * sumx);
    *b = (sumy - *m * sumx) / (long double)size;
}

void linear_regression_double(
    double *x,
    double *y,
    unsigned int size,
    double *m,
    double *b
){
    double sumx, sumy, sumxy, sumx2;
    sumx = sumy = sumxy = sumx2 = 0;
    for (unsigned int i = 0; i < size; i++){
        sumx += x[i];
        sumy += y[i];
        sumxy += x[i] * y[i];
        sumx2 += x[i] * x[i];
    }
    
    *m = ((double)size * sumxy - sumx * sumy) / ((double)size * sumx2 - sumx * sumx);
    *b = (sumy - *m * sumx) / (double)size;
}

void quadratic_regression(
    long double xdat[],
    long double ydat[],
    unsigned int n,
    long double ci[]
){
    long double Sx = 0.0, Sy = 0.0, Sxx = 0.0, Sxy = 0.0, Sxxx = 0.0, Sxxxx = 0.0, Sxxy = 0.0;
    for (unsigned int i = 0; i < n; i++) {
        long double x = xdat[i];
        long double y = ydat[i];
        Sx += x;
        Sy += y;
        Sxx += x * x;
        Sxy += x * y;
        Sxxx += x * x * x;
        Sxxxx += x * x * x * x;
        Sxxy += x * x * y;
    }

    long double det = Sxxxx * (Sxx * n - Sx * Sx) - Sxxx * (Sxxx * n - Sx * Sxx) + Sxx * (Sxxx * Sx - Sxx * Sxx);
    long double detA = Sxxy * (Sxx * n - Sx * Sx) - Sxxx * (Sxy * n - Sx * Sy) + Sxx * (Sxy * Sx - Sxx * Sy);
    long double detB = Sxxxx * (Sxy * n - Sx * Sy) - Sxxy * (Sxxx * n - Sx * Sxx) + Sxx * (Sxxx * Sy - Sxx * Sxy);
    long double detC = Sxxxx * (Sxx * Sy - Sx * Sxy) - Sxxx * (Sxxx * Sy - Sx * Sxxy) + Sxxy * (Sxxx * Sx - Sxx * Sxx);

    ci[0] = detA / det; // coefficient for x^2 (a)
    ci[1] = detB / det; // coefficient for x (b)
    ci[2] = detC / det; // constant term (c)
}

void bubble_sort(
    long double *arr,
    int n
) {
    int i, j;
    long double temp;
    for (i = 0; i < n-1; i++) {
        for (j = 0; j < n-i-1; j++) {
            if (arr[j] > arr[j+1]) {
                // Swap arr[j] and arr[j+1]
                temp = arr[j];
                arr[j] = arr[j+1];
                arr[j+1] = temp;
            }
        }
    }
}

void bubble_sort_double_unsigned(
    double *arr,
    unsigned int n
) {
    unsigned int i, j;
    double temp;
    for (i = 0; i < n-1; i++) {
        for (j = 0; j < n-i-1; j++) {
            if (arr[j] > arr[j+1]) {
                // Swap arr[j] and arr[j+1]
                temp = arr[j];
                arr[j] = arr[j+1];
                arr[j+1] = temp;
            }
        }
    }
}

//===================================================
// Swap function for quick sort algorithm.
//===================================================
void swap(
    double *a,
    double *b
)
{
    double temp = *a;
    *a = *b;
    *b = temp;
}

void swap_long(
    long double *a,
    long double *b
)
{
    long double temp = *a;
    *a = *b;
    *b = temp;
}

//===================================================
// Partition function for quick sort algorithm.
//===================================================
int partition(
    double arr[],
    int low,
    int high
)
{
    double pivot = arr[high];
    int i = low - 1;
    for (int j = low; j <= high - 1; j++) {
        if (arr[j] < pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

int partition_long(
    long double arr[],
    int low,
    int high
)
{
    long double pivot = arr[high];
    int i = low - 1;
    for (int j = low; j <= high - 1; j++) {
        if (arr[j] < pivot) {
            i++;
            swap_long(&arr[i], &arr[j]);
        }
    }
    swap_long(&arr[i + 1], &arr[high]);
    return (i + 1);
}

//===================================================
// Quicksort algorithm function.
//===================================================
void quicksort(
    double arr[],
    int low,
    int high
)
{
    if (low < high) {
        int pi = partition(arr, low, high);
        quicksort(arr, low, pi - 1);
        quicksort(arr, pi + 1, high);
    }
}

void quicksort_long(
    long double arr[],
    int low,
    int high
)
{
    if (low < high) {
        int pi = partition_long(arr, low, high);
        quicksort_long(arr, low, pi - 1);
        quicksort_long(arr, pi + 1, high);
    }
}

//===================================================
// Wrapper function to call quicksort.
//===================================================
void quicksort_double_unsigned(
    double *arr,
    unsigned int n
)
{
    quicksort(arr, 0, (int)(n - 1));
}

void quicksort_long_unsigned(
    long double *arr,
    unsigned int n
)
{
    quicksort_long(arr, 0, (int)(n - 1));
}

//===================================================
// MONTECARLO integrator.
//===================================================
double montecarlo_integration(
    double x[],
    double fx[],
    unsigned int n,
    unsigned int n_samples
)
{
    // Initialize random number generator
    srand((unsigned int)time(NULL));

    // Determine the range of integration
    double x_min = x[0];
    double x_max = x[n - 1];

     // Sum of function values at random points
    double sum_fx = 0.0;

    for (unsigned int i = 0; i < n_samples; i++) {
        // Generate a random x value within the range
        double x_random = x_min + (x_max - x_min) * ((double)rand() / RAND_MAX);

        // Interpolate the function value at the random x
        double f_random = 0.0;
        for (unsigned int j = 0; j < n - 1; j++) {
            if (x_random >= x[j] && x_random <= x[j + 1]) {
                double t = (x_random - x[j]) / (x[j + 1] - x[j]);
                f_random = (1 - t) * fx[j] + t * fx[j + 1];
                break;
            }
        }

        // Accumulate the function values
        if (isnanl(f_random) != 1 && isinfl(f_random) != 1){
            sum_fx += f_random;
        }
    }

    // Estimate the integral
    double integral = (x_max - x_min) * sum_fx / n_samples;

    return integral;
}

long double montecarlo_integration_long(
    long double x[],
    long double fx[],
    unsigned int n,
    unsigned int n_samples
)
{
    // Initialize random number generator
    srand((unsigned int)time(NULL));

    // Determine the range of integration
    long double x_min = x[0];
    long double x_max = x[n - 1];

    // Sum of function values at random points
    long double sum_fx = 0.0;

    unsigned int actualcount = 0;

    for (unsigned int i = 0; i < n_samples; i++) {
        // Generate a random x value within the range
        long double x_random = x_min + (x_max - x_min) * ((long double)rand() / RAND_MAX);

        // Interpolate the function value at the random x
        long double f_random = 0.0;
        for (unsigned int j = 0; j < n - 1; j++) {
            if (x_random >= x[j] && x_random <= x[j + 1]) {
                long double t = (x_random - x[j]) / (x[j + 1] - x[j]);
                f_random = (1 - t) * fx[j] + t * fx[j + 1];
                break;
            }
        }

        // Accumulate the function values

        if (f_random < 1e7){
            sum_fx += f_random;
            actualcount++;
        }
    }

    // Estimate the integral
    long double integral = (x_max - x_min) * sum_fx / (long double)actualcount;

    return integral;
}