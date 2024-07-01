#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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