#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "calculus.h"
#include "stats.h"

double *local_slope(
	double *y,
	double *x,
	unsigned int size,
	unsigned int sample
)
{	
	unsigned int size_slopes = size;

	if (sample % 2 == 0){
		size_slopes = size - sample + 2;
	}
	else {
		size_slopes = size - sample + 1;
	}

	double *slopes = malloc(size_slopes * sizeof(double));
	double m, b;
	for (unsigned int i = 0; i < size_slopes; i++){
		double *yfit, *xfit;
		yfit = malloc(sample * sizeof(double));
		xfit = malloc(sample * sizeof(double));
		for (unsigned int j = 0; j < sample; j++) {
			xfit[j] = x[i + j];
			yfit[j] = y[i + j];
		}
		linear_regression_double(xfit, yfit, sample, &m, &b);
		slopes[i] = m;
		free(yfit);
		free(xfit);
	}

	return slopes;
}