#ifndef MULTIPURPOSE_ROUTINES_H
#define MULTIPURPOSE_ROUTINES_H

/*
This function takes an array of data x, and an upper
and lower limits (x_min, x_max), and returns and array
of all the x_i such that x_min <= x_i <= x_max, and a 
pointer to its size n_xi.
*/
void get_region_array(
    long double *x_data,
    long double x_min,
    long double x_max,
    unsigned int size_data,
    long double **xi,
    unsigned int *n_xi
);

unsigned int get_value_loc(
    long double *xdata,
    long double x,
    unsigned int n
);

#endif