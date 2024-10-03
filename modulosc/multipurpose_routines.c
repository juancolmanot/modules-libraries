#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
This function takes an array of sorted data x, and an upper
and lower limits (x_min, x_max), and returns and array
of all the x_i such that x_min <= x_i <= x_max.
*/
void get_region_array(
    long double *x_data,
    long double x_min,
    long double x_max,
    unsigned int size_data,
    long double **xi,
    unsigned int *n_xi
)
{   
    unsigned int count = 0;
    unsigned int i_start = 0, i_end = 0;


    for (unsigned int i = 0; i < size_data; i++){
        if (x_data[i] >= x_min && x_data[i] < x_max){
            if (count == 0) {
                i_start = i;
            }
            i_end = i;
            count++;
        }
    }

    *n_xi = count;

    long double *xireall = realloc(*xi, count * sizeof(long double));
    if (xireall == NULL){
        fprintf(stderr, "Error reallocating memory.\n");
    }
    
    *xi = xireall;
    
    for (unsigned int i = i_start; i < i_end + 1; i++){
        (*xi)[i - i_start] = x_data[i];
    }

}

void get_region_array_scrambled(
    long double *x_data,
    long double x_min,
    long double x_max,
    unsigned int size_data,
    long double **xi,
    unsigned int *n_xi
)
{   
    unsigned int count = 0;

    for (unsigned int i = 0; i < size_data; i++){
        if (x_data[i] >= x_min && x_data[i] < x_max){
            if (count == 0) {
            }
            count++;
        }
    }

    *n_xi = count;
    long double *xireall = realloc(*xi, count * sizeof(long double));
    if (xireall == NULL){
        fprintf(stderr, "Error reallocating memory.\n");
    }
    
    *xi = xireall;
    
    count = 0;
    for (unsigned int i = 0; i < size_data; i++){
        if (x_data[i] >= x_min && x_data[i] < x_max){
            (*xi)[count] = x_data[i];
            count++;
        }
    }

}

void get_region_indeces(
    long double *x_data,
    long double x_min,
    long double x_max,
    unsigned int size_data,
    unsigned int **indeces,
    unsigned int *n_idx
){
    unsigned int count = 0;

    for (unsigned int i = 0; i < size_data; i++){
        if (x_data[i] >= x_min && x_data[i] < x_max){
            if (count == 0) {
            }
            count++;
        }
    }

    *n_idx = count;
    unsigned int *indecesreall = realloc(*indeces, count * sizeof(unsigned int));
    if (indecesreall == NULL){
        fprintf(stderr, "Error reallocating memory.\n");
    }
    
    *indeces = indecesreall;
    
    count = 0;
    for (unsigned int i = 0; i < size_data; i++){
        if (x_data[i] >= x_min && x_data[i] < x_max){
            (*indeces)[count] = i;
            count++;
        }
    }
}

unsigned int get_value_loc(
    long double *xdata,
    long double x,
    unsigned int n
)
{
    unsigned int index = 0;
    for (unsigned int i = 0; i < n + 1; i++){
        if (i < n) {
            if (xdata[i] == x) {
                index = i;
                break;
            }    
        }
        else {
            index = (unsigned int)(-10);
            return index;
        }
    }

    return index;
}