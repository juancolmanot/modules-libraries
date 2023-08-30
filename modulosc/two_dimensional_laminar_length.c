#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "mt64.h"
#include "linear_algebra.h"

double laminar_length_avg(unsigned int N, double c, double epsilon, double *params) {

    /* seed for initial state generator */
    unsigned long long seed = (unsigned long long) rand();
    /* count of processes with no zero laminar regions */
    unsigned int valid_count = 0;
    /* Number of map to evolve */
    unsigned int n_map = 14;
    /* Number of experiments */
    unsigned int n_tries = 50;
    
    /* Parameters of the map */
    double alpha_c = params[0];
    double beta = params[1];
    double alpha_i = alpha_c - epsilon;

    /* statistical variables */
    unsigned int reinjection_counter;
    unsigned int ejection_counter;
    unsigned int iterations;
    unsigned int start_laminar;

    /* laminar length variables */
    double l_avg;

    double l_exp = 0;

    /* State variables */
    double* xn = calloc(2 , sizeof(double));
    double* xn1 = calloc(2, sizeof(double));
    double* xerr = calloc(2, sizeof(double));
    double* x1err = calloc(2, sizeof(double));

    init_genrand64(seed);

    for (unsigned int i = 0; i < n_tries; i++) {

        reinjection_counter = 0;
        ejection_counter = 0;
        iterations = 0;
        start_laminar = 0;

        l_avg = 0;

        xn[0] = genrand64_real3();
        xn[1] = genrand64_real3();
        xn1[0] = xn1[1] = 0;
        xerr[0] = xerr[1] = x1err[0] = x1err[1] = 0;

        for (unsigned int j = 0; j < N; j++) {

            map_2d_n(xn1, xn, n_map, alpha_i, beta);
            rel_err(xn, xn1, x1err);

            if (xerr[0] > c && x1err[0] < c && iterations > 2) {
                reinjection_counter++;
                start_laminar = iterations;
            }
            else if (xerr[0] < c && x1err[0] > c && iterations > 2 && reinjection_counter > 0) {
                ejection_counter++;
                l_avg += (double) (iterations - start_laminar);
            }
            iterations++;

            xn[0] = xn1[0];
            xn[1] = xn1[1];

            xerr[0] = x1err[0];
            xerr[1] = x1err[1];
        }
        
        if (ejection_counter > 0) {
            valid_count++;
            l_avg = l_avg / (double) ejection_counter;
            l_exp += l_avg;
        }
    }
    l_exp = l_exp / (double)valid_count;

    return l_exp;
}