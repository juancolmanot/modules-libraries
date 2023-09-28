#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>

gsl_matrix *file_data(char filename[1024], int *shape){

    FILE *file;

    file = fopen(filename, "r");

    if (file == NULL) {
        printf("Error opening file.\n");
        return 1;
    }

    gsl_matrix *data = gsl_matrix_alloc(shape[0], shape[1]);
    float a, b;

    for (unsigned int i = 0; i < shape[0]; i++){
        fscanf(file, "%f %f", &a, &b);
        gsl_matrix_set(data, i, 0, (double)a);
        gsl_matrix_set(data, i, 1, (double)b);
    }

    return data;
}