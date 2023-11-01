#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

gsl_matrix *file_data(char filename[1024], int *shape){

    FILE *file;

    file = fopen(filename, "r");

    if (file == NULL) {
        printf("Error opening file.\n");
        return 1;
    }

    gsl_matrix *data = gsl_matrix_alloc(shape[0], shape[1]);
    float a, b, c;

    for (unsigned int i = 0; i < shape[0]; i++){
        fscanf(file, "%f %f %f", &a, &b, &c);
        printf("%f %f %f\n", a, b, c);
        gsl_matrix_set(data, i, 0, (double)a);
        gsl_matrix_set(data, i, 1, (double)b);
    }

    return data;
}

gsl_vector *file_data_vector(char filename[1024], int length) {

    FILE *file;

    file = fopen(filename, "r");

    if (file == NULL) {
        printf("Error opening file.\n");
        return 1;
    }

    gsl_vector *data = gsl_vector_alloc(length);
    float a;

    for (unsigned int i = 0; i < length; i++){
        fscanf(file, "%f", &a);
        gsl_vector_set(data, i, (double)a);
    }

    return data;
}
