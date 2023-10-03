#ifndef READ_FILE_H
#define READ_FILE_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

gsl_matrix *file_data(char filename[1024], int *shape);

gsl_vector *file_data_vector(char filename[1024], int length);

#endif