#ifndef READ_FILE_H
#define READ_FILE_H

#include "gsl/gsl_matrix.h"

gsl_matrix *file_data(char filename[1024], int *shape);

#endif