#include <stdio.h>
#include <stdlib.h>

int *get_n_lines(char filename[1024]){

    FILE *file;
    int *shape = calloc(2, sizeof(double));


    file = fopen(filename, "r");

    if (file == NULL) {
        printf("Error opening file.\n");
        return 1;
    }

    int rows = 1;
    int cols = 0;
    char c;
    char line[1024];

    if (fgets(line, sizeof(line), file) != NULL){
        char *token = strtok(line, " \t");
        while (token != NULL) {
            cols++;
            token = strtok(NULL, " \t");
        }
    }
    else {
        printf("File is empty.\n");
        fclose(file);
        return 1;
    }

    do {
        c = fgetc(file);
        if (c == '\n') rows++;

    } while (c != EOF);

    fclose(file);

    shape[0] = rows;
    shape[1] = cols;

    return shape;
}