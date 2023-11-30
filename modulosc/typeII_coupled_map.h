#ifndef TYPEII_COUPLED
#define TYPEII_COUPLED

typedef struct {
    double mu;
    double gamma;
    double b;
    double eps;
} Parameters;

double *LMBK(double *x, void *params);

#endif