#ifndef TYPEII_COUPLED
#define TYPEII_COUPLED

typedef struct {
    double mu;
    double gamma;
    double b;
} Parameters;

double LMBKx(double x, void *params);
double LMBKy(double y, void *params);

#endif