#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H
#include <gsl/gsl_vector.h>

gsl_vector *gsl_linspace(double x0, double xN, double N);
double* linspace(double x0, double xn, unsigned int n);
double* logspace(double exp0, double expn, unsigned int n);
double* linspace_discrete(double x0, double xn, unsigned int n);
void rotation_operator(double *xrot, double* x, double theta);
void antirotation_operator(double *xrot, double* x, double theta);
double min(double* x, unsigned int length);
double max(double* x, unsigned int length);
double distance(double* x);
void rel_err(double* x, double* x1, double* xerr);
double rel_err_scalar(double x, double x1);
void distance_line(double *x, double *x1, double *xerr);

#endif