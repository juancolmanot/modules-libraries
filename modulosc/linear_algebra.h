#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H


double* linspace(double x0, double xn, unsigned int n);

double* linspace_discrete(double x0, double xn, unsigned int n);

double min(double* x);

double max(double* x);

double distance(double* x);

void rel_err(double* x, double* x1, double* xerr);

double rel_err_scalar(double x, double x1);

#endif