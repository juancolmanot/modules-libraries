#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H


double* linspace(double x0, double xn, unsigned int n);
double* logspace(double exp0, double expn, unsigned int n);
double* linspace_discrete(double x0, double xn, unsigned int n);
double* rotation_operator(double* x, double theta);
double min(double* x, unsigned int length);
double max(double* x, unsigned int length);
double distance(double* x);
void rel_err(double* x, double* x1, double* xerr);
double rel_err_scalar(double x, double x1);

#endif