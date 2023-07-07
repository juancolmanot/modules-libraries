#ifndef INTERMITTENCY_H
#define INTERMITTENCY_H

int sys1(double t, double x[], double f[], void *params);

int jac1(double t, double x[], double *dfdx, double dfdt[], void *params);

int sys2(double t, double x[], double f[], void *params);

int jac2(double t, double x[], double *dfdx, double dfdt[], void *params);

int sys3(double t, double x[], double f[], void *params);

int jac3(double t, double x[], double *dfdx, double dfdt[], void *params);

int sys4(double t, double x[], double f[], void *params);

int jac4(double t, double x[], double *dfdx, double dfdt[], void *params);

int sys5(double t, double x[], double f[], void *params);

int jac5(double t, double x[], double *dfdx, double dfdt[], void *params);

#endif