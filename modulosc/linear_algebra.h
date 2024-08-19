#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

long double **jacobian_numerical(
    void (*func_n)(
        void (*)(
            long double *,
            long double *,
            void *
        ),
        long double *,
        long double *,
        int ,
        void *
    ),
    void (*func)(
            long double *,
            long double *,
            void *
        ),
    long double *x,
    void *params,
    int n,
    int n_map,
    long double epsilon
);

void qr_decomposition(
    long double **A,
    long double **Q,
    long double **R,
    int rows,
    int cols    
);

// Get min value from array
long double la_min(
	long double* x,
	unsigned int length
);

double la_min_d(
    double* x,
    unsigned int length
);

// Get max value from array
long double la_max(
	long double* x,
	unsigned int length
);

double la_max_d(
    double* x,
    unsigned int length
);

// Rotate two dimensional state an angle theta.
void rotate_state(
    long double *x,
    long double *xrot,
    long double theta
);

// Rotate two vectors of two dimensional states an angle theta.
void rotate_vectors(
    long double *x,
    long double *y,
    long double *xr,
    long double *yr,
    long double theta,
    unsigned int size
);

// Find angle that aligns x data with xaxis.
void align_axes(
    long double *x,
    long double *y,
    long double *xal,
    long double *yal,
    unsigned int size,
    long double theta0,
    long double tol,
    unsigned int max_iter,
    long double *err,
    unsigned int *errsize
);

// Compute distance from 4 dimensional state to hyper-surface.
long double distance(
	long double *xn,
	long double *xn1
);

#endif