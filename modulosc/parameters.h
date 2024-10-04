#ifndef PARAMETERS1_H
#define PARAMETERS1_H

typedef struct {
    long unsigned int N;
    unsigned int ntarget;
    double threshold;
    unsigned int transient;
    unsigned int n_map;
    long double p_alpha;
    double p_beta;
} Parameters1;

typedef struct {
    long unsigned int N;
    unsigned int transient;
    unsigned int n_map;
    long double p_alpha;
    double p_beta;
} Parameters2;

typedef struct {
    unsigned int n_steps;
    int n;
    int n_map;
    unsigned int transient;
    long double p_alpha;
    long double p_beta;
    long double epsilon;
} Parameters3;

typedef struct {
    unsigned int nevol;
    int n;
    int n_map;
    unsigned int transient;
    long double alpha0;
    long double alphan;
    long double p_beta;
    unsigned int nalphas;
    long double theta;
} Parameters4;

typedef struct {
    long unsigned int max_steps;
    double t_stationary;
    double t_transient;
    double sigma;
    double rho;
    double beta;
    double h;
    double atol;
    double rtol;
} Parameters5;

typedef struct {
    long unsigned int max_steps;
    double t_stationary;
    double t_transient;
    double sigma;
    double rho;
    double beta;
    double h;
    double atol;
    double rtol;
    double xmin;
    double xmax;
    unsigned int rtarget;
} Parameters6;

typedef struct {
    int n_regions;
    int n_centers;
    unsigned int npoints;
    int current_zone;
    double *xmins;
    double *xmaxs;
    double *xc;
    unsigned int *region_bins;
} Parameters7;

typedef struct {
    int n_regions;
    int n_centers;
    double wi;
    double *xmins;
    double *xmaxs;
    double *xc;
} Parameters7a;

typedef struct {
    int n_regions;
    int n_centers;
    unsigned int npoints;
    int current_zone;
    double *xmins_fit;
    double *xmaxs_fit;
    double *xmins_domain;
    double *xmaxs_domain;
    double *xc;
    unsigned int *region_bins;
} Parameters8;

/* Type of parameters struct for separating
regions of map x(n - 1) vs x(n) and isolating
effects on M and RPD functions.*/
typedef struct {
    int n_regions;
    double *xmins_regions;
    double *xmaxs_regions;
    int *bins;
} Parameters9;

typedef struct {
    int n_regions;
    double *xmins_regions;
    double *xmaxs_regions;
    int *bins;
    long double *wi;
} Parameters10;

int handler1(
    void *user,
    const char* section,
    const char* name,
    const char* value
);

int handler2(
    void *user,
    const char* section,
    const char* name,
    const char* value
);

int handler3(
    void *user,
    const char* section,
    const char* name,
    const char* value
);

int handler4(
    void *user,
    const char* section,
    const char* name,
    const char* value
);

int handler5(
    void *user,
    const char* section,
    const char* name,
    const char* value
);

int handler6(
    void *user,
    const char* section,
    const char* name,
    const char* value
);

int handler7(
    void *user,
    const char* section,
    const char* name,
    const char* value
);

int handler7a(
    void *user,
    const char* section,
    const char* name,
    const char* value
);

int handler8(
    void *user,
    const char* section,
    const char* name,
    const char* value
);

/* Deals with loading parameters of type 9 */
int handler9(
    void *user,
    const char* section,
    const char* name,
    const char* value
);

int handler10(
    void *user,
    const char* section,
    const char* name,
    const char* value
);

void load_parameters_from_file(
    const char* filename,
    void* params,
    int handler(
        void * ,
        const char* ,
        const char* ,
        const char* 
    )
);

#endif
