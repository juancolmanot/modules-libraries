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
    unsigned int n_steps;
    double t_stationary;
    double t_transient;
    double sigma;
    double rho;
    double beta;
    double h;
    double atol;
    double rtol;
} Parameters5;

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
