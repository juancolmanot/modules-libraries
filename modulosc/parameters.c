#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parameters.h"
#include "ini.h"

int handler1(
    void *user,
    const char* section,
    const char* name,
    const char* value
){
    Parameters1 *p = (Parameters1*)user;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("general", "N")) {
 		p->N = strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "ntarget")) {
    	p->ntarget = (unsigned int)strtoul(value, NULL, 10);
    }
	else if (MATCH("general", "threshold")) {
        p->threshold = strtod(value, NULL);
    } else if (MATCH("general", "transient")) {
        p->transient = (unsigned int)strtoul(value, NULL, 10);
    } else if (MATCH("general", "n_map")) {
        p->n_map = (unsigned int)strtoul(value, NULL, 10);
    } else if (MATCH("general", "p_alpha")) {
        p->p_alpha = strtold(value, NULL);
    } else if (MATCH("general", "p_beta")) {
        p->p_beta = strtod(value, NULL);
    } else {
        return 0;  // unknown section/name, error
    }
    return 1;
}

int handler2(
    void *user,
    const char* section,
    const char* name,
    const char* value
){
    Parameters2 *p = (Parameters2*)user;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("general", "N")) {
        p->N = strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "transient")) {
        p->transient = (unsigned int)strtoul(value, NULL, 10);
    } else if (MATCH("general", "n_map")) {
        p->n_map = (unsigned int)strtoul(value, NULL, 10);
    } else if (MATCH("general", "p_alpha")) {
        p->p_alpha = strtold(value, NULL);
    } else if (MATCH("general", "p_beta")) {
        p->p_beta = strtod(value, NULL);
    } else {
        return 0;  // unknown section/name, error
    }
    return 1;
}

int handler3(
    void *user,
    const char* section,
    const char* name,
    const char* value
){
    Parameters3 *p = (Parameters3*)user;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("general", "n_steps")) {
        p->n_steps = (unsigned int)strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "n")) {
        p->n = (int)strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "n_map")) {
        p->n_map = (int)strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "transient")){
        p->transient = (unsigned int)strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "p_alpha")) {
        p->p_alpha = strtold(value, NULL);
    }
    else if (MATCH("general", "p_beta")) {
        p->p_beta = strtod(value, NULL);
    }
    else if (MATCH("general", "epsilon")) {
        p->epsilon = strtod(value, NULL);
    } else {
        return 0;  // unknown section/name, error
    }
    return 1;
}

int handler4(
    void *user,
    const char* section,
    const char* name,
    const char* value
){
    Parameters4 *p = (Parameters4*)user;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("general", "nevol")) {
        p->nevol = (unsigned int)strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "n")) {
        p->n = (int)strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "n_map")) {
        p->n_map = (int)strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "transient")){
        p->transient = (unsigned int)strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "alpha0")) {
        p->alpha0 = strtold(value, NULL);
    }
    else if (MATCH("general", "alphan")) {
        p->alphan = strtold(value, NULL);
    }    
    else if (MATCH("general", "p_beta")) {
        p->p_beta = strtod(value, NULL);
    }
    else if (MATCH("general", "nalphas")) {
        p->nalphas = (unsigned int)strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "theta")) {
        p->theta = strtold(value, NULL);
    } else {
        return 0;  // unknown section/name, error
    }
    return 1;
}

int handler5(
    void *user,
    const char* section,
    const char* name,
    const char* value
){
    Parameters5 *p = (Parameters5*)user;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("general", "max_steps")) {
        p->max_steps = strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "t_stationary")) {
        p->t_stationary = strtod(value, NULL);
    }
    else if (MATCH("general", "t_transient")) {
        p->t_transient = strtod(value, NULL);
    }
    else if (MATCH("general", "sigma")) {
        p->sigma = strtod(value, NULL);
    }
    else if (MATCH("general", "rho")) {
        p->rho = strtod(value, NULL);
    }    
    else if (MATCH("general", "beta")) {
        p->beta = strtod(value, NULL);
    }
    else if (MATCH("general", "h")) {
        p->h = strtod(value, NULL);
    }
    else if (MATCH("general", "atol")) {
        p->atol = strtod(value, NULL);
    }
    else if (MATCH("general", "rtol")) {
        p->rtol = strtod(value, NULL);
    } else {
        return 0;  // unknown section/name, error
    }
    return 1;
}

int handler6(
    void *user,
    const char* section,
    const char* name,
    const char* value
){
    Parameters6 *p = (Parameters6*)user;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("general", "max_steps")) {
        p->max_steps = strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "t_stationary")) {
        p->t_stationary = strtod(value, NULL);
    }
    else if (MATCH("general", "t_transient")) {
        p->t_transient = strtod(value, NULL);
    }
    else if (MATCH("general", "sigma")) {
        p->sigma = strtod(value, NULL);
    }
    else if (MATCH("general", "rho")) {
        p->rho = strtod(value, NULL);
    }    
    else if (MATCH("general", "beta")) {
        p->beta = strtod(value, NULL);
    }
    else if (MATCH("general", "h")) {
        p->h = strtod(value, NULL);
    }
    else if (MATCH("general", "atol")) {
        p->atol = strtod(value, NULL);
    }
    else if (MATCH("general", "rtol")) {
        p->rtol = strtod(value, NULL);
    }
    else if (MATCH("general", "xmin")) {
        p->xmin = strtod(value, NULL);
    }
    else if (MATCH("general", "xmax")) {
        p->xmax = strtod(value, NULL);
    }
    else if (MATCH("general", "rtarget")) {
        p->rtarget = (unsigned int)strtoul(value, NULL, 10);
    }
    else {
        return 0;  // unknown section/name, error
    }
    return 1;
}

int handler7(
    void *user,
    const char* section,
    const char* name,
    const char* value
){
    Parameters7 *p = (Parameters7*)user;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("general", "n_regions")) {
        p->n_regions = (int)strtoul(value, NULL, 10);
        p->xmins = (double*)malloc((size_t)p->n_regions * sizeof(double));
        p->xmaxs = (double*)malloc((size_t)p->n_regions * sizeof(double));
        p->region_bins = (unsigned int*)malloc((size_t)p->n_regions * sizeof(unsigned int));
        if (p->xmins == NULL || p->xmaxs == NULL || p->region_bins == NULL) {
            perror("xmins, xmaxs or region_bins memory allocation failed\n");
            exit(1);
        }
    }
    else if (MATCH("general", "n_centers")) {
        p->n_centers = (int)strtoul(value, NULL, 10);
        p->xc = (double*)malloc((size_t)p->n_centers * sizeof(double));
        if (p->xc == NULL) {
            perror("xc memory allocation failed\n");
            exit(1);
        }
    }
    else if (MATCH("general", "current_zone")) {
        p->current_zone = (int)strtoul(value, NULL, 10);
    }
    else {
        char region_name[64];
        char bins_name[64];
        char center_name[64];
        for (int i = 0; i < p->n_regions; i++) {
            snprintf(region_name, sizeof(region_name), "xmin%u", i + 1);
            if (MATCH("regions", region_name)) {
                p->xmins[i] = strtod(value, NULL);
                return 1;
            }
            snprintf(region_name, sizeof(region_name), "xmax%u", i + 1);
            if (MATCH("regions", region_name)) {
                p->xmaxs[i] = strtod(value, NULL);
                return 1;
            }
            snprintf(bins_name, sizeof(bins_name), "bins%u", i + 1);
            if (MATCH("bins", bins_name)) {
                p->region_bins[i] = (unsigned int)strtoul(value, NULL, 10);
                return 1;
            }
        }

        for (int i = 0; i < p->n_centers; i++) {
            snprintf(center_name, sizeof(center_name), "xc%u", i + 1);
            if (MATCH("centers", center_name)) {
                p->xc[i] = strtod(value, NULL);
                return 1;
            }
        }
        return 0; // unknown section/name, error
    }
    return 1;
}

int handler8(
    void *user,
    const char* section,
    const char* name,
    const char* value
){
    Parameters8 *p = (Parameters8*)user;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("general", "n_regions")) {
        p->n_regions = (int)strtoul(value, NULL, 10);
        p->xmins_fit = (double*)malloc((size_t)p->n_regions * sizeof(double));
        p->xmaxs_fit = (double*)malloc((size_t)p->n_regions * sizeof(double));
        p->xmins_domain = (double*)malloc((size_t)p->n_regions * sizeof(double));
        p->xmaxs_domain = (double*)malloc((size_t)p->n_regions * sizeof(double));
        p->region_bins = (unsigned int*)malloc((size_t)p->n_regions * sizeof(unsigned int));
        if (p->xmins_fit == NULL || p->xmaxs_fit == NULL || p->region_bins == NULL) {
            perror("xmins, xmaxs or region_bins memory allocation failed\n");
            exit(1);
        }
        if (p->xmins_domain == NULL || p->xmaxs_domain == NULL) {
            perror("xmins or xmaxs of domain memory allocation failed\n");
        }
    }
    else if (MATCH("general", "n_centers")) {
        p->n_centers = (int)strtoul(value, NULL, 10);
        p->xc = (double*)malloc((size_t)p->n_centers * sizeof(double));
        if (p->xc == NULL) {
            perror("xc memory allocation failed\n");
            exit(1);
        }
    }
    else if (MATCH("general", "current_zone")) {
        p->current_zone = (int)strtoul(value, NULL, 10);
    }
    else {
        char region_fit_name[64];
        char region_domain_name[64];
        char bins_name[64];
        char center_name[64];
        for (int i = 0; i < p->n_regions; i++) {
            snprintf(region_fit_name, sizeof(region_fit_name), "xmin_fit%u", i + 1);
            if (MATCH("regions_fit", region_fit_name)) {
                p->xmins_fit[i] = strtod(value, NULL);
                return 1;
            }
            snprintf(region_fit_name, sizeof(region_fit_name), "xmax_fit%u", i + 1);
            if (MATCH("regions_fit", region_fit_name)) {
                p->xmaxs_fit[i] = strtod(value, NULL);
                return 1;
            }
            snprintf(region_domain_name, sizeof(region_domain_name), "xmin_domain%u", i + 1);
            if (MATCH("regions_domain", region_domain_name)) {
                p->xmins_domain[i] = strtod(value, NULL);
                return 1;
            }
            snprintf(region_domain_name, sizeof(region_domain_name), "xmax_domain%u", i + 1);
            if (MATCH("regions_domain", region_domain_name)) {
                p->xmaxs_domain[i] = strtod(value, NULL);
                return 1;
            }
            snprintf(bins_name, sizeof(bins_name), "bins%u", i + 1);
            if (MATCH("bins", bins_name)) {
                p->region_bins[i] = (unsigned int)strtoul(value, NULL, 10);
                return 1;
            }
        }

        for (int i = 0; i < p->n_centers; i++) {
            snprintf(center_name, sizeof(center_name), "xc%u", i + 1);
            if (MATCH("centers", center_name)) {
                p->xc[i] = strtod(value, NULL);
                return 1;
            }
        }
        return 0; // unknown section/name, error
    }
    return 1;
}

void load_parameters_from_file(
    const char* filename,
    void *params,
    int handler(
        void * ,
        const char* ,
        const char* ,
        const char* 
    )
){
	if (ini_parse(filename, handler, params) < 0){
		printf("Can't load '%s'\n", filename);
		exit(1);
	}
}
