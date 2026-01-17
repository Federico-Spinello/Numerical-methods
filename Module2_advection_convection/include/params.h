#ifndef PARAMS_H
#define PARAMS_H

// Struttura per parametri globali
typedef struct {
    int nx;
    double L;
    double CFL_diff;
    double CFL_adv;
    double safety;
} GlobalParams;

// Struttura per parametri advection-convection
typedef struct {
    int nt;
    int print_step;
    double nu_factor;
} AdvConvParams;

// Struttura per parametri shock analysis
typedef struct {
    int nt;
    double nu_min_factor;
    double nu_max_factor;
    int n_nu;
} ShockAnalysisParams;

// Funzioni per leggere i parametri
int read_global_params(const char *filename, GlobalParams *params);
int read_advconv_params(const char *filename, AdvConvParams *params);
int read_shock_params(const char *filename, ShockAnalysisParams *params);

#endif