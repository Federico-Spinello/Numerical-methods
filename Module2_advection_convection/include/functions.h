#ifndef FUNCTIONS_H
#define FUNCTIONS_H


//condizione iniziale
void init_sin(double u[], double x[], int nx, double L);
//derivata seconda con schema centrato
void second_derivative(double u[], double ddu[], double dx, int nx);

// upwind per c > 0 ovvero calcola la derivata decentrata prendendo il termine prima
void first_derivative_upwind(double u[], double du[], double dx, int nx);
void first_derivative_4points(double u[], double du[], double dx, int nx );

void RHS(double *u, double *F, int nx, double dx, double nu);
void step(double *u, double *un, int N, double dt, double dx, double nu);

void save_data_fft(double x[], double u[], int N, int step, const char *data_path);

#endif
