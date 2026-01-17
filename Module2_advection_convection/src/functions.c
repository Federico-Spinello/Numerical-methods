#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include <functions.h>


//condizione iniziale sinusoidale
void init_sin(double u[], double x[], int nx, double L){
    double pi = 4.0 * atan(1.0);

    for (int i = 0; i<=nx; i++){
        u[i] = sin(2.0 * pi * x[i] / L);
    }
}

//derivata seconda con schema centrato
void second_derivative(double u[], double ddu[], double dx, int nx) {
    int N = nx + 1;  // Numero totale di punti
    double inv_12dx2 = 1.0 / (12.0 * dx * dx);

    for (int i = 0; i <= nx; i++) {
        // CORREZIONE CRITICA: modulo deve essere N = nx+1, non nx!
        int ip2 = (i + 2) % N;
        int ip1 = (i + 1) % N;
        int im1 = (i - 1 + N) % N;
        int im2 = (i - 2 + N) % N;

        ddu[i] = (-u[ip2] + 16.0*u[ip1] - 30.0*u[i] + 16.0*u[im1] - u[im2]) * inv_12dx2;
    }
}
// upwind per c > 0 ovvero calcola la derivata decentrata prendendo il termine prima
void first_derivative_upwind(double u[], double du[], double dx, int nx) {
    int N = nx + 1;  // Numero totale di punti
    for (int i = 0; i <= nx; i++) {
        // CORREZIONE CRITICA: modulo deve essere N = nx+1
        int im1 = (i - 1 + N) % N;
        du[i] = (u[i] - u[im1]) / dx;
    }
}

 void first_derivative_4points(double u[], double du[], double dx, int nx ){
    int N = nx + 1;  // Numero totale di punti
    for (int i = 0; i <= nx; i++) {
        // CORREZIONE CRITICA: modulo deve essere N = nx+1
        int ip2 = (i + 2) % N;
        int ip1 = (i + 1) % N;
        int im1 = (i - 1 + N) % N;
        int im2 = (i - 2 + N) % N;
        du[i] = (u[im2] - 8.0 * u[im1] + 8.0 * u[ip1] - u[ip2]) / (12.0 * dx);
    }
 }

// RHS: F = nu * u_xx - u * u_x
// u and F are arrays of length N = nx+1
void RHS(double *u, double *F, int nx, double dx, double nu) {
    int N = nx + 1;
    // allocate temporaries on heap to avoid huge stack on large N
    double *du = malloc(N * sizeof(double));
    double *ddu = malloc(N * sizeof(double));
    if (!du || !ddu) {
        fprintf(stderr, "RHS: malloc failed\n");
        if (du) free(du);
        if (ddu) free(ddu);
        return;
    }

    // choose derivative; you can swap with first_derivative_4points if you want
    first_derivative_4points(u, du, dx, nx);
    second_derivative(u, ddu, dx, nx);

    for (int i = 0; i < N; i++) {
        F[i] = nu * ddu[i] - u[i] * du[i];
    }

    free(du);
    free(ddu);
}
void step(double *u, double *un, int N, double dt, double dx, double nu) {
    // N = nx+1 passato dal main
    int nx = N - 1;

    double *k1 = malloc(N * sizeof(double));
    double *k2 = malloc(N * sizeof(double));
    double *k3 = malloc(N * sizeof(double));
    double *k4 = malloc(N * sizeof(double));
    double *utmp = malloc(N * sizeof(double));

    if(!k1 || !k2 || !k3 || !k4 || !utmp){
        fprintf(stderr,"step (RK4): malloc failed\n");
        return;
    }

    // un ← u^n
    memcpy(un, u, N * sizeof(double));

    // ---- k1 = F(u^n)
    RHS(u, k1, nx, dx, nu);

    // ---- k2 = F(u^n + dt/2 k1)
    for(int i = 0; i < N; i++)
        utmp[i] = un[i] + 0.5 * dt * k1[i];
    RHS(utmp, k2, nx, dx, nu);

    // ---- k3 = F(u^n + dt/2 k2)
    for(int i = 0; i < N; i++)
        utmp[i] = un[i] + 0.5 * dt * k2[i];
    RHS(utmp, k3, nx, dx, nu);

    // ---- k4 = F(u^n + dt k3)
    for(int i = 0; i < N; i++)
        utmp[i] = un[i] + dt * k3[i];
    RHS(utmp, k4, nx, dx, nu);

    // ---- u^{n+1}
    for(int i = 0; i < N; i++) {
        u[i] = un[i] + (dt/6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }

    free(k1); free(k2); free(k3); free(k4); free(utmp);
}


void spectrum_fftw_real(const double *u, int M, double L, int step, const char *data_path) {
    // FFTW r2c output length = M/2 + 1 complex numbers
    fftw_plan plan;
    double *in = fftw_alloc_real(M);
    fftw_complex *out = fftw_alloc_complex(M/2 + 1);
    if (!in || !out) { fprintf(stderr, "FFTW alloc failed\n"); if(in) fftw_free(in); if(out) fftw_free(out); return; }

    // copy data
    for (int j = 0; j < M; j++) in[j] = u[j];

    plan = fftw_plan_dft_r2c_1d(M, in, out, FFTW_ESTIMATE);
    fftw_execute(plan);

    // normalization: FFTW computes unnormalized transform; we normalize by M.
    // compute power spectral density P(k) = |û_k|^2 / M^2 (so Parseval consistent)
    char filename[256];
    if (data_path[strlen(data_path)-1] == '/')
        sprintf(filename, "%sspectrum_%05d.dat", data_path, step);
    else
        sprintf(filename, "%s/spectrum_%05d.dat", data_path, step);

    FILE *fp = fopen(filename, "w");
    if (!fp) { fprintf(stderr, "Error opening %s\n", filename); fftw_destroy_plan(plan); fftw_free(in); fftw_free(out); return; }
    double pi = 4.0 * atan(1.0);
    double two_pi_over_L = 2.0 * M*pi / L;
    int K = M/2; // max positive freq index
    for (int k = 0; k <= K; k++) {
        double re = out[k][0];
        double im = out[k][1];
        double mag2 = re*re + im*im;
        double psd = mag2 / ((double)M * (double)M);
        double k_phys = two_pi_over_L * (double)k;

        // Errore sulla PSD: per una singola realizzazione FFT, ogni modo
        // segue una distribuzione chi-quadro con 2 DoF, quindi sigma_PSD ≈ PSD
        double sigma_psd = psd;

        fprintf(fp, "%g %g %g\n", k_phys, psd, sigma_psd);
    }

    fclose(fp);
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
}

void save_data_fft(double x[], double u[], int N, int step, const char *data_path) {
    char filename[200];
    char fftname[200];

    if (data_path[strlen(data_path)-1] == '/') {
        sprintf(filename, "%sdata_%d.dat", data_path, step);
        sprintf(fftname, "%sfft_%d.dat", data_path, step);
    } else {
        sprintf(filename, "%s/data_%d.dat", data_path, step);
        sprintf(fftname, "%s/fft_%d.dat", data_path, step);
    }

    // Salva dati spaziali
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        printf("Errore nell'apertura del file %s\n", filename);
        return;
    }
    for (int i = 0; i < N; i++) {
        fprintf(fp, "%g %g\n", x[i], u[i]);
    }
    fclose(fp);

    // FFT
    int Nfft = N-1;  // escludo l’ultimo punto periodico
    double *in = (double*) malloc(sizeof(double) * Nfft);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nfft/2+1));

    for (int i = 0; i < Nfft; i++) in[i] = u[i];

    fftw_plan p = fftw_plan_dft_r2c_1d(Nfft, in, out, FFTW_ESTIMATE);
    fftw_execute(p);

    double dx = x[1] - x[0];
    double L = dx * Nfft;

    // Salva FFT
    FILE *ffp = fopen(fftname, "w");
    if (!ffp) {
        printf("Errore nell'apertura del file %s\n", fftname);
        fftw_destroy_plan(p);
        free(in);
        fftw_free(out);
        return;
    }
    double pi = 4.0 * atan(1.0);
    
    for (int i = 0; i <= Nfft/2; i++) {
        double k = 2*pi*i / L;
        double psd = (out[i][0]*out[i][0] + out[i][1]*out[i][1]) / (Nfft*Nfft);

        // Errore sulla PSD: per una singola realizzazione FFT, ogni modo
        // segue una distribuzione chi-quadro con 2 DoF, quindi sigma_PSD ≈ PSD
        double sigma_psd = psd;

        fprintf(ffp, "%g %g %g\n", k, psd, sigma_psd);
    }

    fclose(ffp);
    fftw_destroy_plan(p);
    free(in);
    fftw_free(out);
}

