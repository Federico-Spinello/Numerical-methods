#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "functions.h"
#include "params.h"

int sim_shock_analysis(void) {
    // Leggi parametri
    GlobalParams gp;
    ShockAnalysisParams sap;
    
    if (!read_global_params("params.txt", &gp)) return 1;
    if (!read_shock_params("params.txt", &sap)) return 1;
    
    // Calcola parametri derivati
    int N = gp.nx + 1;
    double dx = gp.L / gp.nx;
    double pi = 4.0 * atan(1.0);
    
    double nu_min = dx * dx * pi * sap.nu_min_factor / 2.0;
    double nu_max = dx * dx * pi * sap.nu_max_factor / 2.0;
    
    char *output_file = "./data/shock_values.txt";
    
    // Indice corrispondente a x = 0.5
    int idx_shock = (int)(0.5 / dx);
    if (idx_shock >= N) idx_shock = N - 1;
    
    // Allocazione
    double *x = malloc(N * sizeof(double));
    double *u = malloc(N * sizeof(double));
    double *un = malloc(N * sizeof(double));
    double *du = malloc(N * sizeof(double));
    
    if (!x || !u || !un || !du) {
        fprintf(stderr, "malloc failed in sim_shock_analysis\n");
        return 1;
    }
    
    for (int i = 0; i < N; i++) {
        x[i] = i * dx;
    }
    
    printf("\n=== ANALISI SHOCK VS VISCOSITÀ ===\n");
    printf("Parametri da params.txt:\n");
    printf("  nx = %d, nt = %d\n", gp.nx, sap.nt);
    printf("  nu_min_factor = %.2f -> nu_min = %.6e\n", sap.nu_min_factor, nu_min);
    printf("  nu_max_factor = %.2f -> nu_max = %.6e\n", sap.nu_max_factor, nu_max);
    printf("  n_nu = %d\n", sap.n_nu);
    printf("  Scala: logaritmica\n");
    printf("  Posizione shock: x = %.2f (indice %d)\n", x[idx_shock], idx_shock);
    printf("  Condizione iniziale: sin(2πx/L)\n");
    printf("  Output: %s\n\n", output_file);
    
    // Apri file in append mode
    FILE *fp = fopen(output_file, "a");
    if (!fp) {
        fprintf(stderr, "Errore apertura file %s\n", output_file);
        free(x); free(u); free(un); free(du);
        return 1;
    }
    
    // Scrivi header solo se file vuoto
    fseek(fp, 0, SEEK_END);
    if (ftell(fp) == 0) {
        fprintf(fp, "# Viscosity  Min_Derivative_at_x0.5  TimeStep_Number\n");
    }
    
    // Loop sui valori di nu (scala logaritmica)
    for (int i_nu = 0; i_nu < sap.n_nu; i_nu++) {
        double log_nu_min = log10(nu_min);
        double log_nu_max = log10(nu_max);
        double log_nu = log_nu_min + (log_nu_max - log_nu_min) * i_nu / (sap.n_nu - 1);
        double nu = pow(10.0, log_nu);
        
        printf("Simulazione %d/%d: nu = %.6e", i_nu + 1, sap.n_nu, nu);
        fflush(stdout);
        
        // Condizione iniziale
        init_sin(u, x, gp.nx, gp.L);
        
        double dt_diff = gp.CFL_diff * dx * dx / nu;
        
        // Variabili per tracciare il minimo NEL TEMPO in x=0.5
        double min_du_at_x05 = 0.0;  // Minimo nel tempo
        
        // Time loop con controllo stabilità
        int diverged = 0;
        for (int tstep = 0; tstep < sap.nt; tstep++) {
            double umax = 0.0;
            for (int i = 0; i < N; i++) {
                umax = fmax(umax, fabs(u[i]));
            }
            if (umax < 1e-12) umax = 1e-12;

            // Controllo divergenza
            if (!isfinite(umax) || umax > 1e6) {
                printf(" [DIVERGENZA a step %d, skipping]", tstep);
                diverged = 1;
                break;
            }

            double dt_adv = gp.CFL_adv * dx / umax;
            double dt = gp.safety * fmin(dt_adv, dt_diff);

            step(u, un, N, dt, dx, nu);

            // Calcola derivata prima
            first_derivative_4points(u, du, dx, gp.nx);

            // Valore della derivata in x = 0.5
            double du_at_shock = du[idx_shock];

            // Aggiorna se troviamo un minimo più basso NEL TEMPO
            if (du_at_shock < min_du_at_x05) {
                min_du_at_x05 = du_at_shock;
            }
        }

        if (diverged) {
            printf(" -> SKIPPED (diverged)\n");
            continue;  // Salta questo valore di nu
        }
        
        printf(" -> min(du/dx)|_{x=0.5} = %.6e \n", min_du_at_x05);
        
        // Salva risultato (append)
        fprintf(fp, "%.10e  %.10e \n", nu, min_du_at_x05);
        fflush(fp);
    }
    
    fclose(fp);
    
    printf("\nAnalisi completata! Risultati in: %s\n\n", output_file);
    
    free(x);
    free(u);
    free(un);
    free(du);
    
    return 0;
}