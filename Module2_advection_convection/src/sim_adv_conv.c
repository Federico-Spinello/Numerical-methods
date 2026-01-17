#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "functions.h"
#include "params.h"

int sim_adv_conv(void) {
    // Leggi parametri
    GlobalParams gp;
    AdvConvParams acp;
    
    if (!read_global_params("params.txt", &gp)) return 1;
    if (!read_advconv_params("params.txt", &acp)) return 1;
    
    // Calcola parametri derivati
    int N = gp.nx + 1;
    double dx = gp.L / gp.nx;
    double pi = 4.0 * atan(1.0);
    double nu = dx * dx * pi * acp.nu_factor / 2.0;
    double dt_diff = gp.CFL_diff * dx * dx / nu;
    
    // Allocazione
    double *x = malloc(N * sizeof(double));
    double *u = malloc(N * sizeof(double));
    double *un = malloc(N * sizeof(double));
    
    if (!x || !u || !un) {
        fprintf(stderr, "malloc failed in sim_adv_conv\n");
        return 1;
    }
    
    char *data_path = "./data/";
    
    // Griglia
    for (int i = 0; i < N; i++) {
        x[i] = i * dx;
    }
    
    // Condizione iniziale: sinusoide
    init_sin(u, x, gp.nx, gp.L);
    
    // Salva condizione iniziale
    save_data_fft(x, u, N, 0, data_path);
    
    printf("\n=== SIMULAZIONE ADVECTION-CONVECTION ===\n");
    printf("Parametri da params.txt:\n");
    printf("  nx = %d, nt = %d\n", gp.nx, acp.nt);
    printf("  L = %.2f, dx = %.6f\n", gp.L, dx);
    printf("  nu_factor = %.2f -> nu = %.6e\n", acp.nu_factor, nu);
    printf("  CFL_diff = %.2f, CFL_adv = %.2f, safety = %.2f\n", 
           gp.CFL_diff, gp.CFL_adv, gp.safety);
    printf("  dt_diff = %.6e\n", dt_diff);
    printf("  Condizione iniziale: sin(2πx/L)\n");
    printf("\nInizio time loop...\n");
    
    // Loop temporale con controlli di stabilità
    int divergence_detected = 0;
    for (int tstep = 0; tstep < acp.nt; tstep++) {
        double umax = 0.0;
        double energy = 0.0;  // Controllo energia
        for (int i = 0; i < N; i++) {
            umax = fmax(umax, fabs(u[i]));
            energy += u[i] * u[i];
        }
        energy = sqrt(energy / N);  // Norma L2

        if (umax < 1e-12) umax = 1e-12;

        // Controllo divergenza
        if (!isfinite(umax) || umax > 1e6 || energy > 1e6) {
            fprintf(stderr, "\n*** DIVERGENZA RILEVATA a step %d ***\n", tstep);
            fprintf(stderr, "    umax = %e, energia = %e\n", umax, energy);
            divergence_detected = 1;
            break;
        }

        double dt_adv = gp.CFL_adv * dx / umax;
        double dt = gp.safety * fmin(dt_adv, dt_diff);

        step(u, un, N, dt, dx, nu);

        if (tstep % acp.print_step == 0) {
            save_data_fft(x, u, N, tstep, data_path);
            printf("  Step %d/%d (%.1f%%) | umax=%.4f | E=%.4e | dt=%.2e\n",
                   tstep, acp.nt, 100.0*tstep/acp.nt, umax, energy, dt);
        }
    }

    if (divergence_detected) {
        printf("\nSimulazione interrotta per divergenza!\n");
        free(x); free(u); free(un);
        return 1;
    }
    
    printf("Simulazione completata!\n");
    printf("File salvati in %s\n\n", data_path);
    
    free(x); free(u); free(un);
    return 0;
}
