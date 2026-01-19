/*
 * thermalization_test.c - Programma per studiare la termalizzazione
 *
 * Esegue DUE simulazioni (cold start e hot start) e salva E(t), m(t)
 * ad ogni step durante la termalizzazione per analisi.
 *
 * Parametri letti da params.txt:
 *   THERM_L      - Dimensione reticolo
 *   THERM_T      - Temperatura
 *   THERM_NSTEPS - Numero di cluster updates
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "../include/geometry.h"
#include "../include/random.h"

#define DIM 2
#define OUTPUT_FILE_COLD "data/thermalization_cold.dat"
#define OUTPUT_FILE_HOT "data/thermalization_hot.dat"
#define PARAMS_FILE "params.txt"

/* Funzione per leggere parametri da params.txt */
static int read_params(int *L, double *T, int *nsteps) {
    FILE *fp = fopen(PARAMS_FILE, "r");
    if (!fp) {
        fprintf(stderr, "ERRORE: impossibile aprire %s\n", PARAMS_FILE);
        return 0;
    }

    char line[256];
    int found_L = 0, found_T = 0, found_N = 0;

    while (fgets(line, sizeof(line), fp)) {
        // Salta commenti e righe vuote
        if (line[0] == '#' || line[0] == '\n') continue;

        // Rimuovi newline
        line[strcspn(line, "\n")] = 0;

        // Parsing
        if (strstr(line, "THERM_L") && strstr(line, "=")) {
            char *eq = strchr(line, '=');
            if (eq) {
                *L = atoi(eq + 1);
                found_L = 1;
            }
        }
        else if (strstr(line, "THERM_T") && strstr(line, "=")) {
            char *eq = strchr(line, '=');
            if (eq) {
                *T = atof(eq + 1);
                found_T = 1;
            }
        }
        else if (strstr(line, "THERM_NSTEPS") && strstr(line, "=")) {
            char *eq = strchr(line, '=');
            if (eq) {
                *nsteps = atoi(eq + 1);
                found_N = 1;
            }
        }
    }

    fclose(fp);

    if (!found_L || !found_T || !found_N) {
        fprintf(stderr, "ERRORE: parametri THERM_L, THERM_T, THERM_NSTEPS non trovati in %s\n", PARAMS_FILE);
        return 0;
    }

    return 1;
}

/* Calcola magnetizzazione per sito */
static double magn(int const * const restrict lattice, long int volume) {
    long int r, sum = 0;
    for(r=0; r<volume; r++) {
        sum += lattice[r];
    }
    return (double)sum / (double)volume;
}

/* Calcola energia per sito */
static double energy(int const * const restrict lattice,
                     long int const * const restrict nnp,
                     long int volume) {
    long int r, sum = 0;
    for(r=0; r<volume; r++) {
        for(int i=0; i<DIM; i++) {
            sum += -lattice[r] * lattice[nnp[dirgeo(r, i, volume)]];
        }
    }
    return (double)sum / (double)volume;
}

/* Costruzione cluster (versione semplificata non-ricorsiva) */
static void build_cluster_norec(int const * const restrict lattice,
                                 long int r_seed,
                                 int * restrict occup,
                                 long int * restrict pointtoocc,
                                 long int * restrict clustersize,
                                 long int const * const restrict nnp,
                                 long int const * const restrict nnm,
                                 long int volume,
                                 double prob) {
    (void) r_seed;
    long int oldcs = 0;
    long int oldcsnew, index, r1;

    while(*clustersize > oldcs) {
        oldcsnew = *clustersize;

        for(index=oldcs; index<oldcsnew; index++) {
            r1 = pointtoocc[index];

            for(int i=0; i<DIM; i++) {
                // Avanti
                long int nn_pos = nnp[dirgeo(r1, i, volume)];
                if(occup[nn_pos]==0 && lattice[r1]*lattice[nn_pos]==1) {
                    if(myrand() < prob) {
                        occup[nn_pos] = 1;
                        pointtoocc[*clustersize] = nn_pos;
                        (*clustersize)++;
                    }
                }

                // Indietro
                long int nn_neg = nnm[dirgeo(r1, i, volume)];
                if(occup[nn_neg]==0 && lattice[r1]*lattice[nn_neg]==1) {
                    if(myrand() < prob) {
                        occup[nn_neg] = 1;
                        pointtoocc[*clustersize] = nn_neg;
                        (*clustersize)++;
                    }
                }
            }
        }
        oldcs = oldcsnew;
    }
}

/* Funzione per eseguire una simulazione di termalizzazione */
static void run_thermalization(int L, double T, int nsteps, int hot_start, const char *output_file,
                               long int *nnp, long int *nnm) {
    long int volume = (long int)L * (long int)L;
    double beta = 1.0 / T;
    double prob = 1.0 - exp(-2.0 * beta);

    // Alloca memoria per il reticolo e array ausiliari
    int *lattice = (int *)malloc(volume * sizeof(int));
    int *occup = (int *)malloc(volume * sizeof(int));
    long int *pointtoocc = (long int *)malloc(volume * sizeof(long int));

    if(!lattice || !occup || !pointtoocc) {
        fprintf(stderr, "ERRORE: allocazione memoria fallita\n");
        return;
    }

    // Inizializza reticolo
    if (hot_start) {
        // Hot start: spin casuali
        printf("Inizializzazione: HOT START (spin casuali)\n");
        for(long int r=0; r<volume; r++) {
            lattice[r] = (myrand() < 0.5) ? 1 : -1;
        }
    } else {
        // Cold start: tutti spin +1
        printf("Inizializzazione: COLD START (tutti spin +1)\n");
        for(long int r=0; r<volume; r++) {
            lattice[r] = 1;
        }
    }

    // Apri file output
    FILE *fp = fopen(output_file, "w");
    if(!fp) {
        fprintf(stderr, "ERRORE: impossibile aprire %s\n", output_file);
        free(lattice);
        free(occup);
        free(pointtoocc);
        return;
    }

    fprintf(fp, "# Dati termalizzazione: L=%d, T=%.4f, %s\n", L, T, hot_start ? "HOT START" : "COLD START");
    fprintf(fp, "# Colonne: step E m\n");

    printf("Esecuzione %d cluster updates...\n", nsteps);

    // Loop termalizzazione CON salvataggio dati
    for(int iter=0; iter<nsteps; iter++) {
        // Azzera occupazione
        for(long int r=0; r<volume; r++) {
            occup[r] = 0;
        }
        long int clustersize = 0;

        // Scegli seed random
        long int r = (long int)((double)volume * myrand());
        occup[r] = 1;
        pointtoocc[clustersize] = r;
        clustersize++;

        // Costruisci cluster
        build_cluster_norec(lattice, r, occup, pointtoocc, &clustersize, nnp, nnm, volume, prob);

        // Flippa cluster
        for(long int i=0; i<clustersize; i++) {
            lattice[pointtoocc[i]] = -lattice[pointtoocc[i]];
        }

        // Calcola osservabili
        double E = energy(lattice, nnp, volume);
        double m = magn(lattice, volume);

        // Salva dati
        fprintf(fp, "%d %.8f %.8f\n", iter, E, fabs(m));

        // Progress indicator
        if((iter+1) % 1000 == 0) {
            printf("  Step %d/%d completato\n", iter+1, nsteps);
        }
    }

    fclose(fp);
    printf("Dati salvati in: %s\n\n", output_file);

    // Cleanup
    free(lattice);
    free(occup);
    free(pointtoocc);
}

int main(void) {
    // Leggi parametri da params.txt
    int L, nsteps;
    double T;

    if (!read_params(&L, &T, &nsteps)) {
        fprintf(stderr, "ERRORE: impossibile leggere parametri da %s\n", PARAMS_FILE);
        return 1;
    }

    long int volume = (long int)L * (long int)L;

    // Cold start: 410000 steps (per grafico completo con zone rosse/verdi)
    // Hot start: usa nsteps da params.txt (per grafico zoom overlay)
    int nsteps_cold = 410000;
    int nsteps_hot = nsteps;  // tipicamente 15000

    printf("========================================\n");
    printf("STUDIO TERMALIZZAZIONE\n");
    printf("========================================\n");
    printf("Parametri letti da %s:\n", PARAMS_FILE);
    printf("  L = %d\n", L);
    printf("  T = %.4f\n", T);
    printf("  Volume = %ld\n", volume);
    printf("  Steps COLD: %d (per grafico completo)\n", nsteps_cold);
    printf("  Steps HOT:  %d (per grafico zoom)\n", nsteps_hot);
    printf("  Output COLD: %s\n", OUTPUT_FILE_COLD);
    printf("  Output HOT:  %s\n\n", OUTPUT_FILE_HOT);

    // Alloca memoria per geometria (condivisa tra le due simulazioni)
    long int *nnp = (long int *)malloc(DIM * volume * sizeof(long int));
    long int *nnm = (long int *)malloc(DIM * volume * sizeof(long int));

    if(!nnp || !nnm) {
        fprintf(stderr, "ERRORE: allocazione memoria fallita\n");
        return 1;
    }

    // Inizializza geometria
    init_neighbors(nnp, nnm, L, DIM);

    // Inizializza RNG
    myrand_init((unsigned long int)time(NULL), 12345);

    // ========== COLD START (410000 steps) ==========
    printf("----------------------------------------\n");
    printf("SIMULAZIONE 1: COLD START (%d steps)\n", nsteps_cold);
    printf("----------------------------------------\n");
    run_thermalization(L, T, nsteps_cold, 0, OUTPUT_FILE_COLD, nnp, nnm);

    // Re-inizializza RNG per avere sequenza statistica indipendente
    myrand_init((unsigned long int)time(NULL) + 1000, 67890);

    // ========== HOT START (nsteps steps) ==========
    printf("----------------------------------------\n");
    printf("SIMULAZIONE 2: HOT START (%d steps)\n", nsteps_hot);
    printf("----------------------------------------\n");
    run_thermalization(L, T, nsteps_hot, 1, OUTPUT_FILE_HOT, nnp, nnm);

    // Cleanup
    free(nnp);
    free(nnm);

    printf("========================================\n");
    printf("COMPLETATO!\n");
    printf("========================================\n");

    return 0;
}
