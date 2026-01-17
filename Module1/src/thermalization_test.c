/*
 * thermalization_test.c - Programma per studiare la termalizzazione
 *
 * Esegue UNA simulazione e salva E(t), m(t) ad ogni step
 * durante la termalizzazione per analisi.
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
#define OUTPUT_FILE "data/thermalization_data.dat"
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

int main(void) {
    // Leggi parametri da params.txt
    int L, nsteps;
    double T;

    if (!read_params(&L, &T, &nsteps)) {
        fprintf(stderr, "ERRORE: impossibile leggere parametri da %s\n", PARAMS_FILE);
        return 1;
    }

    long int volume = (long int)L * (long int)L;
    double beta = 1.0 / T;
    double prob = 1.0 - exp(-2.0 * beta);

    printf("========================================\n");
    printf("STUDIO TERMALIZZAZIONE\n");
    printf("========================================\n");
    printf("Parametri letti da %s:\n", PARAMS_FILE);
    printf("  L = %d\n", L);
    printf("  T = %.4f\n", T);
    printf("  Volume = %ld\n", volume);
    printf("  Steps = %d\n", nsteps);
    printf("  Output: %s\n\n", OUTPUT_FILE);

    // Alloca memoria
    int *lattice = (int *)malloc(volume * sizeof(int));
    long int *nnp = (long int *)malloc(DIM * volume * sizeof(long int));
    long int *nnm = (long int *)malloc(DIM * volume * sizeof(long int));
    int *occup = (int *)malloc(volume * sizeof(int));
    long int *pointtoocc = (long int *)malloc(volume * sizeof(long int));

    if(!lattice || !nnp || !nnm || !occup || !pointtoocc) {
        fprintf(stderr, "ERRORE: allocazione memoria fallita\n");
        return 1;
    }

    // Inizializza RNG
    myrand_init((unsigned long int)time(NULL), 12345);

    // Inizializza geometria
    init_neighbors(nnp, nnm, L, DIM);

    // Cold start: tutti spin +1
    printf("Inizializzazione: cold start (tutti spin +1)\n");
    for(long int r=0; r<volume; r++) {
        lattice[r] = 1;
    }

    // Apri file output
    FILE *fp = fopen(OUTPUT_FILE, "w");
    if(!fp) {
        fprintf(stderr, "ERRORE: impossibile aprire %s\n", OUTPUT_FILE);
        return 1;
    }

    fprintf(fp, "# Dati termalizzazione: L=%d, T=%.4f\n", L, T);
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
    printf("\n✓ Dati salvati in: %s\n", OUTPUT_FILE);

    // Cleanup
    free(lattice);
    free(nnp);
    free(nnm);
    free(occup);
    free(pointtoocc);

    return 0;
}
