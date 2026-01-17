/*
 * main.c - Wrapper per singola simulazione Ising 2D
 *
 * Versione snellita per esecuzione parallela tramite parallel_run.py.
 * Ogni invocazione esegue UNA sola simulazione (L, T) specificata in params.txt.
 *
 * Workflow:
 *   1. Legge parametri singolo job da params.txt
 *   2. Esegue run_ising_simulation() una volta
 *   3. Termina
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "../include/ising_wolff.h"

#define MAX_LINE_LENGTH 256

/*
 * Parametri per singolo job (no array, no range)
 */
typedef struct {
    int L;                    // Dimensione reticolo (singolo valore)
    double T;                 // Temperatura (singolo valore)
    long thermalization;      // Numero cluster updates termalizzazione
    long measurements;        // Numero cluster updates misure
    char output_file[256];    // Path output file
} SingleJobParams;

/*
 * Parsing di una riga "KEY = value"
 * Ritorna 1 se la chiave corrisponde, 0 altrimenti
 */
static int parse_line(const char *line, const char *key, char *value) {
    // Salta commenti e linee vuote
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\0') {
        return 0;
    }

    // Cerca '='
    char *equals = strchr(line, '=');
    if (equals == NULL) {
        return 0;
    }

    // Estrae chiave
    char line_key[128];
    size_t key_len = equals - line;
    if (key_len >= 128) return 0;

    strncpy(line_key, line, key_len);
    line_key[key_len] = '\0';

    // Trim spazi dalla chiave
    char *k = line_key;
    while (*k == ' ' || *k == '\t') k++;
    char *k_end = line_key + key_len - 1;
    while (k_end > k && (*k_end == ' ' || *k_end == '\t')) {
        *k_end = '\0';
        k_end--;
    }

    // Confronta chiave
    if (strcmp(k, key) != 0) {
        return 0;
    }

    // Estrae valore
    char *val = equals + 1;
    while (*val == ' ' || *val == '\t') val++;

    // Rimuove newline e spazi finali
    strcpy(value, val);
    char *newline = strchr(value, '\n');
    if (newline) *newline = '\0';

    size_t len = strlen(value);
    while (len > 0 && (value[len-1] == ' ' || value[len-1] == '\t')) {
        value[--len] = '\0';
    }

    return 1;
}

/*
 * Carica parametri per singolo job da params.txt
 */
static int load_single_job_params(const char *filename, SingleJobParams *params) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "ERRORE: Impossibile aprire %s\n", filename);
        return 0;
    }

    char line[MAX_LINE_LENGTH];
    char value[MAX_LINE_LENGTH];

    // Valori di default
    params->L = 0;
    params->T = 0.0;
    params->thermalization = 10000;
    params->measurements = 100000;
    strcpy(params->output_file, "output.dat");

    // Leggi file riga per riga
    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) {
        if (parse_line(line, "L", value)) {
            params->L = atoi(value);
        }
        else if (parse_line(line, "T", value)) {
            params->T = atof(value);
        }
        else if (parse_line(line, "THERMALIZATION", value)) {
            params->thermalization = atol(value);
        }
        else if (parse_line(line, "MEASUREMENTS", value)) {
            params->measurements = atol(value);
        }
        else if (parse_line(line, "OUTPUT_FILE", value)) {
            strncpy(params->output_file, value, 255);
            params->output_file[255] = '\0';
        }
    }

    fclose(fp);

    // Validazione
    if (params->L <= 0) {
        fprintf(stderr, "ERRORE: L non valido (%d)\n", params->L);
        return 0;
    }
    if (params->T <= 0.0) {
        fprintf(stderr, "ERRORE: T non valido (%.4f)\n", params->T);
        return 0;
    }

    return 1;
}

/*
 * Crea directory se non esiste
 */
static void ensure_output_dir(const char *filepath) {
    char dirpath[256];
    strncpy(dirpath, filepath, 255);
    dirpath[255] = '\0';

    // Trova l'ultima '/'
    char *last_slash = strrchr(dirpath, '/');
    if (last_slash != NULL) {
        *last_slash = '\0';

        struct stat st = {0};
        if (stat(dirpath, &st) == -1) {
            mkdir(dirpath, 0755);
        }
    }
}

/*
 * Main - Esegue singola simulazione
 */
int main(int argc, char **argv) {
    SingleJobParams job;
    SimulationParams sim_params;
    SimulationResults results;

    // Carica parametri singolo job
    if (!load_single_job_params("params.txt", &job)) {
        fprintf(stderr, "ERRORE: Impossibile caricare parametri\n");
        return EXIT_FAILURE;
    }

    // Info esecuzione
    printf("Ising 2D - Singola simulazione\n");
    printf("  L = %d\n", job.L);
    printf("  T = %.4f\n", job.T);
    printf("  Termalizzazione: %ld cluster updates\n", job.thermalization);
    printf("  Misurazioni: %ld cluster updates\n", job.measurements);
    printf("  Output: %s\n", job.output_file);

    // Assicura che la directory di output esista
    ensure_output_dir(job.output_file);

    // Prepara parametri simulazione
    sim_params.L = job.L;
    sim_params.T = job.T;
    sim_params.thermalization = job.thermalization;
    sim_params.measurements = job.measurements;
    strncpy(sim_params.output_file, job.output_file, 255);
    sim_params.output_file[255] = '\0';

    // Esegui simulazione
    printf("\nEsecuzione simulazione...\n");
    if (run_ising_simulation(&sim_params, &results) != 0) {
        fprintf(stderr, "ERRORE durante simulazione L=%d T=%.4f\n", job.L, job.T);
        return EXIT_FAILURE;
    }

    printf("✓ Simulazione completata\n");
    printf("  Dati salvati: %s\n", job.output_file);

    return EXIT_SUCCESS;
}
