#ifndef ISING_WOLFF_H
#define ISING_WOLFF_H

/*
 * ising_wolff.h - Dichiarazioni per simulazione Ising 2D con algoritmo di Wolff
 *
 * Questo modulo implementa una singola simulazione Monte Carlo del modello di
 * Ising 2D utilizzando l'algoritmo di Wolff (cluster updates).
 */

/*
 * Parametri simulazione
 */
typedef struct {
    int L;                    // Dimensione reticolo (L×L)
    double T;                 // Temperatura (in unità J/k_B)
    long int thermalization;  // Numero update termalizzazione
    long int measurements;    // Numero update per misure
    char output_file[256];    // File output dati
} SimulationParams;

/*
 * Risultati simulazione
 */
typedef struct {
    double mag_mean;          // <|m|>
    double mag2_mean;         // <m²>
    double mag4_mean;         // <m⁴>
    double energy_mean;       // <E>
    double energy2_mean;      // <E²>
    double mag_err;           // errore standard su <|m|>
    double energy_err;        // errore standard su <E>
    double chi_err;           // errore standard su chi
    double C_err;             // errore standard su C
} SimulationResults;

/*
 * Esegue una singola simulazione Ising 2D con Wolff
 *
 * Parametri:
 *   params: Struttura con parametri simulazione
 *   results: Struttura dove salvare risultati (output)
 *
 * Ritorna:
 *   0 se successo, -1 se errore
 */
int run_ising_simulation(const SimulationParams *params, SimulationResults *results);

#endif // ISING_WOLFF_H
