#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Prototipi delle simulazioni
int sim_adv_conv(void);
int sim_shock_analysis(void);

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Uso: %s <mode>\n", argv[0]);
        fprintf(stderr, "  mode = sim      : Simulazione advection-diffusion standard\n");
        fprintf(stderr, "  mode = shock    : Analisi pendenza shock vs viscosità\n");
        return 1;
    }

    if (strcmp(argv[1], "sim") == 0) {
        return sim_adv_conv();
    } 
    else if (strcmp(argv[1], "shock") == 0) {
        return sim_shock_analysis();
    } 
    else {
        fprintf(stderr, "Modalità sconosciuta: %s\n", argv[1]);
        fprintf(stderr, "Usa 'sim' o 'shock'\n");
        return 1;
    }

    return 0;
}