#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "params.h"

// Funzione helper per trovare una sezione
static int goto_section(FILE *fp, const char *section) {
    char line[256];
    char target[256];
    snprintf(target, sizeof(target), "[%s]", section);
    
    rewind(fp);
    while (fgets(line, sizeof(line), fp)) {
        // Rimuovi spazi bianchi
        char *ptr = line;
        while (*ptr == ' ' || *ptr == '\t') ptr++;
        if (strncmp(ptr, target, strlen(target)) == 0) {
            return 1;
        }
    }
    return 0;
}

// Funzione helper per leggere un parametro
static int read_param(FILE *fp, const char *param_name, char *value) {
    char line[256];
    long pos = ftell(fp);
    
    while (fgets(line, sizeof(line), fp)) {
        // Stop se inizia una nuova sezione
        if (line[0] == '[') {
            fseek(fp, pos, SEEK_SET);
            return 0;
        }
        
        // Salta commenti e righe vuote
        char *ptr = line;
        while (*ptr == ' ' || *ptr == '\t') ptr++;
        if (*ptr == '#' || *ptr == '\n') {
            pos = ftell(fp);
            continue;
        }
        
        // Cerca il parametro
        char name[64];
        if (sscanf(line, "%s = %s", name, value) == 2) {
            if (strcmp(name, param_name) == 0) {
                // Rimuovi eventuali commenti inline
                char *comment = strchr(value, '#');
                if (comment) *comment = '\0';
                return 1;
            }
        }
        pos = ftell(fp);
    }
    return 0;
}

int read_global_params(const char *filename, GlobalParams *params) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Errore: impossibile aprire %s\n", filename);
        return 0;
    }
    
    if (!goto_section(fp, "GLOBAL")) {
        fprintf(stderr, "Errore: sezione [GLOBAL] non trovata\n");
        fclose(fp);
        return 0;
    }
    
    char value[256];
    
    if (read_param(fp, "nx", value)) params->nx = atoi(value);
    if (read_param(fp, "L", value)) params->L = atof(value);
    if (read_param(fp, "CFL_diff", value)) params->CFL_diff = atof(value);
    if (read_param(fp, "CFL_adv", value)) params->CFL_adv = atof(value);
    if (read_param(fp, "safety", value)) params->safety = atof(value);
    
    fclose(fp);
    return 1;
}

int read_advconv_params(const char *filename, AdvConvParams *params) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Errore: impossibile aprire %s\n", filename);
        return 0;
    }
    
    if (!goto_section(fp, "ADVECTION_CONVECTION")) {
        fprintf(stderr, "Errore: sezione [ADVECTION_CONVECTION] non trovata\n");
        fclose(fp);
        return 0;
    }
    
    char value[256];
    
    if (read_param(fp, "nt", value)) params->nt = atoi(value);
    if (read_param(fp, "print_step", value)) params->print_step = atoi(value);
    if (read_param(fp, "nu_factor", value)) params->nu_factor = atof(value);
    
    fclose(fp);
    return 1;
}

int read_shock_params(const char *filename, ShockAnalysisParams *params) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Errore: impossibile aprire %s\n", filename);
        return 0;
    }
    
    if (!goto_section(fp, "SHOCK_ANALYSIS")) {
        fprintf(stderr, "Errore: sezione [SHOCK_ANALYSIS] non trovata\n");
        fclose(fp);
        return 0;
    }
    
    char value[256];
    
    if (read_param(fp, "nt", value)) params->nt = atoi(value);
    if (read_param(fp, "nu_min_factor", value)) params->nu_min_factor = atof(value);
    if (read_param(fp, "nu_max_factor", value)) params->nu_max_factor = atof(value);
    if (read_param(fp, "n_nu", value)) params->n_nu = atoi(value);
    
    fclose(fp);
    return 1;
}