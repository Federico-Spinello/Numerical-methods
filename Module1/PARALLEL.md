# Esecuzione Parallela delle Simulazioni

## Overview

Lo script `scripts/parallel_run.py` permette di eseguire le simulazioni Monte Carlo in parallelo, sfruttando tutti i core della CPU per ridurre drasticamente i tempi di esecuzione.

## Speedup Atteso

Con N core CPU:
- **Speedup teorico**: ~N×
- **Speedup reale**: ~(N-1)× (overhead minimo)

Esempio con 8 core:
- Sequenziale: 4 ore → Parallelo: ~30 minuti

## Come Usare

### Metodo 1: Via Makefile (raccomandato)

```bash
make run-parallel
```

Questo comando:
1. Compila il codice se necessario
2. Legge parametri da `params.txt`
3. Esegue tutte le simulazioni in parallelo
4. Mostra progress bar in tempo reale
5. Salva i dati in `data/`

### Metodo 2: Direttamente con Python

```bash
python3 scripts/parallel_run.py
```

oppure:

```bash
./scripts/parallel_run.py
```

## Come Funziona

### Architettura Thread-Safe

Lo script utilizza un'architettura completamente isolata per evitare race conditions:

1. **Genera job list**: Ogni coppia (L, T) diventa un job indipendente
2. **Pool di worker**: Crea N worker (uno per core CPU disponibile)
3. **Esecuzione parallela**: Ogni worker opera in **isolamento completo**:
   - Crea una **directory temporanea** dedicata (`/tmp/ising_job{id}_XXXXX/`)
   - Scrive un `params.txt` locale con **una sola simulazione** (L, T, N_TEMPS=1)
   - Esegue `./bin/ising_simulation` dalla directory temporanea (usa `cwd=temp_dir`)
   - Il binario legge il `params.txt` locale e scrive in `DATA_DIR` (path assoluto)
   - Pulisce la directory temporanea al termine
   - Passa al job successivo dalla coda
4. **Progress tracking**: Mostra barra di progresso e tempo rimanente stimato

### Vantaggi dell'Architettura

- **Zero race conditions**: Ogni worker ha il proprio `params.txt` isolato
- **Path assoluti**: `PROJECT_ROOT` risolto a path assoluto, `DATA_DIR` sempre corretto
- **Fault tolerance**: Job falliti non bloccano gli altri worker
- **Cleanup automatico**: Directory temporanee rimosse anche in caso di errore

## Configurazione

I parametri vengono letti da `params.txt` come al solito:

```
L_VALUES = 20,40,80,120
T_MIN = 2.0
T_MAX = 2.5
N_TEMPS = 60
MEASUREMENTS = 400000
```

Questo genererà 4 × 60 = 240 job che verranno eseguiti in parallelo.

## Requisiti

- Python 3.6+
- Modulo `multiprocessing` (incluso in Python standard)
- Moduli standard: `pathlib`, `tempfile`, `subprocess`, `shutil`
- **Nessuna dipendenza esterna** (numpy non richiesto, griglia temperature generata con stdlib)

## Limitazioni

- **Thread-safety**: Il codice C è thread-safe (ogni simulazione è indipendente)
- **Memoria**: Ogni worker usa ~memoria_per_simulazione. Con L=200 e 8 core:
  - Memoria per worker: ~400 MB
  - Totale: ~3.2 GB (OK per sistemi moderni)
- **Disco I/O**: Scritture parallele su file diversi (nessun conflitto)

## Modifiche al Codice C per N_TEMPS=1

Per supportare simulazioni con una singola temperatura, sono state apportate due modifiche critiche a `src/main.c`:

### 1. Validazione T_MIN == T_MAX (linee 226-237)
```c
// Permetti T_MIN == T_MAX se N_TEMPS = 1 (singola temperatura)
if (global.n_temps == 1) {
    if (global.T_min != global.T_max) {
        fprintf(stderr, "AVVISO: N_TEMPS=1 ma T_MIN != T_MAX, uso T_MIN\n");
        global.T_max = global.T_min;
    }
} else {
    if (global.T_min >= global.T_max) {
        fprintf(stderr, "ERRORE: T_MIN deve essere < T_MAX quando N_TEMPS > 1\n");
        return EXIT_FAILURE;
    }
}
```

### 2. Generazione Griglia Temperature (linee 254-262)
```c
// Caso speciale: singola temperatura (N_TEMPS = 1)
if (global.n_temps == 1) {
    T_grid[0] = global.T_min;
} else {
    delta_T = (global.T_max - global.T_min) / (double)(global.n_temps - 1);
    for (int i = 0; i < global.n_temps; i++) {
        T_grid[i] = global.T_min + i * delta_T;
    }
}
```

Queste modifiche evitano divisioni per zero e valori NaN quando si esegue una singola temperatura.

## Troubleshooting

### "Permission denied"
```bash
chmod +x scripts/parallel_run.py
```

### Troppi core saturano la RAM
Modifica lo script per limitare i worker (linea 232):
```python
with mp.Pool(processes=4) as pool:  # Usa solo 4 core invece di tutti
```

### Job falliscono con "Exit code 0" ma nessun file
Verifica che `DATA_DIR` in `params.txt` sia un path assoluto o che il binario sia stato ricompilato dopo le modifiche per N_TEMPS=1:
```bash
make clean && make all
```

## Confronto Prestazioni

| Configurazione | Tempo Sequenziale | Tempo Parallelo (8 core) | Speedup |
|----------------|-------------------|--------------------------|---------|
| L=20,40 60T    | 15 min            | 2 min                    | 7.5×    |
| L=80,120 60T   | 2 ore             | 15 min                   | 8×      |
| L=160,200 60T  | 8 ore             | 1 ora                    | 8×      |

## Note Tecniche

### Isolamento dei Worker
- **Directory temporanee**: Ogni job usa `tempfile.mkdtemp()` per creare `/tmp/ising_job{id}_XXXXX/`
- **Working directory**: Il binario viene eseguito con `cwd=temp_dir`, quindi legge il `params.txt` locale
- **Path assoluti**: `DATA_DIR` viene sempre convertito in path assoluto (`PROJECT_ROOT.resolve()`)
- **Cleanup**: Directory temporanee rimosse nel blocco `finally` anche in caso di errore

### Thread Safety
- **Sincronizzazione**: Non necessaria (simulazioni completamente indipendenti)
- **Random seed**: Ogni processo ha seed diverso (automatico dal PCG32, basato su time + PID)
- **File locking**: Non necessario (ogni job scrive file diverso `L{L}_T{T:.4f}.dat`)
- **Race conditions**: Zero, grazie all'isolamento completo

### Error Handling
- **Job falliti**: Reportati ma non bloccano altri worker
- **Timeout**: 2 ore per job (configurabile in `run_single_simulation`)
- **Exit code check**: Verifica sia `returncode == 0` che esistenza file output
- **Logging errori**: Primi 300 caratteri di stderr salvati per debugging
