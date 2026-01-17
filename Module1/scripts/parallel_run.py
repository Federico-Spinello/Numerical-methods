#!/usr/bin/env python3
"""
parallel_run.py - Esecuzione parallela delle simulazioni Ising 2D

Questo script permette di eseguire le simulazioni Monte Carlo in parallelo,
sfruttando tutti i core della CPU per ridurre drasticamente i tempi di esecuzione.

Architettura Thread-Safe:
-------------------------
1. Genera job list: Ogni coppia (L, T) diventa un job indipendente
2. Pool di worker: Crea N worker (uno per core CPU disponibile)
3. Esecuzione parallela: Ogni worker opera in isolamento completo:
   - Crea una directory temporanea dedicata (/tmp/ising_job{id}_XXXXX/)
   - Scrive un params.txt locale con una sola simulazione (L, T, N_TEMPS=1)
   - Esegue ./bin/ising_simulation dalla directory temporanea
   - Pulisce la directory temporanea al termine
4. Progress tracking: Mostra barra di progresso e tempo rimanente stimato

Vantaggi:
---------
- Zero race conditions: Ogni worker ha il proprio params.txt isolato
- Path assoluti: DATA_DIR sempre corretto
- Fault tolerance: Job falliti non bloccano gli altri worker
- Cleanup automatico: Directory temporanee rimosse anche in caso di errore

Uso:
----
    python3 scripts/parallel_run.py
    # oppure
    make run-parallel
"""

import os
import sys
import time
import shutil
import tempfile
import subprocess
import multiprocessing as mp
from pathlib import Path
from datetime import datetime


# ============================================================================
# CONFIGURAZIONE
# ============================================================================

# Directory del progetto (risolto a path assoluto)
PROJECT_ROOT = Path(__file__).parent.parent.resolve()
BIN_PATH = PROJECT_ROOT / "bin" / "ising_simulation"
PARAMS_FILE = PROJECT_ROOT / "params.txt"
DATA_DIR = PROJECT_ROOT / "data"

# Timeout per singola simulazione (2 ore)
JOB_TIMEOUT = 7200


# ============================================================================
# UTILITA'
# ============================================================================

from pathlib import Path

def parse_params():
    """
    Legge i parametri dal file params.txt esistente e li converte nei tipi corretti.

    Returns:
        dict: dizionario con i parametri
    """

    # Percorso del file dei parametri
    PARAMS_FILE = PROJECT_ROOT / "params.txt"

    # Dizionario dove memorizzare i parametri
    params = {}

    # Apro il file in lettura
    with open(PARAMS_FILE, 'r') as f:
        for line in f:
            # Rimuovo spazi iniziali e finali
            line_stripped = line.strip()

            # Ignoro righe vuote o commenti
            if line_stripped == "" or line_stripped.startswith("#"):
                continue

            # Controllo che la riga contenga '='
            if '=' not in line_stripped:
                print(f"ATTENZIONE: riga malformata ignorata: {line_stripped}")
                continue

            # Separazione chiave-valore
            key_raw, value_raw = line_stripped.split('=', 1)

            # Pulizia spazi
            key = key_raw.strip()
            value = value_raw.strip()

            # Conversione dei valori in base al tipo atteso
            if key == 'L_VALUES':
                # Lista di interi separati da virgola
                values_list = value.split(',')
                int_list = []
                for v in values_list:
                    v_clean = v.strip()
                    int_list.append(int(v_clean))
                params[key] = int_list

            elif key in ('T_MIN', 'T_MAX'):
                # Float
                params[key] = float(value)

            elif key in ('N_TEMPS', 'THERMALIZATION', 'MEASUREMENTS'):
                # Interi
                params[key] = int(value)

            elif key == 'DATA_DIR':
                # Percorso, convertito in path assoluto
                data_path = Path(value)
                if not data_path.is_absolute():
                    data_path = PROJECT_ROOT / data_path
                params[key] = str(data_path.resolve())

            else:
                # Chiave sconosciuta
                print(f"ATTENZIONE: chiave sconosciuta '{key}' ignorata")

    # Ritorno il dizionario con tutti i parametri
    return params



def generate_temperature_grid(T_min, T_max, n_temps):
    """
    Genera una griglia uniforme di temperature.

    Args:
        T_min: temperatura minima
        T_max: temperatura massima
        n_temps: numero di temperature

    Returns:
        list: lista di temperature
    """

    # Calcolo del passo di temperatura
    delta_T = (T_max - T_min) / (n_temps - 1)

    # Costruzione della lista di temperature
    temperatures = []
    for i in range(n_temps):
        T = T_min + i * delta_T
        temperatures.append(T)

    return temperatures


def generate_jobs(params):
    """
    Genera la lista di job (L, T) da eseguire.

    Args:
        params: dizionario parametri

    Returns:
        list: lista di tuple (job_id, L, T)
    """
    L_values = params['L_VALUES']
    temperatures = generate_temperature_grid(
        params['T_MIN'],
        params['T_MAX'],
        params['N_TEMPS']
    )

    jobs = []
    job_id = 0
    for L in L_values:
        for T in temperatures:
            jobs.append((job_id, L, T))
            job_id += 1

    return jobs


def write_job_params(temp_dir, L, T, params):
    """
    Scrive il file params.txt per un singolo job.

    Args:
        temp_dir: directory temporanea del job
        L: dimensione reticolo
        T: temperatura
        params: dizionario parametri globali
    """
    output_file = f"{params['DATA_DIR']}/L{L}_T{T:.4f}.dat"
    params_content = f"""# Parametri per singola simulazione (generato da parallel_run.py)
L = {L}
T = {T}
THERMALIZATION = {params['THERMALIZATION']}
MEASUREMENTS = {params['MEASUREMENTS']}
OUTPUT_FILE = {output_file}
"""
    params_file = Path(temp_dir) / "params.txt"
    with open(params_file, 'w') as f:
        f.write(params_content)


def run_single_simulation(job_args):
    """
    Esegue una singola simulazione in una directory temporanea isolata.

    Args:
        job_args: tuple (job_id, L, T, params, bin_path)

    Returns:
        dict: risultato con job_id, success, message, duration
    """
    job_id, L, T, params, bin_path = job_args

    start_time = time.time()
    temp_dir = None

    try:
        # Crea directory temporanea isolata
        temp_dir = tempfile.mkdtemp(prefix=f"ising_job{job_id}_")

        # Scrivi params.txt locale
        write_job_params(temp_dir, L, T, params)

        # Nome file output atteso
        output_file = Path(params['DATA_DIR']) / f"L{L}_T{T:.4f}.dat"

        # Esegui simulazione dalla directory temporanea
        result = subprocess.run(
            [str(bin_path)],
            cwd=temp_dir,
            capture_output=True,
            text=True,
            timeout=JOB_TIMEOUT
        )

        duration = time.time() - start_time

        # Verifica successo
        if result.returncode != 0:
            # Estrai il messaggio di errore (tronca se troppo lungo)
            if result.stderr:
                error_msg = result.stderr[:300]
            else:
                error_msg = "Unknown error"
            return {
                'job_id': job_id,
                'L': L,
                'T': T,
                'success': False,
                'message': f"Exit code {result.returncode}: {error_msg}",
                'duration': duration
            }

        # Verifica che il file output esista
        if not output_file.exists():
            return {
                'job_id': job_id,
                'L': L,
                'T': T,
                'success': False,
                'message': f"Output file not found: {output_file}",
                'duration': duration
            }

        return {
            'job_id': job_id,
            'L': L,
            'T': T,
            'success': True,
            'message': "OK",
            'duration': duration
        }

    except subprocess.TimeoutExpired:
        duration = time.time() - start_time
        return {
            'job_id': job_id,
            'L': L,
            'T': T,
            'success': False,
            'message': f"Timeout after {JOB_TIMEOUT}s",
            'duration': duration
        }

    except Exception as e:
        duration = time.time() - start_time
        return {
            'job_id': job_id,
            'L': L,
            'T': T,
            'success': False,
            'message': str(e),
            'duration': duration
        }

    finally:
        # Cleanup directory temporanea
        if temp_dir and os.path.exists(temp_dir):
            try:
                shutil.rmtree(temp_dir)
            except Exception:
                pass


def print_progress_bar(current, total, width=50, prefix='', suffix=''):
    """
    Stampa una barra di progresso sulla stessa riga.

    Args:
        current: numero di elementi completati
        total: numero totale di elementi
        width: larghezza della barra in caratteri
        prefix: testo da mostrare prima della barra
        suffix: testo da mostrare dopo la barra
    """
    # Calcola la percentuale di completamento
    percent = current / total

    # Calcola quanti caratteri della barra riempire
    filled = int(width * percent)

    # Costruisce la barra con '=' per la parte completata e '-' per quella rimanente
    bar_filled = '=' * filled
    bar_empty = '-' * (width - filled)
    bar = bar_filled + bar_empty

    # Calcola la percentuale come stringa
    percent_str = f"{percent * 100:.1f}%"

    # Stampa la barra sovrascrivendo la riga precedente
    output = f'\r{prefix} [{bar}] {current}/{total} ({percent_str}) {suffix}'
    sys.stdout.write(output)
    sys.stdout.flush()


def format_duration(seconds):
    """
    Formatta una durata in secondi in un formato leggibile.

    Args:
        seconds: durata in secondi

    Returns:
        str: durata formattata (es. "45.2s", "3m 20s", "1h 15m")
    """
    # Caso 1: meno di un minuto -> mostra secondi con decimale
    if seconds < 60:
        return f"{seconds:.1f}s"

    # Caso 2: meno di un'ora -> mostra minuti e secondi
    elif seconds < 3600:
        minuti = int(seconds // 60)
        secondi = int(seconds % 60)
        return f"{minuti}m {secondi}s"

    # Caso 3: un'ora o più -> mostra ore e minuti
    else:
        ore = int(seconds // 3600)
        secondi_rimanenti = seconds % 3600
        minuti = int(secondi_rimanenti // 60)
        return f"{ore}h {minuti}m"


# ============================================================================
# MAIN
# ============================================================================

def main():
    """Funzione principale."""
    print("=" * 70)
    print("  SIMULAZIONE PARALLELA ISING 2D - ALGORITMO DI WOLFF")
    print("=" * 70)
    print(f"  Data: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Project root: {PROJECT_ROOT}")
    print()

    # Verifica che il binario esista
    if not BIN_PATH.exists():
        print(f"ERRORE: Binario non trovato: {BIN_PATH}")
        print("Esegui 'make all' per compilare il codice.")
        sys.exit(1)

    # Carica parametri
    params = parse_params()


    # Crea directory data se non esiste
    data_dir = Path(params['DATA_DIR'])
    data_dir.mkdir(parents=True, exist_ok=True)

    # Genera lista job
    jobs = generate_jobs(params)
    n_jobs = len(jobs)

    # Info simulazione
    n_cpus = mp.cpu_count()
    print(f"  Dimensioni L: {params['L_VALUES']}")
    print(f"  Temperature: {params['N_TEMPS']} punti in [{params['T_MIN']}, {params['T_MAX']}]")
    print(f"  Termalizzazione: {params['THERMALIZATION']} steps")
    print(f"  Misure: {params['MEASUREMENTS']} steps")
    print(f"  Directory output: {params['DATA_DIR']}")
    print()
    print(f"  Job totali: {n_jobs}")
    print(f"  CPU disponibili: {n_cpus}")
    print(f"  Worker paralleli: {n_cpus}")
    print()

    # Prepara argomenti per i worker
    job_args = []
    for job_id, L, T in jobs:
        args_tuple = (job_id, L, T, params, str(BIN_PATH))
        job_args.append(args_tuple)

    # Esegui in parallelo
    print("Avvio simulazioni...")
    print("-" * 70)

    start_time = time.time()
    completed = 0
    failed = 0
    failed_jobs = []

    # Usa Pool per esecuzione parallela
    with mp.Pool(processes=n_cpus) as pool:
        # imap_unordered per risultati in ordine di completamento
        for result in pool.imap_unordered(run_single_simulation, job_args):
            completed += 1

            # Controlla se il job è riuscito o fallito
            if result['success']:
                status = "OK"
            else:
                status = "FAILED"
                failed += 1
                failed_jobs.append(result)

            # Stima tempo rimanente
            elapsed = time.time() - start_time
            avg_time = elapsed / completed
            remaining = avg_time * (n_jobs - completed)

            print_progress_bar(
                completed, n_jobs,
                prefix='  Progresso:',
                suffix=f"ETA: {format_duration(remaining)}    "
            )

    # Fine
    total_time = time.time() - start_time
    print()
    print("-" * 70)
    print()
    print(f"  Simulazioni completate: {completed - failed}/{n_jobs}")
    print(f"  Simulazioni fallite: {failed}")
    print(f"  Tempo totale: {format_duration(total_time)}")

    # Calcola e mostra lo speedup stimato
    if n_jobs > 0 and completed > 0:
        tempo_medio_per_job = total_time / completed
        tempo_sequenziale_stimato = tempo_medio_per_job * n_jobs
        speedup = tempo_sequenziale_stimato / total_time
        print(f"  Speedup stimato: ~{speedup:.1f}x")

    # Report job falliti
    if failed_jobs:
        print()
        print("  Job falliti:")

        # Mostra al massimo 10 job falliti per non riempire lo schermo
        max_jobs_to_show = 10
        jobs_to_display = failed_jobs[:max_jobs_to_show]

        for job in jobs_to_display:
            L_value = job['L']
            T_value = job['T']
            message = job['message'][:50]  # Tronca messaggi lunghi
            print(f"    L={L_value}, T={T_value:.4f}: {message}")

        # Se ci sono altri job falliti, indica quanti
        remaining_failed = len(failed_jobs) - max_jobs_to_show
        if remaining_failed > 0:
            print(f"    ... e altri {remaining_failed} job")

    print()
    print("=" * 70)

    # Exit code basato su successo
    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
