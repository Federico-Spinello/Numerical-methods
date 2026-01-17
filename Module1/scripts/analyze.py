#!/usr/bin/env python3
"""
================================================================================
analyze.py - Analisi dati simulazione Ising 2D
================================================================================
Questo script legge i dati prodotti dalla simulazione Monte Carlo del modello
di Ising 2D e genera i grafici per analizzare:
- Magnetizzazione vs temperatura
- Suscettività magnetica vs temperatura
- Scaling della suscettività con L
- Finite Size Scaling (data collapse)
- Binder cumulant per determinare Tc
- Energia e calore specifico
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import glob
from scipy.optimize import curve_fit
from datetime import datetime

# ============================================================================
# CONFIGURAZIONE
# ============================================================================

# Directory dei dati, grafici e risultati
DATA_DIR = Path("data")
PLOTS_DIR = Path("plots")
PLOTS_DIR.mkdir(exist_ok=True)  # Crea la directory se non esiste

# File dei parametri
PARAMS_FILE = Path("params.txt")

# Parametri fisici noti per Ising 2D
T_C_EXACT = 2.26918  # Temperatura critica esatta (soluzione di Onsager)
GAMMA_NU = 1.75    # Esponente critico γ/ν per Ising 2D
NU = 1.0           # Esponente critico ν per Ising 2D

# Stile dei grafici
plt.style.use('seaborn-v0_8-darkgrid')
plt.rcParams['font.size'] = 10
# Palette di colori estesa (supporta fino a 20 dimensioni L)
COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
          '#1a55FF', '#FF6B35', '#2ECC71', '#E74C3C', '#9B59B6',
          '#34495E', '#F39C12', '#16A085', '#C0392B', '#8E44AD']

# ============================================================================
# UTILITÀ
# ============================================================================

def load_params():
    """
    Carica i parametri dal file params.txt.

    Returns:
        dict: dizionario con i parametri letti dal file
    """
    params = {}
    if not PARAMS_FILE.exists():
        print(f"ATTENZIONE: File {PARAMS_FILE} non trovato")
        return params

    with open(PARAMS_FILE, 'r') as f:
        for line in f:
            line = line.strip()
            # Salta commenti e righe vuote
            if not line or line.startswith('#'):
                continue
            # Parse formato: NOME = valore
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()
                params[key] = value

    return params


def get_color(index):
    """
    Ritorna un colore dalla palette in modo ciclico.
    Se l'indice supera il numero di colori disponibili, ricomincia dall'inizio.

    Args:
        index: indice del colore richiesto

    Returns:
        str: codice colore esadecimale
    """
    return COLORS[index % len(COLORS)]


def bootstrap_peak_error(T_data, O_data, err_O_data, n_bootstrap=1000, window_size=5):
    """
    Stima l'errore sulla temperatura del picco usando il metodo bootstrap non parametrico.

    Algoritmo (da Appendice 6.2 del paper):
    1. Identifica il picco: T_max = T[i_max] dove i_max = argmax O(T)
    2. Estrae una finestra di (2*window_size + 1) punti attorno al picco
    3. Esegue n_bootstrap iterazioni:
       - Genera un campione bootstrap aggiungendo rumore gaussiano N(0, sigma_O(T_i))
       - Trova il massimo del campione bootstrap
    4. L'errore è la deviazione standard della distribuzione bootstrap

    Args:
        T_data: array delle temperature
        O_data: array dell'osservabile (chi o C)
        err_O_data: array degli errori sull'osservabile
        n_bootstrap: numero di iterazioni bootstrap (default 1000)
        window_size: numero di punti da considerare a sinistra e destra del picco (default 5)

    Returns:
        T_max: temperatura del picco
        T_max_err: errore bootstrap sulla temperatura del picco
    """
    # Step 1: Trova l'indice del picco
    idx_max = np.argmax(O_data)
    T_max_original = T_data[idx_max]

    # Step 2: Estrai finestra attorno al picco (±window_size punti)
    idx_start = max(0, idx_max - window_size)
    idx_end = min(len(T_data), idx_max + window_size + 1)

    T_window = T_data[idx_start:idx_end]
    O_window = O_data[idx_start:idx_end]
    err_window = err_O_data[idx_start:idx_end]

    # Step 3: Bootstrap
    T_max_bootstrap = []

    for r in range(n_bootstrap):
        # Genera campione bootstrap aggiungendo rumore gaussiano N(0, sigma_O)
        noise = np.random.normal(0, err_window)
        O_bootstrap = O_window + noise

        # Trova il massimo del campione bootstrap
        idx_max_boot = np.argmax(O_bootstrap)
        T_max_boot = T_window[idx_max_boot]
        T_max_bootstrap.append(T_max_boot)

    T_max_bootstrap = np.array(T_max_bootstrap)

    # Step 4: L'errore è la deviazione standard della distribuzione bootstrap
    T_max_err = np.std(T_max_bootstrap)

    return T_max_original, T_max_err


def save_results(results_dict):
    """
    Salva i risultati dell'analisi in formato testo (risultati.txt).

    Args:
        results_dict: dizionario con tutti i risultati numerici
    """
    # Aggiungi timestamp
    results_dict['date'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # Salva solo versione testuale leggibile
    txt_file = Path("risultati.txt")
    with open(txt_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("RISULTATI ANALISI ISING 2D - ALGORITMO DI WOLFF\n")
        f.write("="*70 + "\n")
        f.write(f"Data analisi: {results_dict['date']}\n")
        f.write(f"Dimensioni L simulate: {results_dict['L_values']}\n")
        f.write(f"Numero di temperature: {results_dict['n_temps']}\n")
        f.write("\n")

        f.write("-"*70 + "\n")
        f.write("TEMPERATURA CRITICA\n")
        f.write("-"*70 + "\n")
        f.write(f"Tc teorico (Onsager):          {results_dict['Tc_theory']:.5f}\n")
        f.write(f"Tc da picchi chi(L):           {results_dict['Tc_peaks']:.4f} +/- {results_dict['Tc_peaks_err']:.4f}\n")
        f.write(f"  chi^2_red (FSS fit):         {results_dict.get('chi2_red_chi', 0.0):.4f}\n")
        f.write(f"Tc da picchi C(L):             {results_dict['Tc_heat']:.4f} +/- {results_dict['Tc_heat_err']:.4f}\n")
        f.write(f"  chi^2_red (FSS fit):         {results_dict.get('chi2_red_C', 0.0):.4f}\n")
        f.write(f"Tc da Binder crossing:         {results_dict['Tc_binder']:.4f} +/- {results_dict['Tc_binder_err']:.4f}\n")
        f.write("\n")

        f.write("-"*70 + "\n")
        f.write("ESPONENTI CRITICI\n")
        f.write("-"*70 + "\n")
        f.write(f"gamma/nu misurato:                {results_dict['gamma_nu']:.4f} +/- {results_dict['gamma_nu_err']:.4f}\n")
        f.write(f"gamma/nu teorico:                 {results_dict['gamma_nu_theory']:.4f}\n")
        f.write(f"chi^2_red (collapse chi/L^gamma/nu):   {results_dict.get('chi2_red_scaling', results_dict['gamma_nu_chi2']):.4f}\n")

    print(f"  Risultati salvati in: {txt_file}")

# ============================================================================
# CARICAMENTO DATI
# ============================================================================

def load_all_data():
    """
    Carica tutti i file dati dalla directory DATA_DIR e li organizza per L.

    Ritorna:
        dict: dizionario con chiavi L, contenente array numpy per:
              T, m, E, chi, C, binder
    """
    data = {}

    # Trova tutti i file .dat nella directory data/
    files = sorted(glob.glob(str(DATA_DIR / "L*.dat")))

    if not files:
        print(f"ERRORE: Nessun file trovato in {DATA_DIR}/")
        return None

    # Carica ogni file
    for filepath in files:
        # Parse del nome file: L20_T2.1000.dat
        filename = Path(filepath).stem
        parts = filename.split('_')
        L = int(parts[0][1:])  # Rimuovi la 'L' e converti a int

        # Inizializza dizionario per questo L se non esiste
        if L not in data:
            data[L] = {'T': [], 'm': [], 'err_m': [], 'E': [], 'err_E': [],
                       'chi': [], 'err_chi': [], 'C': [], 'err_C': [],
                       'binder': [], 'err_binder': []}

        # Leggi i dati dal file (formato: T m err_m E err_E chi err_chi C err_C binder err_binder)
        values = np.loadtxt(filepath)
        # Assicurati che values sia un array 1D (una sola riga)
        if values.ndim == 1:
            data[L]['T'].append(values[0])
            data[L]['m'].append(values[1])
            data[L]['err_m'].append(values[2])
            data[L]['E'].append(values[3])
            data[L]['err_E'].append(values[4])
            data[L]['chi'].append(values[5])
            data[L]['err_chi'].append(values[6])
            data[L]['C'].append(values[7])
            data[L]['err_C'].append(values[8])
            data[L]['binder'].append(values[9])
            data[L]['err_binder'].append(values[10])
            

    # Converti liste in array numpy e ordina per temperatura
    for L in data:
        for key in data[L]:
            data[L][key] = np.array(data[L][key])

        # Ordina tutti gli array per temperatura crescente
        idx = np.argsort(data[L]['T'])
        for key in data[L]:
            data[L][key] = data[L][key][idx]

    return data

# ============================================================================
# GRAFICO 1: Magnetizzazione vs Temperatura
# ============================================================================

def plot_magnetization(data):
    """
    Grafico della magnetizzazione <|m|> in funzione della temperatura
    per diverse dimensioni del reticolo.
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    L_values = sorted(data.keys())

    # Plot per ogni dimensione L
    for i, L in enumerate(L_values):
        color = get_color(i)
        # plot
        ax.plot(data[L]['T'], data[L]['m'],
                '-', label=f'L={L}', color=color,
                markersize=5, linewidth=1.5, zorder=3)


    # Linea verticale alla temperatura critica
    ax.axvline(T_C_EXACT, color='red', linestyle='--',
              alpha=0.6, linewidth=2, label=f'$T_c$ = {T_C_EXACT}')

    ax.set_xlabel('Temperatura $T$ [$J/k_B$]', fontsize=13)
    ax.set_ylabel('Magnetizzazione $\\langle |m| \\rangle$', fontsize=13)
    ax.set_title('Magnetizzazione vs Temperatura (Algoritmo di Wolff)', fontsize=14, pad=15)
    ax.legend(fontsize=11, framealpha=0.9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'magnetization.png', dpi=300, bbox_inches='tight')
    print(" Grafico magnetization.png salvato")
    plt.close()

# ============================================================================
# GRAFICO 2: Suscettività vs Temperatura
# ============================================================================

def plot_susceptibility(data):
    """
    Grafico della suscettività magnetica χ in funzione della temperatura.
    Il picco della suscettività indica la temperatura critica.
    Usa FSS per estrarre Tc: T_max(L) ≈ Tc + a*L^(-1/ν)
    Gli errori su T_max(L) sono stimati con bootstrap non parametrico (1000 iterazioni).
    Ritorna i valori dei picchi per uso successivo.
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    L_values = sorted(data.keys())
    T_max_values = []
    T_max_errors = []
    L_array = []

    # Plot per ogni dimensione L e trova i picchi con errori bootstrap
    for i, L in enumerate(L_values):
        color = get_color(i)
        # Linea principale
        ax.plot(data[L]['T'], data[L]['chi'],
                '-', label=f'L={L}', color=color,
                markersize=5, linewidth=1.5, zorder=3)
        # Banda di errore
        ax.fill_between(data[L]['T'],
                        data[L]['chi'] - data[L]['err_chi'],
                        data[L]['chi'] + data[L]['err_chi'],
                        color=color, alpha=0.2, zorder=1)

        # Trova temperatura del picco con errore bootstrap per questo L
        T_max, T_max_err = bootstrap_peak_error(
            data[L]['T'],
            data[L]['chi'],
            data[L]['err_chi'],
            n_bootstrap=1000,
            window_size=5
        )
        T_max_values.append(T_max)
        T_max_errors.append(T_max_err)
        L_array.append(L)

    # FIT FSS: T_max(L) = Tc + a * L^(-1/nu)
    # Con nu = 1 per Ising 2D: T_max(L) = Tc + a/L
    L_array = np.array(L_array)
    T_max_values = np.array(T_max_values)
    T_max_errors = np.array(T_max_errors)

    def fss_shift(L, Tc, a):
        """T_max(L) = Tc + a * L^(-1/nu) con nu=1"""
        return Tc + a / L

    # Fit pesato con gli errori bootstrap
    # sigma=T_max_errors fornisce i pesi 1/sigma^2 per il fit
    # absolute_sigma=True indica che gli errori sono assoluti (non relativi)
    popt, pcov = curve_fit(
        fss_shift,
        L_array,
        T_max_values,
        p0=[T_C_EXACT, 0.5],
        sigma=T_max_errors,
        absolute_sigma=True
    )
    Tc_from_peaks, a_fit = popt
    perr = np.sqrt(np.diag(pcov))
    Tc_from_peaks_err = perr[0]

    # Calcolo del chi^2 ridotto per il fit FSS dei picchi chi
    # Usando gli errori bootstrap come sigma_i
    T_max_predicted = fss_shift(L_array, Tc_from_peaks, a_fit)
    residuals_chi = T_max_values - T_max_predicted

    # Numero di parametri del fit (Tc e a)
    n_params_chi = 2

    # Gradi di liberta: numero di punti - numero di parametri
    dof_chi = len(T_max_values) - n_params_chi

    # Calcolo del chi^2 ridotto usando gli errori bootstrap
    # chi^2 = sum((y_obs - y_fit)^2 / sigma_i^2)
    if dof_chi > 0:
        chi2_chi = np.sum((residuals_chi / T_max_errors) ** 2)
        chi2_red_chi = chi2_chi / dof_chi
    else:
        chi2_red_chi = 0.0

    print(f"  -> Fit FSS picchi chi: Tc = {Tc_from_peaks:.4f} +/- {Tc_from_peaks_err:.4f}")
    print(f"                         chi^2_red = {chi2_red_chi:.4f}")
    print(f"                         a = {a_fit:.4f} +/- {perr[1]:.4f}")

    # Linea verticale alla temperatura critica teorica
    ax.axvline(T_C_EXACT, color='red', linestyle='--',
              alpha=0.4, linewidth=2, label=f'$T_c$ teorico = {T_C_EXACT}')

    # Linea verticale alla temperatura critica stimata dai picchi
    ax.axvline(Tc_from_peaks, color='blue', linestyle='-',
              alpha=0.6, linewidth=2, label=f'$T_c$ (FSS picchi χ) = {Tc_from_peaks:.3f}')

    ax.set_xlabel('Temperatura $T$ [$J/k_B$]', fontsize=13)
    ax.set_ylabel('Suscettività $\\chi$', fontsize=13)
    ax.set_title('Suscettività Magnetica vs Temperatura', fontsize=14, pad=15)
    ax.legend(fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'susceptibility.png', dpi=300, bbox_inches='tight')
    print(" Grafico susceptibility.png salvato")
    plt.close()

    return Tc_from_peaks, Tc_from_peaks_err, chi2_red_chi

# ============================================================================
# GRAFICO 3: Scaling del picco di suscettività
# ============================================================================

def plot_chi_scaling(data):
    """
    Verifica la legge di scaling: chi_max ~ L^(gamma/nu)
    Calcola chi^2 ridotto per la bonta del fit.
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    L_values = sorted(data.keys())
    chi_max = []
    T_max = []

    # Trova il picco di chi per ogni L
    for L in L_values:
        idx_max = np.argmax(data[L]['chi'])
        chi_max.append(data[L]['chi'][idx_max])
        T_max.append(data[L]['T'][idx_max])

    L_values = np.array(L_values)
    chi_max = np.array(chi_max)

    # Funzione per il fit: chi_max = A * L^(gamma/nu)
    def power_law(L, A, gamma_nu):
        return A * L**gamma_nu

    # Fit con curve_fit (senza pesi)
    popt, pcov = curve_fit(power_law, L_values, chi_max, p0=[1.0, GAMMA_NU])
    A_fit, exponent = popt
    perr = np.sqrt(np.diag(pcov))
    A_err, exponent_err = perr

    # Calcola chi^2 ridotto
    # Metodo: stima sigma dai residui (fit non pesato)
    residuals = chi_max - power_law(L_values, *popt)
    n_params = 2  # A e gamma/nu
    dof = len(chi_max) - n_params  # gradi di liberta

    # Stima sigma dai residui
    sigma = np.std(residuals)

    # Calcolo chi^2 ridotto
    chi2_gamma = np.sum((residuals / sigma)**2)
    chi2_red_gamma = chi2_gamma / dof

    # Plot dei dati
    ax.plot(L_values, chi_max, 'o', markersize=10, label='Dati', color='#1f77b4')

    # Plot del fit
    L_fit = np.linspace(L_values[0]*0.9, L_values[-1]*1.1, 100)
    chi_fit = power_law(L_fit, A_fit, exponent)
    ax.plot(L_fit, chi_fit, '--', linewidth=2,
           label=f'Fit: $\\chi_{{max}} = {A_fit:.2f} \\cdot L^{{{exponent:.3f}}}$ ($\\chi^2_{{red}}={chi2_red_gamma:.4f}$)',
           color='#ff7f0e')

    # Curva teorica
    chi_theo = A_fit * L_fit**GAMMA_NU
    ax.plot(L_fit, chi_theo, ':', linewidth=2,
           label=f'Teoria: $\\gamma/\\nu = {GAMMA_NU}$',
           color='red', alpha=0.7)

    ax.set_xlabel('Dimensione reticolo $L$', fontsize=13)
    ax.set_ylabel('$\\chi_{max}$', fontsize=13)
    ax.set_title('Scaling del Picco di Suscettivita', fontsize=14, pad=15)
    ax.legend(fontsize=10, framealpha=0.9)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'chi_scaling.png', dpi=300, bbox_inches='tight')
    print("Grafico chi_scaling.png salvato")
    print(f"  -> gamma/nu = {exponent:.4f} +/- {exponent_err:.4f} (teorico: {GAMMA_NU})")
    print(f"  -> A = {A_fit:.2f} +/- {A_err:.2f}")
    print(f"  -> chi^2_red = {chi2_red_gamma:.4f}")
    plt.close()

    return exponent, T_max, exponent_err, chi2_red_gamma

# ============================================================================
# GRAFICO 4: Finite Size Scaling (Data Collapse)
# ============================================================================

def plot_fss(data):
    """
    Data collapse usando Finite Size Scaling.
    Se la temperatura critica Tc è corretta, i dati di tutte le dimensioni
    dovrebbero collassare su una curva universale quando plottati come:
        y = χ / L^(γ/ν)  vs  x = L^(1/ν) * (T - Tc)
    """
    
    Tc = T_C_EXACT

    fig, ax = plt.subplots(figsize=(9, 6))

    L_values = sorted(data.keys())

    # Plot per ogni L con variabili scalate
    for i, L in enumerate(L_values):
        T = data[L]['T']
        chi = data[L]['chi']
        err_chi = data[L]['err_chi']

        # Variabili scalate secondo Finite Size Scaling
        x_scaled = L**(1/NU) * (T - Tc)  # Distanza ridotta da Tc
        y_scaled = chi / L**GAMMA_NU      # Suscettività ridotta

        color = get_color(i)
        # Grafico
        ax.plot(x_scaled, y_scaled, 'o',
                label=f'L={L}', color=color,
                markersize=6, zorder=3)

    # Linea verticale a x=0 (che corrisponde a T=Tc)
    ax.axvline(0, color='red', linestyle='--', alpha=0.3, linewidth=2)
    ax.set_xlabel('$L^{1/\\nu}(T - T_c)$', fontsize=14)
    ax.set_ylabel('$\\chi / L^{\\gamma/\\nu}$', fontsize=14)
    ax.set_title(f'Data Collapse (Finite Size Scaling, $T_c = {Tc:.3f}$)',
                fontsize=14, pad=15)
    ax.legend(fontsize=11, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-8, 8)

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'fss_collapse.png', dpi=300, bbox_inches='tight')
    print(" Grafico fss_collapse.png salvato")
    plt.close()

# ============================================================================
# GRAFICO 5: Binder Cumulant
# ============================================================================

def calculate_binder_crossing(data):
    """
    Calcola Tc dal crossing del Binder cumulant usando la dispersione minima.

    METODO:
    1. Restringe il calcolo al range T ∈ [2.2, 2.3] vicino a Tc
    2. Crea una griglia fine di temperature in questo range
    3. Per ogni T, interpola U_L(T) per tutti gli L
    4. Calcola σ(T) = dispersione dei valori U_L tra tutti gli L
    5. Trova T dove σ(T) è minimo -> questo è il crossing point
    6. Stima l'errore dal range 

    Ritorna:
        float: Stima di Tc dal crossing
        float: Errore stimato dalla larghezza del minimo
    """
    L_values = sorted(data.keys())

    if len(L_values) < 2:
        return T_C_EXACT, 0.0

    # -------------------------------------------------------------------------
    # STEP 1: Definisce il range di temperatura ristretto per il crossing
    # -------------------------------------------------------------------------
    # Restringe l'analisi a [2.2, 2.3], regione vicina a Tc_teorico = 2.269
    T_min_crossing = 2.2
    T_max_crossing = 2.3

    # Determina il range di temperature disponibile

    T_min = data[40]['T'].min()
    T_max = data[40]['T'].max()

    # Usa l'intersezione: il range richiesto [2.2, 2.3] intersecato con
    # il range disponibile nei dati
    T_min = max(T_min, T_min_crossing)
    T_max = min(T_max, T_max_crossing)

    # -------------------------------------------------------------------------
    # STEP 2: Crea griglia fine di temperature nel range ristretto
    # -------------------------------------------------------------------------
    # 50000 punti -> risoluzione ΔT ≈ 0.000002 nel range [2.2, 2.3]
    T_range = np.linspace(T_min, T_max, 50000)

    # -------------------------------------------------------------------------
    # STEP 3: Interpola U_L(T) su griglia comune per tutti gli L
    # -------------------------------------------------------------------------
    # Crea matrice: righe = temperature, colonne = dimensioni L
    n_temps = len(T_range)
    n_L = len(L_values)
    U_matrix = np.zeros((n_temps, n_L))

    # Riempi la matrice colonna per colonna (una colonna per ogni L)
    for j in range(n_L):
        L = L_values[j]
        T_data = data[L]['T']
        U_data = data[L]['binder']
        # Interpola tutti i punti T_range per questo L
        for i in range(n_temps):
            U_matrix[i, j] = np.interp(T_range[i], T_data, U_data)

    # -------------------------------------------------------------------------
    # STEP 4: Calcola la dispersione sigma(T) tra i diversi L
    # -------------------------------------------------------------------------
    # Per ogni T, calcola la deviazione standard dei valori U_L
    # sigma(T) = std({U_L1(T), U_L2(T), U_L3(T), ...})
    # Al crossing tutte le curve convergono -> sigma(Tc) e' minimo
    std_devs = np.std(U_matrix, axis=1)

    # -------------------------------------------------------------------------
    # STEP 5: Trova il minimo della dispersione -> Tc
    # -------------------------------------------------------------------------
    idx_min = np.argmin(std_devs)
    Tc_binder = T_range[idx_min]
    std_min = std_devs[idx_min]

    # -------------------------------------------------------------------------
    # STEP 6: Stima l'errore dalla larghezza del minimo
    # -------------------------------------------------------------------------
    # L'errore viene definito come la semi-ampiezza del range di temperature
    # dove σ(T) ≤ 1.05·σ_min (entro 5% del minimo).
    # Questo rappresenta l'incertezza sulla posizione del crossing dovuta
    # alla larghezza finita del minimo e alle fluttuazioni statistiche.

    # Definisco la soglia di accettazione attorno al minimo
    threshold = 1.05 * std_min

    # Lista (temporanea) dei valori di T che soddisfano il criterio
    T_near_min = []

    # Scorro esplicitamente tutti i punti
    for i in range(len(std_devs)):
        current_std = std_devs[i]
        current_T = T_range[i]

        # Accetto il punto se la deviazione standard
        # non supera la soglia fissata
        if current_std <= threshold:
            T_near_min.append(current_T)

    # Converto in array NumPy dopo aver selezionato i punti
    T_near_min = np.array(T_near_min)

    # Ora decido come stimare l'errore
    if len(T_near_min) > 1:
        # Errore = semi-ampiezza dell'intervallo compatibile col minimo
        T_max_near = np.max(T_near_min)
        T_min_near = np.min(T_near_min)
        T_error = 0.5 * (T_max_near - T_min_near)
    else:
        # Minimo troppo stretto: uso la risoluzione della griglia
        grid_spacing = (T_max - T_min) / len(T_range)
        T_error = grid_spacing


    # Ritorna anche T_range e std_devs per il plot
    return Tc_binder, T_error, T_range, std_devs


def plot_binder_dispersion(T_range, std_devs):
    """
    Grafico della dispersione sigma(T) del Binder cumulant.
    Usa i dati gia calcolati da calculate_binder_crossing.
    """
    # Trova il minimo
    idx_min = np.argmin(std_devs)
    Tc_crossing = T_range[idx_min]
    sigma_min = std_devs[idx_min]

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(T_range, std_devs, 'b-', linewidth=2, label='Dispersione $U_L(T)$', zorder=3)
    ax.axvline(Tc_crossing, color='red', linestyle='--', linewidth=2,
               label=f'$T_c$ (minimo) = {Tc_crossing:.4f}', zorder=2)
    ax.axvline(T_C_EXACT, color='darkgreen', linestyle=':', linewidth=2,
               label=f'$T_c$ teorico = {T_C_EXACT}', zorder=2)

    # Evidenzia il minimo
    ax.plot(Tc_crossing, sigma_min, 'ro', markersize=10,
            label='Minimo dispersione', zorder=4)

    ax.set_xlabel('Temperatura $T$ [$J/k_B$]', fontsize=13)
    ax.set_ylabel('Dispersione $\\sigma(U_L)$', fontsize=13)
    ax.set_title('Dispersione del Binder Cumulant (metodo del crossing)', fontsize=14, pad=15)
    ax.legend(fontsize=11, framealpha=0.9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'binder_dispersion.png', dpi=300, bbox_inches='tight')
    print("Grafico binder_dispersion.png salvato")
    plt.close()

def plot_binder(data, Tc_binder, Tc_error):
    """
    Binder cumulant: U = 1 - <m^4>/(3*<m^2>^2)
    Le curve per diversi L si intersecano alla temperatura critica Tc.
    Usa Tc_binder e Tc_error gia calcolati da calculate_binder_crossing.
    """

    fig, ax = plt.subplots(figsize=(8, 6))

    L_values = sorted(data.keys())

    # Plot per ogni L
    for i, L in enumerate(L_values):
        color = get_color(i)
        ax.plot(data[L]['T'], data[L]['binder'],
               '-', label=f'L={L}', color=get_color(i),
               markersize=5, linewidth=1.5)
        ax.fill_between(data[L]['T'],
                            data[L]['binder'] - data[L]['err_binder'],
                            data[L]['binder'] + data[L]['err_binder'],
                            color=color, alpha=0.2, zorder=1)

    # Linee di riferimento
    ax.axvline(T_C_EXACT, color='red', linestyle='--',
              alpha=0.4, linewidth=2, label=f'$T_c$ teorico = {T_C_EXACT}')
    ax.axvline(Tc_binder, color='blue', linestyle='-',
              alpha=0.6, linewidth=2, label=f'$T_c$ (crossing) = {Tc_binder:.3f}')

    ax.set_xlabel('Temperatura $T$ [$J/k_B$]', fontsize=13)
    ax.set_ylabel('Binder Cumulant $U_L$', fontsize=13)
    ax.set_title('Binder Cumulant (crossing point -> $T_c$)', fontsize=14, pad=15)
    ax.legend(fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 3.5)

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'binder_cumulant.png', dpi=300, bbox_inches='tight')
    print("Grafico binder_cumulant.png salvato")
    print(f"  -> Tc dal crossing Binder: {Tc_binder:.4f} +/- {Tc_error:.4f}")
    plt.close()

# ============================================================================
# GRAFICO 6: Energia e Calore Specifico
# ============================================================================

def plot_energy_heat(data):
    """
    Grafici di energia e calore specifico in funzione della temperatura.
    Usa FSS per estrarre Tc dai picchi: T_max(L) ≈ Tc + b*L^(-1/ν)
    Gli errori su T_max(L) sono stimati con bootstrap non parametrico (1000 iterazioni).
    Ritorna Tc dal calore specifico e i dettagli dei picchi per ogni L.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    L_values = sorted(data.keys())
    C_max_values = []
    T_at_C_max = []
    T_at_C_max_errors = []
    L_array = []
    heat_peaks_details = {}

    # --- Grafico Energia ---
    for i, L in enumerate(L_values):
        color = get_color(i)
        # Linea principale
        ax1.plot(data[L]['T'], data[L]['E'],
                 '-', label=f'L={L}', color=color,
                 markersize=5, linewidth=1.5, zorder=3)
        # Banda di errore
        ax1.fill_between(data[L]['T'],
                         data[L]['E'] - data[L]['err_E'],
                         data[L]['E'] + data[L]['err_E'],
                         color=color, alpha=0.2, zorder=1)

    ax1.axvline(T_C_EXACT, color='red', linestyle='--', alpha=0.6, linewidth=2)
    ax1.set_xlabel('Temperatura $T$ [$J/k_B$]', fontsize=12)
    ax1.set_ylabel('Energia per spin $E/N$', fontsize=12)
    ax1.set_title('Energia vs Temperatura', fontsize=13)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    # --- Grafico Calore Specifico ---
    for i, L in enumerate(L_values):
        color = get_color(i)
        # Linea principale
        ax2.plot(data[L]['T'], data[L]['C'],
                 '-', label=f'L={L}', color=color,
                 markersize=5, linewidth=1.5, zorder=3)

        # Trova il picco del calore specifico con errore bootstrap per questo L
        T_max, T_max_err = bootstrap_peak_error(
            data[L]['T'],
            data[L]['C'],
            data[L]['err_C'],
            n_bootstrap=1000,
            window_size=5
        )

        # Trova il valore di C_max corrispondente
        idx_max = np.argmax(data[L]['C'])
        C_max = data[L]['C'][idx_max]

        C_max_values.append(C_max)
        T_at_C_max.append(T_max)
        T_at_C_max_errors.append(T_max_err)
        L_array.append(L)

        # Salva i dettagli per ogni L
        heat_peaks_details[L] = {
            'T_max': float(T_max),
            'T_max_err': float(T_max_err),
            'C_max': float(C_max)
        }

    # FIT FSS: T_max(L) = Tc + b * L^(-1/nu)
    # Con nu = 1 per Ising 2D: T_max(L) = Tc + b/L
    L_array = np.array(L_array)
    T_at_C_max = np.array(T_at_C_max)
    T_at_C_max_errors = np.array(T_at_C_max_errors)

    def fss_shift(L, Tc, b):
        return Tc + b / L

    # Fit pesato con gli errori bootstrap
    # sigma=T_at_C_max_errors fornisce i pesi 1/sigma^2 per il fit
    # absolute_sigma=True indica che gli errori sono assoluti (non relativi)
    popt, pcov = curve_fit(
        fss_shift,
        L_array,
        T_at_C_max,
        p0=[T_C_EXACT, 0.5],
        sigma=T_at_C_max_errors,
        absolute_sigma=True
    )
    Tc_from_heat, b_fit = popt
    perr = np.sqrt(np.diag(pcov))
    Tc_from_heat_err = perr[0]

    # Calcolo del chi^2 ridotto per il fit FSS dei picchi del calore specifico
    # Usando gli errori bootstrap come sigma_i
    T_at_C_max_predicted = fss_shift(L_array, Tc_from_heat, b_fit)
    residuals_C = T_at_C_max - T_at_C_max_predicted

    # Numero di parametri del fit (Tc e b)
    n_params_C = 2

    # Gradi di liberta: numero di punti - numero di parametri
    dof_C = len(T_at_C_max) - n_params_C

    # Calcolo del chi^2 ridotto usando gli errori bootstrap
    # chi^2 = sum((y_obs - y_fit)^2 / sigma_i^2)
    if dof_C > 0:
        chi2_C = np.sum((residuals_C / T_at_C_max_errors) ** 2)
        chi2_red_C = chi2_C / dof_C
    else:
        chi2_red_C = 0.0

    print(f"  -> Fit FSS picchi C: Tc = {Tc_from_heat:.4f} +/- {Tc_from_heat_err:.4f}")
    print(f"                       chi^2_red = {chi2_red_C:.4f}")
    print(f"                       b = {b_fit:.4f} +/- {perr[1]:.4f}")

    # Aggiungi linea verticale per Tc teorico
    ax2.axvline(T_C_EXACT, color='red', linestyle='--', alpha=0.4, linewidth=2,
                label=f'$T_c$ teorico = {T_C_EXACT}')

    # Aggiungi linea verticale per Tc da picchi calore specifico
    ax2.axvline(Tc_from_heat, color='green', linestyle='-', alpha=0.6, linewidth=2,
                label=f'$T_c$ (FSS picchi C) = {Tc_from_heat:.3f}')

    ax2.set_xlabel('Temperatura $T$ [$J/k_B$]', fontsize=12)
    ax2.set_ylabel('Calore Specifico $C$', fontsize=12)
    ax2.set_title('Calore Specifico vs Temperatura', fontsize=13)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'energy_heat.png', dpi=300, bbox_inches='tight')
    print(" Grafico energy_heat.png salvato")
    plt.close()

    return Tc_from_heat, Tc_from_heat_err, chi2_red_C, heat_peaks_details

# ============================================================================
# MAIN
# ============================================================================

def main():
    """
    Funzione principale: carica i dati e genera tutti i grafici.
    """
    print("="*70)
    print("  ANALISI ISING 2D - ALGORITMO DI WOLFF")
    print("="*70)
    print(f"Directory dati: {DATA_DIR}")
    print(f"Directory plot: {PLOTS_DIR}\n")

    # Carica i dati
    print("Caricamento dati...")
    data = load_all_data()

    if data is None or len(data) == 0:
        print("ERRORE: Nessun dato trovato!")
        return

    # Informazioni sui dati caricati
    L_values = sorted(data.keys())
    n_temps = len(data[L_values[0]]['T'])
    print(f" Dati caricati: {len(L_values)} dimensioni, {n_temps} temperature ciascuna")
    print(f"  L = {L_values}\n")

    # Genera tutti i grafici
    print("Generazione grafici...")
    print("-" * 70)

    # Calcola Tc dal Binder crossing
    Tc_binder, Tc_error, T_range_binder, std_devs_binder = calculate_binder_crossing(data)

    # Genera i grafici
    plot_magnetization(data)
    Tc_peaks_chi, Tc_peaks_chi_err, chi2_red_chi = plot_susceptibility(data)
    exponent, _, exponent_err, chi2_red_gamma = plot_chi_scaling(data)
    plot_fss(data)
    plot_binder_dispersion(T_range_binder, std_devs_binder)
    plot_binder(data, Tc_binder, Tc_error)
    Tc_from_heat, Tc_from_heat_err, chi2_red_C, heat_peaks_details = plot_energy_heat(data)
    plot_heat_with_errors(data)

    print("-" * 70)

    # Riepilogo risultati
    print("\n" + "="*70)
    print("  RISULTATI")
    print("="*70)
    print(f"Esponente critico γ/ν:")
    print(f"  Misurato:  {exponent:.3f}")
    print(f"  Teorico:   {GAMMA_NU}")
    print(f"  Differenza: {abs(exponent - GAMMA_NU)/GAMMA_NU * 100:.1f}%")
    print(f"\nTemperatura critica Tc:")
    print(f"  Teorico (Onsager):     {T_C_EXACT}")
    print(f"  Dai picchi χ(L):       {Tc_peaks_chi:.4f} +/- {Tc_peaks_chi_err:.4f}")
    print(f"  Dai picchi C(L):       {Tc_from_heat:.4f} +/- {Tc_from_heat_err:.4f}")
    print(f"  Dal crossing Binder:   {Tc_binder:.4f} +/- {Tc_error:.4f}")
    print(f"\n  -> Miglior stima Tc = {Tc_binder:.4f} (dal Binder crossing)")
    print(f"  -> Errore rispetto teorico: {abs(Tc_binder - T_C_EXACT)*1000:.1f} mK")
    print(f"\nPicchi del calore specifico per ogni L:")
    for L in sorted(heat_peaks_details.keys()):
        peak_data = heat_peaks_details[L]
        print(f"  L={L:3d}:  T_max = {peak_data['T_max']:.4f},  C_max = {peak_data['C_max']:.6f}")
    print("\n Analisi completata!")
    print(f"Grafici salvati in: {PLOTS_DIR}/\n")

    # Prepara dizionario risultati per salvare
    results = {
        'L_values': L_values,
        'n_temps': n_temps,
        'Tc_theory': T_C_EXACT,
        'Tc_peaks': float(Tc_peaks_chi),
        'Tc_peaks_err': float(Tc_peaks_chi_err),
        'chi2_red_chi': float(chi2_red_chi),
        'Tc_heat': float(Tc_from_heat),
        'Tc_heat_err': float(Tc_from_heat_err),
        'chi2_red_C': float(chi2_red_C),
        'Tc_binder': float(Tc_binder),
        'Tc_binder_err': float(Tc_error),
        'Tc_binder_error_mK': float(abs(Tc_binder - T_C_EXACT) * 1000),
        'Tc_binder_error_pct': float(abs(Tc_binder - T_C_EXACT) / T_C_EXACT * 100),
        'gamma_nu': float(exponent),
        'gamma_nu_err': float(exponent_err),
        'gamma_nu_theory': float(GAMMA_NU),
        'gamma_nu_chi2': float(chi2_red_gamma),
        'U_star_theory': 1.6,
        'heat_peaks': heat_peaks_details
    }

    # Salva risultati in file
    save_results(results)

def plot_thermalization():
    """
    Genera grafico della termalizzazione e zoom.

    Legge dati reali da data/thermalization_data.dat se disponibili
    (generati con 'make thermalization'), altrimenti usa dati simulati.
    """

    therm_data_file = DATA_DIR / "thermalization_data.dat"

    # Leggi parametri da params.txt
    params = load_params()
    L_value = int(params.get('THERM_L', 100))
    T_value = float(params.get('THERM_T', 2.269))
    thermalization_steps = int(params.get('THERMALIZATION', 10000))

    print(f"  -> Parametri da params.txt: L={L_value}, T={T_value}, THERMALIZATION={thermalization_steps}")

    # Controlla se ci sono dati reali
    if therm_data_file.exists():
        print("  -> Usando dati reali da thermalization_data.dat")
        try:
            # Carica dati: step, E, m
            data = np.loadtxt(therm_data_file, comments='#')
            steps = data[:, 0].astype(int)
            E = data[:, 1]
            m = data[:, 2]
            n_steps = len(steps)
            is_real_data = True
        except Exception as e:
            print(f"  ! Errore lettura dati reali: {e}")
            print("  -> Usando dati simulati")
            is_real_data = False
    else:
        print("  -> Dati reali non trovati (esegui 'make thermalization' per generarli)")
        print("  -> Usando dati simulati basati su comportamento tipico")
        is_real_data = False

    # Se non ci sono dati reali, genera dati simulati
    if not is_real_data:
        # Genera dati simulati per dimostrazione
        n_steps = 410000  # Simula dati completi
        steps = np.arange(n_steps)

        # Energia: parte da ~-1.6, raggiunge equilibrio ~-1.205 dopo termalizzazione
        E = -1.6 + 0.4 * (1 - np.exp(-steps/500))
        E += 0.013 * np.random.randn(n_steps)  # aggiungi rumore

        # Magnetizzazione: parte alta ~0.8, scende a ~0.05 (vicino a Tc)
        m = 0.8 * np.exp(-steps/200) + 0.05
        m += 0.038 * np.random.randn(n_steps)  # aggiungi rumore
        m = np.abs(m)  # assicurati sia positiva

    # Calcola statistiche sulla fase di equilibrio
    E_eq = E[thermalization_steps:]
    m_eq = m[thermalization_steps:]
    E_mean = np.mean(E_eq)
    E_std = np.std(E_eq)
    m_mean = np.mean(m_eq)
    m_std = np.std(m_eq)

    # ====== GRAFICO COMPLETO ======
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # --- Energia ---
    ax1.axvspan(0, thermalization_steps, alpha=0.3, color='red',
                label=f'Termalizzazione ({thermalization_steps} steps)', zorder=1)
    ax1.axvspan(thermalization_steps, n_steps, alpha=0.3, color='green',
                label='Misure (equilibrio)', zorder=1)
    ax1.plot(steps, E, 'b-', linewidth=0.5, zorder=2)
    ax1.axhspan(E_mean - E_std, E_mean + E_std,
                xmin=(thermalization_steps/n_steps), xmax=1.0,
                color='gray', alpha=0.3, label=r'$\pm 1\sigma$', zorder=3)
    ax1.axhline(E_mean, color='black', linestyle='--', linewidth=1.5,
                label=f'Media equilibrio: {E_mean:.3f} $\\pm$ {E_std:.3f}', zorder=4)
    ax1.axvline(1000, color='magenta', linestyle=':', linewidth=1.5, label='1000 steps', zorder=5)

    ax1.set_xlabel('Cluster updates', fontsize=12)
    ax1.set_ylabel('Energia per spin $E/N$', fontsize=12)
    ax1.set_title(f'Termalizzazione: Evoluzione dell\'Energia ($L = {L_value}$, $T = {T_value}$)',
                  fontsize=13, pad=10)
    ax1.legend(fontsize=9, loc='upper right')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, n_steps)

    # --- Magnetizzazione ---
    ax2.axvspan(0, thermalization_steps, alpha=0.3, color='red',
                label=f'Termalizzazione ({thermalization_steps} steps)', zorder=1)
    ax2.axvspan(thermalization_steps, n_steps, alpha=0.3, color='green',
                label='Misure (equilibrio)', zorder=1)
    ax2.plot(steps, m, 'r-', linewidth=0.5, zorder=2)
    ax2.axhspan(m_mean - m_std, m_mean + m_std,
                xmin=(thermalization_steps/n_steps), xmax=1.0,
                color='gray', alpha=0.3, label=r'$\pm 1\sigma$', zorder=3)
    ax2.axhline(m_mean, color='black', linestyle='--', linewidth=1.5,
                label=f'Media equilibrio: {m_mean:.3f} $\\pm$ {m_std:.3f}', zorder=4)

    ax2.set_xlabel('Cluster updates', fontsize=12)
    ax2.set_ylabel('Magnetizzazione $|m|$', fontsize=12)
    ax2.set_title('Termalizzazione: Evoluzione della Magnetizzazione', fontsize=13, pad=10)
    ax2.legend(fontsize=9, loc='upper right')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, n_steps)

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'thermalization.png', dpi=300, bbox_inches='tight')
    print("  Grafico thermalization.png salvato")
    plt.close()

    # ====== GRAFICO ZOOM (0-15000 steps) ======
    n_zoom = 15000
    steps_zoom = steps[:n_zoom]
    E_zoom = E[:n_zoom]
    m_zoom = m[:n_zoom]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # --- Energia zoom ---
    ax1.axvspan(0, thermalization_steps, alpha=0.3, color='red',
                label=f'Termalizzazione ({thermalization_steps} steps)', zorder=1)
    ax1.axvspan(thermalization_steps, n_zoom, alpha=0.3, color='green',
                label='Misure (equilibrio)', zorder=1)
    ax1.plot(steps_zoom, E_zoom, 'b-', linewidth=0.5, zorder=2)
    ax1.axhspan(E_mean - E_std, E_mean + E_std,
                xmin=(thermalization_steps/n_zoom), xmax=1.0,
                color='gray', alpha=0.3, label=r'$\pm 1\sigma$', zorder=3)
    ax1.axhline(E_mean, color='black', linestyle='--', linewidth=1.5,
                label=f'Media equilibrio: {E_mean:.3f} $\\pm$ {E_std:.3f}', zorder=4)
    ax1.axvline(1000, color='magenta', linestyle=':', linewidth=1.5, label='1000 steps', zorder=5)
    ax1.axvline(thermalization_steps, color='orange', linestyle=':', linewidth=1.5, label='10000 steps', zorder=5)

    ax1.set_xlabel('Cluster updates', fontsize=12)
    ax1.set_ylabel('Energia per spin $E/N$', fontsize=12)
    ax1.set_title(f'Termalizzazione (Zoom: 0-{n_zoom} steps) - Energia ($L = {L_value}$, $T = {T_value}$)',
                  fontsize=13, pad=10)
    ax1.legend(fontsize=9, loc='upper right')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, n_zoom)

    # --- Magnetizzazione zoom ---
    ax2.axvspan(0, thermalization_steps, alpha=0.3, color='red',
                label=f'Termalizzazione ({thermalization_steps} steps)', zorder=1)
    ax2.axvspan(thermalization_steps, n_zoom, alpha=0.3, color='green',
                label='Misure (equilibrio)', zorder=1)
    ax2.plot(steps_zoom, m_zoom, 'r-', linewidth=0.5, zorder=2)
    ax2.axhspan(m_mean - m_std, m_mean + m_std,
                xmin=(thermalization_steps/n_zoom), xmax=1.0,
                color='gray', alpha=0.3, label=r'$\pm 1\sigma$', zorder=3)
    ax2.axhline(m_mean, color='black', linestyle='--', linewidth=1.5,
                label=f'Media equilibrio: {m_mean:.3f} $\\pm$ {m_std:.3f}', zorder=4)
    ax2.axvline(thermalization_steps, color='orange', linestyle=':', linewidth=1.5, label='10000 steps', zorder=5)

    ax2.set_xlabel('Cluster updates', fontsize=12)
    ax2.set_ylabel('Magnetizzazione $|m|$', fontsize=12)
    ax2.set_title(f'Termalizzazione (Zoom: 0-{n_zoom} steps) - Magnetizzazione', fontsize=13, pad=10)
    ax2.legend(fontsize=9, loc='upper right')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, n_zoom)

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'thermalization_zoom.png', dpi=300, bbox_inches='tight')
    print("  Grafico thermalization_zoom.png salvato")
    plt.close()

    # Stampa statistiche
    data_type_str = "dati reali" if is_real_data else "dati simulati"
    print(f"\n  Statistiche termalizzazione ({data_type_str}):")
    print(f"    Energia a equilibrio ({thermalization_steps}-{n_steps}):    {E_mean:.4f} +/- {E_std:.4f}")
    print(f"    Magnetizzazione a equilibrio:          {m_mean:.4f} +/- {m_std:.4f}")

    if not is_real_data:
        print(f"\n  NOTA: Per usare dati reali, esegui 'make thermalization' per generare")
        print(f"        thermalization_data.dat, poi ri-esegui 'make analyze'.")


def plot_heat_with_errors(data):
    """
    Grafico del calore specifico con barre d'errore visibili.
    """
    fig, ax = plt.subplots(figsize=(8, 5))

    L_values = sorted(data.keys())

    for i, L in enumerate(L_values):
        color = get_color(i)
        ax.errorbar(data[L]['T'], data[L]['C'], yerr=data[L]['err_C'],
                    fmt='-', color=color, linewidth=1, capsize=2,
                    markersize=3, label=f'L={L}')

    ax.axvline(T_C_EXACT, color='red', linestyle='--', linewidth=1.5,
               label=f'$T_c$ teorico = {T_C_EXACT}')

    ax.set_xlabel('Temperatura $T$ [$J/k_B$]', fontsize=12)
    ax.set_ylabel('Calore Specifico $C$', fontsize=12)
    ax.set_title("Calore Specifico vs Temperatura (con barre d'errore)", fontsize=13, pad=10)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'heat_with_errors.png', dpi=300, bbox_inches='tight')
    print("  Grafico heat_with_errors.png salvato")
    plt.close()


# Esegui lo script
if __name__ == '__main__':
    main()

    # Genera anche il grafico di termalizzazione
    print("\n" + "="*70)
    print("  GRAFICO TERMALIZZAZIONE")
    print("="*70)
    plot_thermalization()