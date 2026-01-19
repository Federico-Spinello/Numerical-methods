#!/usr/bin/env python3
"""
Analizza tutti i timestep, trova quello con il miglior fit usando come criterio
la distanza minima da P = 2.0 in unità di sigma (non chi-quadro).
Genera un grafico statico salvato in PDF.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Backend non-interattivo
import matplotlib.pyplot as plt
import glob
from scipy.optimize import curve_fit
import os

# Configurazione globale font size
plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 20,
    'axes.titlesize': 22,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 18,
})

# Funzione di fit: A * k^(-P)
def power_law(k, A, P):
    return A * k**(-P)

# Caricamento file
files = sorted(
    glob.glob("./data/data_*.dat"),
    key=lambda x: int(x.split('_')[1].split('.')[0])
)

print(f"Trovati {len(files)} file da analizzare")

# Variabili per tracciare il miglior fit
# Cerchiamo il fit dove P è più vicino a 2.0 in unità di sigma
best_distance_from_2 = np.inf  # Distanza |P - 2| / sigma_P
best_step = -1
best_data = None
best_fit_params = None
best_fit_errors = None
best_chi2 = None
best_chi2_reduced = None

# Range di fit
k_min = 6
k_max = 900
psd_min = 1e-10

print("\nAnalisi dei fit per ogni timestep...")

for idx, filename in enumerate(files):
    step = int(filename.split('_')[1].split('.')[0])

    try:
        # Caricamento dati
        data = np.loadtxt(filename)
        x = data[:, 0]
        u = data[:, 1]
        dx = x[1] - x[0]
        n = len(x)

        # Rimuoviamo l'ultimo punto se è duplicato (periodicità)
        if np.isclose(u[0], u[-1]):
            u_fft = u[:-1]
            n_fft = len(u_fft)
            L = dx * n_fft
        else:
            u_fft = u
            n_fft = n
            L = dx * n_fft

        # Carica lo spettro FFT dal file generato in C (che include gli errori)
        fft_filename = filename.replace('data_', 'fft_')
        fft_data = np.loadtxt(fft_filename)
        k = fft_data[:, 0]
        psd = fft_data[:, 1]
        sigma = fft_data[:, 2]

        # Range di fit
        k_fit_range = (k > k_min) & (k < k_max) & (psd > psd_min)

        if np.sum(k_fit_range) > 5:  # Almeno 5 punti per fittare
            k_fit = k[k_fit_range]
            psd_fit = psd[k_fit_range]
            sigma_fit = sigma[k_fit_range]

            # Proteggi contro sigma = 0
            sigma_fit[sigma_fit < 1e-15] = 1e-15

            # Fit con curve_fit usando gli errori come pesi
            popt, pcov = curve_fit(power_law, k_fit, psd_fit,
                                p0=[1e-2, 2],
                                sigma=sigma_fit,  # IMPORTANTE: usa gli errori dal C
                                absolute_sigma=True,  # Usa sigma come errori assoluti
                                bounds=([0, 0], [np.inf, np.inf]))

            A_fit, P_fit = popt
            A_err, P_err = np.sqrt(np.diag(pcov))

            psd_model = power_law(k_fit, A_fit, P_fit)
            residuals = psd_fit - psd_model
            P_target = 2.0
            distance_from_2 = abs(P_fit - P_target) / P_err  # Numero di sigma

            # Aggiorna il miglior fit se P è più vicino a 2.0
            if distance_from_2 < best_distance_from_2:
                best_distance_from_2 = distance_from_2
                best_step = step
                best_data = {
                    'x': x,
                    'u': u,
                    'k': k,
                    'psd': psd,
                    'k_fit': k_fit,
                    'psd_fit': psd_fit
                }
                best_fit_params = (A_fit, P_fit)
                best_fit_errors = (A_err, P_err)

            if (idx + 1) % 50 == 0:
                print(f"Processati {idx + 1}/{len(files)} file...")

    except Exception as e:
        # Salta file problematici
        continue

print(f"\n{'='*60}")
print(f"MIGLIOR FIT TROVATO:")
print(f"Criterio: minima distanza da P = 2.0 in unità di sigma")
print(f"{'='*60}")
print(f"Step temporale: {best_step}")
print(f"Parametri fit:")
print(f"  A = {best_fit_params[0]:.4e} ± {best_fit_errors[0]:.4e}")
print(f"  P = {best_fit_params[1]:.6f} ± {best_fit_errors[1]:.6f}")
print(f"\nDistanza da P = 2.0:")
print(f"  |P - 2.0| / σ_P = {best_distance_from_2:.4f} σ")
print(f"  Valore teorico: P = 2.0 (regime viscoso)")

print(f"{'='*60}\n")

# Crea il grafico
if best_data is not None:
    print("Generazione grafico...")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # --- Plot funzione reale ---
    ax1.plot(best_data['x'], best_data['u'], color='b', linewidth=2.5)
    ax1.set_xlabel("x")
    ax1.set_ylabel("u(x)")
    ax1.set_title(f"Funzione u(x) - Step {best_step}")
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(-1.5, 1.5)
    ax1.set_xlim(best_data['x'][0], best_data['x'][-1])

    # --- Plot spettro ---
    ax2.plot(best_data['k'], best_data['psd'], color='r', linewidth=2.5, label='Dati')

    # Plot del fit
    k_fine = np.logspace(np.log10(best_data['k_fit'].min()),
                         np.log10(best_data['k_fit'].max()), 200)
    psd_model = power_law(k_fine, *best_fit_params)

    A_fit, P_fit = best_fit_params
    A_err, P_err = best_fit_errors

    # Calcola distanza da P=2 in sigma
    distance_sigma = best_distance_from_2

    ax2.plot(k_fine, psd_model, 'k--', linewidth=2.5,
            label=f'Fit: $A k^{{-P}}$\n$A = {A_fit:.2e} \\pm {A_err:.1e}$\n$P = {P_fit:.5f} \\pm {P_err:.5f}$\n$|P-2|/\\sigma_P = {distance_sigma:.2f}\\sigma$\n')

    ax2.set_xlabel("k")
    ax2.set_ylabel("|û(k)|²/N²")
    ax2.set_title(f"Spettro FFT - Best Fit (Step {best_step})")
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(1, best_data['k'].max())
    ax2.set_ylim(10**(-10), 10**(0))
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.legend(loc='lower left')

    plt.tight_layout()

    # Salva in PDF
    os.makedirs("./images", exist_ok=True)
    output_path = "./images/fit.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nGrafico salvato in: {output_path}")

    # Salva anche i dati del miglior fit
    fit_info_path = "./best_fit_info.txt"
    with open(fit_info_path, 'w') as f:
        f.write(f"Miglior fit trovato per l'analisi spettrale\n")
        f.write(f"Criterio: minima distanza da P = 2.0 in unità di sigma\n")
        f.write(f"{'='*60}\n")
        f.write(f"Step temporale: {best_step}\n")
        f.write(f"\nParametri fit:\n")
        f.write(f"  A = {best_fit_params[0]:.6e} ± {best_fit_errors[0]:.6e}\n")
        f.write(f"  P = {best_fit_params[1]:.8f} ± {best_fit_errors[1]:.8f}\n")
        f.write(f"\nDistanza da P = 2.0 (valore teorico regime viscoso):\n")
        f.write(f"  |P - 2.0| / sigma_P = {best_distance_from_2:.6f} sigma\n")
        f.write(f"  Deviazione assoluta (|P - 2.0|): {abs(best_fit_params[1] - 2.0):.6e}\n")
        f.write(f"\nRange di fit:\n")
        f.write(f"  k_min = {k_min}\n")
        f.write(f"  k_max = {k_max}\n")
        f.write(f"  psd_min = {psd_min}\n")
        f.write(f"  Numero punti usati: {len(best_data['k_fit'])}\n")
        f.write(f"\nInterpretazione:\n")
        if best_distance_from_2 < 1.0:
            f.write(f"  ECCELLENTE: P compatibile con 2.0 entro 1σ\n")
            f.write(f"  Lo spettro segue la legge k^-2 del regime viscoso\n")
        elif best_distance_from_2 < 2.0:
            f.write(f"  BUONO: P compatibile con 2.0 entro 2σ\n")
            f.write(f"  Lo spettro è consistente con il regime viscoso\n")
        elif best_distance_from_2 < 3.0:
            f.write(f"  ACCETTABILE: P compatibile con 2.0 entro 3σ\n")
            f.write(f"  Possibile regime di transizione\n")
        else:
            f.write(f"  ATTENZIONE: P significativamente diverso da 2.0\n")
            f.write(f"  Regime non-viscoso o transitorio\n")

    print(f"Informazioni salvate in: {fit_info_path}")
else:
    print("Errore: nessun fit valido trovato!")
