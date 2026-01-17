import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from datetime import datetime
import os

# Configurazione globale font size
plt.rcParams.update({
    'font.size': 14,           # Dimensione base
    'axes.labelsize': 20,      # Label assi
    'axes.titlesize': 22,      # Titolo plot
    'xtick.labelsize': 14,     # Tick x
    'ytick.labelsize': 14,     # Tick y
    'legend.fontsize': 18,     # Legenda
})

# File da leggere
shock_file = "./data/shock_values.txt"

# Crea cartella screen se non esiste
os.makedirs("./screen", exist_ok=True)

def on_key(event):
    if event.key == 'p':
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"./screen/shock_plot_{timestamp}.pdf"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"→ Screenshot salvato: {filename}")

# Leggi il file shock_values.txt
if not os.path.exists(shock_file):
    print(f"Errore: file {shock_file} non trovato!")
    print("Esegui prima 'make shock' per generare i dati.")
    exit(1)

data = np.loadtxt(shock_file, comments='#')

if data.size == 0:
    print("File vuoto!")
    exit(1)

# Assicurati che data sia 2D anche con una sola riga
if data.ndim == 1:
    data = data.reshape(1, -1)

nu_all = data[:, 0]
du_all = data[:, 1]

n_points = len(nu_all)
print(f"Trovati {n_points} punti")
print("Premi 'P' per salvare lo screenshot")

# Crea il plot (singolo)
fig, ax = plt.subplots(figsize=(10, 7))
fig.canvas.mpl_connect('key_press_event', on_key)

# --- Plot du vs nu ---
ax.plot(nu_all, -du_all, 'o-', color='b', linewidth=2.5, markersize=10, markerfacecolor='cyan', markeredgewidth=2)
ax.set_xlabel("Viscosità ν")
ax.set_ylabel("-min(du/dx) [Pendenza shock]")
ax.set_title("Pendenza Shock vs Viscosità", fontweight='bold')
ax.grid(True, alpha=0.3)
ax.set_xscale('log')
ax.set_yscale('log')

# Aggiungi fit per verificare scaling con errore
log_nu = np.log10(nu_all)
log_du = np.log10(-du_all)

# Fit lineare
coeffs = np.polyfit(log_nu, log_du, 1)
slope = coeffs[0]
intercept = coeffs[1]

# Calcola errore sulla pendenza
# Residui
log_du_fit = slope * log_nu + intercept
residuals = log_du - log_du_fit
residual_std = np.std(residuals)

# Errore sulla pendenza
n = len(log_nu)
sum_x_squared = np.sum((log_nu - np.mean(log_nu))**2)
slope_error = residual_std / np.sqrt(sum_x_squared)

print(f"\n=== RISULTATI FIT ===")
print(f"Pendenza: {slope:.4f} ± {slope_error:.4f}")
print(f"Intercetta: {intercept:.4f}")
print(f"Deviazione standard residui: {residual_std:.4f}")

# Plot del fit
nu_fit = np.logspace(np.log10(nu_all.min()), np.log10(nu_all.max()), 100)
du_fit = 10**(intercept) * nu_fit**slope

ax.plot(nu_fit, du_fit, 'k--', linewidth=2.5, alpha=0.7,
        label=f'Fit: slope = {slope:.3f} ± {slope_error:.3f}')
ax.legend(loc='best')

plt.tight_layout()
plt.show()