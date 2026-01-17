import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import glob
from datetime import datetime
import os
from scipy.optimize import curve_fit

# Configurazione globale font size
plt.rcParams.update({
    'font.size': 14,           # Dimensione base
    'axes.labelsize': 20,      # Label assi
    'axes.titlesize': 22,      # Titolo plot
    'xtick.labelsize': 14,     # Tick x
    'ytick.labelsize': 14,     # Tick y
    'legend.fontsize': 18,     # Legenda
})

# Funzione di fit: A * k^(-2) + B (offset)
def power_law(k, A, P):
    return A * k**(-P)

plt.ion()

# Caricamento file
files = sorted(
    glob.glob("./data/data_*.dat"),
    key=lambda x: int(x.split('_')[1].split('.')[0])
)
print(f"Trovati {len(files)} file da plottare")
print("Comandi: SPAZIO = pausa/resume, R = restart, Q = quit, S = salva screenshot (in pausa)")

# Lista per salvare i frame della GIF
gif_frames = []

# Crea cartella screen se non esiste
os.makedirs("./screen", exist_ok=True)

# Variabili di controllo
paused = False
restart = False
quit_flag = False

def on_key(event):
    global paused, restart, quit_flag
    if event.key == ' ':
        paused = not paused
        status = "PAUSA" if paused else "RIPRODUZIONE"
        print(f"→ {status}")
    elif event.key == 'r':
        restart = True
        paused = False
        print("→ RESTART")
    elif event.key == 'q':
        quit_flag = True
        paused = False
        print("→ QUIT")
    elif event.key == 'p':
        if paused:
            # Genera nome file con data e ora
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"./screen/screen_{timestamp}.pdf"
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"→ Screenshot salvato: {filename}")
        else:
            print("→ Metti in pausa (SPAZIO) prima di salvare!")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
fig.canvas.mpl_connect('key_press_event', on_key)

while True:
    restart = False
    
    for idx, filename in enumerate(files):
        # Controllo quit
        if quit_flag:
            break
            
        # Controllo restart
        if restart:
            print("Riavvio animazione...")
            break
            
        # Controllo pausa
        while paused:
            plt.pause(0.1)
            if restart or quit_flag:
                break
        
        if restart or quit_flag:
            break
            
        step = int(filename.split('_')[1].split('.')[0])
        
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
        
        # FFT reale
        U = np.fft.rfft(u_fft)
        psd = (np.abs(U)**2) / n_fft**2
        k = 2 * np.pi * np.arange(len(U)) / L
        
        # Pulizia axes
        ax1.clear()
        ax2.clear()

        
        # --- Plot funzione reale ---
        ax1.plot(x, u, color='b', linewidth=2.5)
        ax1.set_xlabel("x")
        ax1.set_ylabel("u(x)")
        ax1.set_title(f"Funzione u(x) - Step {step}")
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(-1.5, 1.5)
        ax1.set_xlim(x[0], x[-1])

        # --- Plot spettro ---
        ax2.plot(k, psd, color='r', linewidth=2.5)
        ax2.set_xlabel("k")
        ax2.set_ylabel("|û(k)|²/N²")
        ax2.set_title(f"Spettro FFT - Step {step}")
        ax2.grid(True, alpha=0.3)
        ax2.set_xlim(1, k.max())
        ax2.set_ylim(10**(-10), 10**(0))
        ax2.set_xscale('log')
        ax2.set_yscale('log')

        plt.tight_layout()
        plt.draw()
        plt.pause(0.01)
    
    # Se quit o se non c'è stato restart, esci dal loop
    if quit_flag or not restart:
        break

if quit_flag:
    print("Uscita dall'animazione")
else:
    print("Animazione completata")

plt.ioff()
plt.show()