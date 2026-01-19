# Simulatore Equazione di Burgers 1D

Simulazione numerica dell'equazione di Burgers in una dimensione con condizioni al contorno periodiche, integrazione temporale Runge-Kutta 4 e analisi spettrale FFT.

## 📋 Descrizione

Questo progetto implementa un risolutore numerico per l'equazione di Burgers:

```
∂u/∂t + u · ∂u/∂x = ν · ∂²u/∂x²
```

## 🗂️ Struttura del Progetto

```
.
├── Makefile              # Build automation
├── README.md
├── requirements.txt      # Dipendenze Python
├── params.txt            # Parametri simulazione
├── sim                   # Programma compilato (generato da make)
├── src/
│   ├── main.c                  # Entry point e dispatcher
│   ├── sim_adv_conv.c          # Simulazione standard
│   ├── sim_shock_analysis.c    # Analisi shock vs viscosità
│   ├── functions.c             # Derivate, RHS, integratore RK4, FFT
│   └── params.c                # Gestione parametri
├── include/
│   ├── functions.h
│   └── params.h
├── scripts/
│   ├── plot.py                 # Visualizzazione risultati simulazione
│   ├── plot_shock.py           # Visualizzazione analisi shock
│   └── find_best_fit.py        # Analisi best fit spettrale
├── data/                 # Output simulazioni (generato)
│   ├── data_*.dat        # Stati u(x) per ogni timestep salvato
│   ├── fft_*.dat         # Spettri di potenza FFT
│   └── shock_values.txt  # Risultati analisi shock
├── screen/               # Screenshot salvati (generato)
│   └── *.pdf             # Plot salvati con tasto 'P'
├── paper/                # Paper LaTeX
│   └── paper.tex
└── venv/                 # Virtual environment Python (generato da make setup)
```

## 🔧 Dipendenze

### C
- **gcc** - Compilatore C
- **libfftw3** - Fast Fourier Transform library
- **libm** - Libreria matematica standard


## 🚀 Compilazione ed Esecuzione

### Setup iniziale (solo prima volta)
```bash
make setup          # Crea venv e installa dipendenze Python
make                # Compila il programma C
```

### Simulazione Standard
Esegue una simulazione completa con salvataggio dei dati temporali e plot animato:
```bash
make run
```
Questo comando:
1. Pulisce i vecchi dati
2. Esegue `./sim sim` (simulazione standard)
3. Trova il best fit spettrale
4. Genera automaticamente i plot animati

### Analisi Shock
Esegue un'analisi parametrica della pendenza dello shock al variare della viscosità:
```bash
make shock
```
Questo comando:
1. Esegue `./sim shock` (analisi shock)
2. Salva i risultati in `./data/shock_values.txt`
3. Genera automaticamente i plot

### Solo plot (senza ricalcolare)
```bash
make plot           # Plot simulazione standard
make pshock         # Plot analisi shock
make bestfit        # Trova best fit spettrale
```

### Pulizia
```bash
make clean          # Rimuove eseguibile e dati
make cleanscreen    # Rimuove screenshot
make help           # Mostra tutti i comandi
```

## 🧮 Metodi Numerici Implementati

### Schema Temporale: Runge-Kutta 4 (RK4)

Integratore esplicito del quarto ordine con quattro valutazioni per step:

```
k₁ = F(uⁿ)
k₂ = F(uⁿ + dt/2 · k₁)
k₃ = F(uⁿ + dt/2 · k₂)
k₄ = F(uⁿ + dt · k₃)

uⁿ⁺¹ = uⁿ + dt/6 · (k₁ + 2k₂ + 2k₃ + k₄)
```

### Derivate Spaziali

**Derivata Prima** - Schema centrato a 4 punti (4° ordine):
```
du/dx ≈ (u[i-2] - 8·u[i-1] + 8·u[i+1] - u[i+2]) / (12·dx)
```

**Derivata Seconda** - Schema centrato a 5 punti (4° ordine):
```
d²u/dx² ≈ (-u[i+2] + 16·u[i+1] - 30·u[i] + 16·u[i-1] - u[i-2]) / (12·dx²)
```

**Alternative disponibili:**
- `first_derivative_upwind()` - Upwind per c > 0 (più stabile, meno accurato)

### Condizioni al Contorno Periodiche

Gli indici vengono wrappati usando l'operatore modulo:
```c
int ip = (i + 1) % nx;  // punto successivo
int im = (i - 1 + nx) % nx;  // punto precedente
```

### Time Step Adattivo

Il time step viene calcolato ad ogni iterazione per garantire stabilità:

```
dt_diff = CFL_diff · dx² / ν        (condizione diffusiva)
dt_adv = CFL_adv · dx / max|u|      (condizione avvettiva)
dt = safety · min(dt_diff, dt_adv)
```

Parametri di default:
- `CFL_diff = 0.5`
- `CFL_adv = 0.4`
- `safety = 0.9`

## ⚙️ Parametri Principali

### Simulazione Standard (sim_adv_conv.c)

| Parametro | Valore | Descrizione |
|-----------|--------|-------------|
| `nx` | 500 | Numero di celle spaziali |
| `nt` | 3000 | Numero di timesteps |
| `L` | 1.0 | Lunghezza dominio |
| `c` | 10.0 | Coefficiente avvezione * |
| `nu` | dx²·π·140/2 | Coefficiente viscosità |
| `print_step` | 5 | Frequenza salvataggio dati |

### Analisi Shock (sim_shock_analysis.c)

| Parametro | Valore | Descrizione |
|-----------|--------|-------------|
| `nu_min` | dx²·π·100/2 | Viscosità minima |
| `nu_max` | dx²·π·600/2 | Viscosità massima |
| `n_nu` | 100 | Numero di simulazioni |
| `use_log_scale` | 1 | Scala logaritmica (1) o lineare (0) |
| `x_shock` | 0.5 | Posizione dove misurare du/dx |

Il coefficiente di avvezione `c` non è utilizzato perché sto facendo una simulazione non lineare, ma se volessi tornare al caso lineare, mi basterebbe sostituire la `c` al posto della `u[i]` quando calcolo il `RHS`.

## 📈 Output e Visualizzazione

### File Dati Generati

#### Simulazione Standard
- **`data/data_XXXXX.dat`** - Profilo spaziale u(x) al timestep XXXXX
  - Colonna 1: posizione x
  - Colonna 2: valore u(x)

- **`data/fft_XXXXX.dat`** - Spettro di potenza al timestep XXXXX
  - Colonna 1: numero d'onda k
  - Colonna 2: densità spettrale |û(k)|²/N²

#### Analisi Shock
- **`data/shock_values.txt`** - Risultati analisi parametrica
  - Colonna 1: Viscosità ν
  - Colonna 2: Max |du/dx| (pendenza shock)
  - Colonna 3: Time step finale Δt

### Visualizzazione plot.py (Simulazione Standard)

Il codice Python genera un'animazione con due subplot:

1. **Sinistra**: Profilo u(x) in funzione dello spazio
2. **Destra**: Spettro di potenza in scala logaritmica

**Comandi interattivi:**
- `SPAZIO` - Pausa/Resume animazione
- `R` - Restart dall'inizio
- `Q` - Quit
- `P` - Salva screenshot in `./screen/` (solo in pausa)

### Visualizzazione plot_shock.py (Analisi Shock)

Genera due plot statici:

1. **Sinistra**: Pendenza shock (max |du/dx|) vs viscosità ν
2. **Destra**: Time step Δt vs viscosità ν

**Comandi:**
- `P` - Salva screenshot in `./screen/shock_plot_TIMESTAMP.pdf`

## 🔬 Analisi FFT

### Implementazione

Utilizza **FFTW3** (Fastest Fourier Transform in the West) per calcolare trasformate reali:
- `fftw_plan_dft_r2c_1d()` - Piano per trasformata reale → complessa
- Ottimizzato per dati reali (dimezza calcoli e memoria)

### Power Spectral Density (PSD)

```
PSD(k) = |û(k)|² / N²
```

dove:
- `û(k)` è la trasformata di Fourier di u(x)
- N è il numero di punti spaziali
- Normalizzazione consistente con Parseval

### Interpretazione Fisica

- **Picco a k basso**: Strutture grandi (lunghezza d'onda lunga)
- **Energia a k alto**: Strutture piccole, dettagli fini
- **Cascata energetica**: In turbolenza, energia trasferita da k bassi a k alti
- **Dissipazione viscosa**: Taglia le alte frequenze (filtro passa-basso)

## 🔬 Analisi Shock vs Viscosità

### Obiettivo

Studiare come la viscosità influenza la pendenza degli shock che si formano nell'equazione di Burgers.

### Metodologia

1. Per ogni valore di ν nel range specificato:
   - Esegue una simulazione completa fino a nt timesteps
   - Calcola la derivata prima du/dx usando lo schema a 4 punti
   - Misura max|du/dx| in x = 0.5 (posizione dello shock)
   - Registra ν, max|du/dx|, e Δt finale

2. I risultati vengono salvati in ordine crescente di ν in `shock_values.txt`

3. Se si esegue di nuovo l'analisi con nuovi valori di ν, questi vengono inseriti automaticamente nella posizione corretta

### Interpretazione Risultati

- **Viscosità bassa**: Shock più ripidi (|du/dx| grande), transizioni brusche
- **Viscosità alta**: Shock più dolci (|du/dx| piccolo), maggiore diffusione
- **Time step**: Diminuisce con viscosità più bassa (condizione CFL più restrittiva)

### Configurazione Range

Per modificare il range di viscosità, edita `sim_shock_analysis.c`:

```c
double nu_min = dx * dx * pi * 10.0 / 2.0;   // Valore minimo
double nu_max = dx * dx * pi * 300.0 / 2.0;  // Valore massimo
int n_nu = 20;                                // Numero di punti
int use_log_scale = 1;                        // 1=log, 0=lineare
```

## 🧪 Casi Test Suggeriti

### 1. Shock Formation (discontinuità)
```c
init_sin(u, x, nx, L);  
nu = 0.000251;  // (= dx * dx * pi * 4000.0 / 2.0) viscosità molto bassa
```
Osserva la formazione di shock e la cascata energetica verso k alti, successivamente la viscosità ammazza tutto, quindi lo spettro ritorna ad essere piccato verso k molto bassi.

### 2. Diffusione Dominante
```c
double x0 = 0.5; double width = 0.03;
init_gauss(u, x, nx, x0, width);
nu = 0.1;  // viscosità alta
```
Osserva la dissipazione energetica e l'allargamento della gaussiana.

### 3. Turbolenza 1D
```c
init_brutto(u, x, nx, x0, width);  // multi-scala
nu = dx*dx*pi*400.0/2.0;  // intermedio
```
Osserva l'interazione tra scale diverse e redistribuzione energia.

### 4. Analisi Parametrica Shock
```bash
make shock
python plot_shock.py
```
Studia la relazione tra viscosità e ripidità degli shock. Aspettati una relazione di power-law: |du/dx| ∝ ν^(-α).

## 🐛 Debug e Troubleshooting

### Simulazione instabile (NaN, overflow)
- Riduci `CFL_adv` e `CFL_diff`
- Aumenta `safety` factor

### Malloc errors
- Riduci `nx` o `nt` se la memoria è insufficiente
- Controlla che tutti i `malloc` abbiano corrispondente `free`

### Analisi shock troppo lenta
- Riduci `n_nu` (numero di simulazioni)
- Riduci `nt` se non serve raggiungere steady state
- Usa scala logaritmica (`use_log_scale = 1`) per campionare meglio

### Plot non si aggiorna
- Verifica che `matplotlib` usi backend `TkAgg`
- Controlla che i file in `./data/` esistano e non siano vuoti

## 📚 Riferimenti Teorici

**Equazione di Burgers**: Modello semplificato delle equazioni di Navier-Stokes, usato per studiare:
- Formazione di shock in fluidi
- Turbolenza 1D
- Metodi numerici per PDEs non lineari

**Importanza fisica**:
- Competizione tra avvezione (steepening) e diffusione (smoothing)
- Conservazione della massa
- Cascata energetica (Kolmogorov in 3D)
- Relazione shock thickness ~ ν/u (teoria di Rankine-Hugoniot)

## 📝 Note Tecniche

- Il codice usa allocazione dinamica per supportare griglie grandi
- Tutti i calcoli in double precision (64-bit)
- Flags di ottimizzazione: `-O2`
- Flag di debug disponibili: `-g`
- Architettura modulare: `main.c` dispatcher + simulazioni separate
- Salvataggio ordinato automatico per analisi shock

## 🎯 Workflow Tipico

### Esplorazione Iniziale
```bash
make run          # Simulazione singola con visualizzazione
# Premi SPAZIO per pausare, P per salvare screenshot
```

### Analisi Parametrica
```bash
make shock        # Esegue sweep di viscosità
python plot_shock.py  # Visualizza risultati
# Premi P per salvare il grafico
```

### Modifica e Test
```bash
# Edita src/sim_adv_conv.c o src/sim_shock_analysis.c
make              # Ricompila
make run          # Testa
```

## 👤 Autore

Progetto di simulazione numerica per lo studio dell'equazione di Burgers con analisi spettrale e parametrica degli shock di Federico Spinello.