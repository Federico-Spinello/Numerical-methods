# Metodi Numerici per la Fisica

Progetti sviluppati per l'esame di **Metodi Numerici per la Fisica**.

---

## Moduli

### [Module 1 - Ising 2D con Algoritmo di Wolff](Module1/)

Simulazione Monte Carlo del **modello di Ising 2D** per lo studio della transizione di fase ferromagnetica.

**Obiettivi:**
- Determinazione della temperatura critica Tc tramite Binder cumulant, picchi di suscettivita e calore specifico
- Stima dell'esponente critico gamma/nu tramite Finite Size Scaling
- Validazione dei risultati con la soluzione esatta di Onsager (1944)

**Metodi:**
- Algoritmo di Wolff (cluster algorithm) per ridurre il critical slowing down
- Analisi Finite Size Scaling con data collapse
- Bootstrap non parametrico per stima errori

---

### [Module 2 - Equazione di Burgers 1D](Module2_advection_convection/)

Simulazione numerica dell'**equazione di Burgers** in una dimensione con analisi spettrale FFT.

**Obiettivi:**
- Studio della formazione di shock in fluidi
- Analisi della competizione tra avvezione (steepening) e diffusione (smoothing)
- Relazione tra viscosita e pendenza degli shock

**Metodi:**
- Integrazione temporale Runge-Kutta 4 (RK4)
- Derivate spaziali con schemi centrati al 4 ordine
- Time step adattivo per stabilita (condizioni CFL)
- Analisi spettrale con FFTW3

**Analisi:**
- Simulazione standard con visualizzazione animata
- Analisi parametrica shock vs viscosita

---

## Quick Start

```bash
# Module 1 - Ising 2D
cd Module1
make setup      # Installa dipendenze Python
make full       # Compila, simula, analizza

# Module 2 - Burgers
cd Module2_advection_convection
make            # Compila
make run        # Simula e visualizza
```

---

**Autore:** Federico Spinello
**Corso:** Metodi Numerici per la Fisica
**Anno:** 2025/2026
