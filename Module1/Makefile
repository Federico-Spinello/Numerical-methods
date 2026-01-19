# Compilatore e flags
CC = gcc
CFLAGS = -O3 -Wall -Wextra -march=native -std=c99 -I$(INCLUDE_DIR)
LDFLAGS = -lm
DEBUGFLAGS = -g -O0 -DDEBUG -I$(INCLUDE_DIR)

# Directory
SRC_DIR = src
INCLUDE_DIR = include
BIN_DIR = bin
DATA_DIR = data
PLOTS_DIR = plots
SCRIPTS_DIR = scripts
PAPER_DIR = paper

# Files
SOURCES = $(SRC_DIR)/main.c $(SRC_DIR)/ising_wolff.c $(SRC_DIR)/geometry.c \
          $(SRC_DIR)/random.c $(SRC_DIR)/pcg32min.c
HEADERS = $(INCLUDE_DIR)/ising_wolff.h $(INCLUDE_DIR)/geometry.h \
          $(INCLUDE_DIR)/random.h $(INCLUDE_DIR)/pcg32min.h
OBJECTS = $(BIN_DIR)/main.o $(BIN_DIR)/ising_wolff.o $(BIN_DIR)/geometry.o \
          $(BIN_DIR)/random.o $(BIN_DIR)/pcg32min.o
EXECUTABLE = $(BIN_DIR)/ising_simulation

# Python e Virtual Environment
VENV_DIR = venv
PYTHON = python3
VENV_PYTHON = $(VENV_DIR)/bin/python3
VENV_PIP = $(VENV_DIR)/bin/pip3

# ============================================================================
# Target principali
# ============================================================================

.PHONY: all clean run analyze full help dirs venv-setup cleanvenv paper cleanpaper update-paper fullpaper thermalization

# Default: compila
all: $(EXECUTABLE)

# Compila l'eseguibile
$(EXECUTABLE): $(OBJECTS) | $(BIN_DIR)
	@echo "Linking $@..."
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)
	@echo "✓ Build completato: $@"

# Crea bin directory se non esiste
$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

# Compila oggetti (da src/ a bin/)
$(BIN_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS) | $(BIN_DIR)
	@echo "Compiling $<..."
	$(CC) $(CFLAGS) -c $< -o $@

# Crea directory necessarie
dirs:
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(DATA_DIR)
	@mkdir -p $(PLOTS_DIR)
	@echo "✓ Directory create"

# ============================================================================
# Esecuzione simulazioni
# ============================================================================

# Esegui simulazioni in PARALLELO (usa tutti i core)
run: $(EXECUTABLE) dirs
	@echo "=== AVVIO SIMULAZIONI PARALLELE ==="
	@echo "Questo userà tutti i core disponibili per velocizzare l'esecuzione"
	@if [ -d $(VENV_DIR) ]; then \
		$(VENV_PYTHON) $(SCRIPTS_DIR)/parallel_run.py; \
	else \
		$(PYTHON) $(SCRIPTS_DIR)/parallel_run.py; \
	fi
	@echo "✓ Simulazioni parallele completate"

# ============================================================================
# Python Virtual Environment
# ============================================================================

# Crea e configura il virtual environment
$(VENV_DIR)/bin/activate:
	@touch $(VENV_DIR)/bin/activate

# Setup completo del venv (target utente)
venv-setup: $(VENV_DIR)/bin/activate
	@echo "✓ Virtual environment configurato"
	@echo "Pacchetti installati:"
	@$(VENV_PIP) list | grep -E '(numpy|matplotlib|scipy)'

# ============================================================================
# Analisi dati
# ============================================================================

# Analizza tutti i dati e genera grafici (gestisce automaticamente il venv)
analyze: $(VENV_DIR)/bin/activate
	@echo "=== ANALISI DATI ==="
	@if [ ! -d "$(DATA_DIR)" ]; then \
		echo "ERRORE: Directory $(DATA_DIR) non trovata. Esegui 'make run' prima."; \
		exit 1; \
	fi
	@echo "Attivazione virtual environment..."
	@$(VENV_PYTHON) $(SCRIPTS_DIR)/analyze.py
	@echo "✓ Analisi completata"
	@echo "Grafici salvati in $(PLOTS_DIR)/"

# ============================================================================
# Workflow completo
# ============================================================================

# Pipeline completa: compila -> simula -> analizza
full: all run analyze
	@echo ""
	@echo "======================================"
	@echo "✓ PIPELINE COMPLETA TERMINATA"
	@echo "======================================"
	@echo "Eseguibile: $(EXECUTABLE)"
	@echo "Dati:       $(DATA_DIR)/"
	@echo "Grafici:    $(PLOTS_DIR)/"
	@echo ""

# ============================================================================
# Debug e sviluppo
# ============================================================================

# Compila con debug symbols
debug: CFLAGS = $(DEBUGFLAGS) -Wall -Wextra -std=c99
debug: clean $(EXECUTABLE)
	@echo "✓ Build debug completato"

# ============================================================================
# Test termalizzazione
# ============================================================================

# Compila ed esegue test termalizzazione (genera dati reali E(t), m(t) per cold e hot start)
thermalization: $(BIN_DIR)
	@echo "=== COMPILAZIONE TEST TERMALIZZAZIONE ==="
	$(CC) $(CFLAGS) -o $(BIN_DIR)/thermalization_test \
		$(SRC_DIR)/thermalization_test.c \
		$(SRC_DIR)/geometry.c \
		$(SRC_DIR)/random.c \
		$(SRC_DIR)/pcg32min.c \
		$(LDFLAGS)
	@echo "✓ Compilazione completata"
	@echo ""
	@echo "=== ESECUZIONE TEST ==="
	@mkdir -p $(DATA_DIR)
	$(BIN_DIR)/thermalization_test
	@echo ""
	@echo "✓ Test completato. Dati salvati in:"
	@echo "  - $(DATA_DIR)/thermalization_cold.dat (cold start)"
	@echo "  - $(DATA_DIR)/thermalization_hot.dat (hot start)"

# ============================================================================
# LaTeX Paper
# ============================================================================

# Compila paper LaTeX 
paper:
	@echo ""
	@echo "=== COMPILAZIONE PAPER ==="
	@if ! command -v pdflatex >/dev/null 2>&1; then \
		echo "ERRORE: pdflatex non installato!"; \
		echo "Installa con: sudo apt install texlive-lang-italian texlive-science texlive-publishers"; \
		exit 1; \
	fi
	@mkdir -p $(PAPER_DIR)
	@echo "Rimozione PDF precedente..."
	@rm -f $(PAPER_DIR)/paper.pdf
	@echo "Compilazione paper.tex (prima passata)..."
	cd $(PAPER_DIR) && pdflatex -interaction=nonstopmode paper.tex > /dev/null
	@echo "Seconda passata (per riferimenti)..."
	cd $(PAPER_DIR) && pdflatex -interaction=nonstopmode paper.tex > /dev/null
	@echo "✓ Paper compilato: $(PAPER_DIR)/paper.pdf"

# Pipeline completa: simula -> analizza -> aggiorna paper -> compila PDF
fullpaper: all run analyze paper
	@echo ""
	@echo "======================================"
	@echo "✓ PAPER COMPLETO GENERATO"
	@echo "======================================"
	@echo "File finale: $(PAPER_DIR)/paper.pdf"
	@echo ""

# Rimuovi file ausiliari LaTeX
cleanpaper:
	@echo "Pulizia file LaTeX ausiliari..."
	rm -f $(PAPER_DIR)/paper.aux $(PAPER_DIR)/paper.log $(PAPER_DIR)/paper.out \
	      $(PAPER_DIR)/paper.toc $(PAPER_DIR)/paper.bbl $(PAPER_DIR)/paper.blg $(PAPER_DIR)/paper.pdf
	@echo "✓ File LaTeX puliti"

# ============================================================================
# Pulizia
# ============================================================================

# Rimuovi file oggetto e eseguibile
clean:
	@echo "Pulizia file compilati..."
	rm -f $(OBJECTS)
	rm -rf $(BIN_DIR)
	@echo "✓ Pulizia completata"

# Rimuovi anche dati e grafici
cleanall: clean
	@echo "Rimozione dati e grafici..."
	rm -rf $(DATA_DIR)
	rm -rf $(PLOTS_DIR)
	@echo "✓ Pulizia totale completata"


# ============================================================================
# Help
# ============================================================================

help:
	@echo "======================================"
	@echo "  Makefile - Simulazione Ising 2D"
	@echo "======================================"
	@echo ""
	@echo "Target disponibili:"
	@echo ""
	@echo "  make all          - Compila il programma"
	@echo "  make run - Esegue simulazioni in PARALLELO (usa tutti i core)"
	@echo "  make test         - Esegue test rapido"
	@echo "  make analyze      - Analizza dati e genera grafici + JSON (auto venv)"
	@echo "  make full         - Pipeline completa: compila->simula->analizza"
	@echo ""
	@echo "  make update	   - Aggiorna paper.tex con risultati da JSON"
	@echo "  make paper        - Aggiorna paper.tex e compila PDF"
	@echo "  make fullpaper    - Pipeline completa: simula->analizza->paper PDF"
	@echo "  make cleanpaper   - Rimuove file ausiliari LaTeX"
	@echo ""
	@echo "  make venv-setup   - Configura virtual environment Python"
	@echo "  make clean        - Rimuove file compilati"
	@echo "  make cleanall     - Rimuove tutto (dati + grafici + compilati)"
	@echo "  make help         - Mostra questo messaggio"