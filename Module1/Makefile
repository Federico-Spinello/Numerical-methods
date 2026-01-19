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
VENV_PIP = $(VENV_DIR)/bin/pip

# ============================================================================
# Target principali
# ============================================================================

.PHONY: all clean run analyze full help setup paper cleanpaper fullpaper thermalization cleanall debug

# Default: compila
all: $(EXECUTABLE)

# Compila l'eseguibile
$(EXECUTABLE): $(OBJECTS) | $(BIN_DIR)
	@echo "Linking $@..."
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)
	@echo "Build completato: $@"

# Crea bin directory se non esiste
$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

# Compila oggetti (da src/ a bin/)
$(BIN_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS) | $(BIN_DIR)
	@echo "Compiling $<..."
	$(CC) $(CFLAGS) -c $< -o $@

# ============================================================================
# Setup ambiente
# ============================================================================

# Setup completo: crea venv e installa dipendenze
setup:
	@echo "=== SETUP AMBIENTE ==="
	@echo "Creazione virtual environment..."
	@rm -rf $(VENV_DIR)
	@$(PYTHON) -m venv $(VENV_DIR)
	@echo "Installazione dipendenze da requirements.txt..."
	@./$(VENV_DIR)/bin/pip install --upgrade pip > /dev/null
	@./$(VENV_DIR)/bin/pip install -r requirements.txt
	@echo ""
	@echo "Setup completato!"
	@echo "Pacchetti installati:"
	@./$(VENV_DIR)/bin/pip list | grep -E '(numpy|matplotlib|scipy)'

# ============================================================================
# Esecuzione simulazioni
# ============================================================================

# Esegui simulazioni in PARALLELO (usa tutti i core)
run: $(EXECUTABLE)
	@echo "=== AVVIO SIMULAZIONI PARALLELE ==="
	@mkdir -p $(DATA_DIR)
	@mkdir -p $(PLOTS_DIR)
	@$(VENV_PYTHON) $(SCRIPTS_DIR)/parallel_run.py
	@echo "Simulazioni completate"

# ============================================================================
# Analisi dati
# ============================================================================

# Analizza tutti i dati e genera grafici
analyze:
	@echo "=== ANALISI DATI ==="
	@if [ ! -d "$(DATA_DIR)" ]; then \
		echo "ERRORE: Directory $(DATA_DIR) non trovata. Esegui 'make run' prima."; \
		exit 1; \
	fi
	@$(VENV_PYTHON) $(SCRIPTS_DIR)/analyze.py
	@echo "Analisi completata"
	@echo "Grafici salvati in $(PLOTS_DIR)/"

# ============================================================================
# Workflow completo
# ============================================================================

# Pipeline completa: compila -> simula -> analizza
full: all run analyze
	@echo ""
	@echo "======================================"
	@echo "PIPELINE COMPLETA TERMINATA"
	@echo "======================================"
	@echo "Eseguibile: $(EXECUTABLE)"
	@echo "Dati:       $(DATA_DIR)/"
	@echo "Grafici:    $(PLOTS_DIR)/"

# ============================================================================
# Test termalizzazione
# ============================================================================

# Compila ed esegue test termalizzazione (genera dati E(t), m(t) per cold e hot start)
thermalization: $(BIN_DIR)
	@echo "=== TEST TERMALIZZAZIONE ==="
	$(CC) $(CFLAGS) -o $(BIN_DIR)/thermalization_test \
		$(SRC_DIR)/thermalization_test.c \
		$(SRC_DIR)/geometry.c \
		$(SRC_DIR)/random.c \
		$(SRC_DIR)/pcg32min.c \
		$(LDFLAGS)
	@mkdir -p $(DATA_DIR)
	$(BIN_DIR)/thermalization_test
	@echo ""
	@echo "Dati salvati in:"
	@echo "  - $(DATA_DIR)/thermalization_cold.dat"
	@echo "  - $(DATA_DIR)/thermalization_hot.dat"

# ============================================================================
# LaTeX Paper
# ============================================================================

# Compila paper LaTeX
paper:
	@echo "=== COMPILAZIONE PAPER ==="
	@if ! command -v pdflatex >/dev/null 2>&1; then \
		echo "ERRORE: pdflatex non installato!"; \
		echo "Installa con: sudo apt install texlive-lang-italian texlive-science texlive-publishers"; \
		exit 1; \
	fi
	@mkdir -p $(PAPER_DIR)
	@rm -f $(PAPER_DIR)/paper.pdf
	@cd $(PAPER_DIR) && pdflatex -interaction=nonstopmode paper.tex > /dev/null
	@cd $(PAPER_DIR) && pdflatex -interaction=nonstopmode paper.tex > /dev/null
	@echo "Paper compilato: $(PAPER_DIR)/paper.pdf"

# Pipeline completa: simula -> analizza -> compila PDF
fullpaper: all run analyze paper
	@echo ""
	@echo "======================================"
	@echo "PAPER COMPLETO GENERATO"
	@echo "======================================"
	@echo "File finale: $(PAPER_DIR)/paper.pdf"

# Rimuovi file ausiliari LaTeX
cleanpaper:
	@echo "Pulizia file LaTeX..."
	@rm -f $(PAPER_DIR)/paper.aux $(PAPER_DIR)/paper.log $(PAPER_DIR)/paper.out \
	       $(PAPER_DIR)/paper.toc $(PAPER_DIR)/paper.bbl $(PAPER_DIR)/paper.blg $(PAPER_DIR)/paper.pdf
	@echo "File LaTeX puliti"

# ============================================================================
# Debug
# ============================================================================

# Compila con debug symbols
debug: CFLAGS = $(DEBUGFLAGS) -Wall -Wextra -std=c99
debug: clean $(EXECUTABLE)
	@echo "Build debug completato"

# ============================================================================
# Pulizia
# ============================================================================

# Rimuovi file compilati
clean:
	@echo "Pulizia file compilati..."
	@rm -f $(OBJECTS)
	@rm -rf $(BIN_DIR)
	@echo "Pulizia completata"

# Rimuovi tutto (dati + grafici + compilati)
cleanall: clean cleanpaper
	@echo "Rimozione dati e grafici..."
	@rm -rf $(DATA_DIR)
	@rm -rf $(PLOTS_DIR)
	@echo "Pulizia totale completata"

# ============================================================================
# Help
# ============================================================================

help:
	@echo "======================================"
	@echo "  Makefile - Simulazione Ising 2D"
	@echo "======================================"
	@echo ""
	@echo "Setup iniziale:"
	@echo "  make setup           Crea venv e installa dipendenze"
	@echo ""
	@echo "Compilazione ed esecuzione:"
	@echo "  make all             Compila il programma"
	@echo "  make run             Esegue simulazioni (parallelo)"
	@echo "  make analyze         Analizza dati e genera grafici"
	@echo "  make full            Pipeline: compila + simula + analizza"
	@echo ""
	@echo "Test e paper:"
	@echo "  make thermalization  Test termalizzazione (cold/hot start)"
	@echo "  make paper           Compila paper LaTeX"
	@echo "  make fullpaper       Pipeline completa + paper"
	@echo ""
	@echo "Pulizia:"
	@echo "  make clean           Rimuove file compilati"
	@echo "  make cleanpaper      Rimuove file ausiliari LaTeX"
	@echo "  make cleanall        Rimuove tutto"
	@echo ""
	@echo "  make help            Mostra questo messaggio"
