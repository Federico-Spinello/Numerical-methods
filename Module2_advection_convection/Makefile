# Compilatore e flags
CC = gcc
CFLAGS = -O2 -lm -g -Iinclude -lfftw3

# Directory e file
TARGET = sim
SRC = ./src/main.c src/sim_adv_conv.c src/sim_shock_analysis.c src/params.c src/functions.c
DATA = ./data
SCRIPTS = ./scripts/

# Python e Virtual Environment
VENV_DIR = venv
PYTHON = python3
VENV_PYTHON = $(VENV_DIR)/bin/python3

.PHONY: all setup run shock plot pshock bestfit paper paper-clean clean cleanscreen help

# Default: compila
all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(SRC) -o $(TARGET) $(CFLAGS)

# Setup: crea venv e installa dipendenze
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

# Simulazione standard
run: $(TARGET)
	rm -f ./$(DATA)/data_*.dat ./$(DATA)/fft_*.dat
	./$(TARGET) sim
	make bestfit
	make plot

# Analisi shock
shock: $(TARGET)
	rm -f ./$(DATA)/shock_values.txt
	./$(TARGET) shock
	make pshock

# Plot
plot:
	$(VENV_PYTHON) $(SCRIPTS)plot.py

pshock:
	$(VENV_PYTHON) $(SCRIPTS)plot_shock.py

bestfit:
	$(VENV_PYTHON) $(SCRIPTS)find_best_fit.py

# Paper LaTeX
paper:
	cd paper && pdflatex -interaction=nonstopmode -file-line-error paper.tex || true
	cd paper && pdflatex -interaction=nonstopmode -file-line-error paper.tex || true

paper-clean:
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.toc

# Pulizia
clean:
	rm -f $(TARGET) ./$(DATA)/data_*.dat ./$(DATA)/fft_*.dat ./$(DATA)/shock_values.txt

cleanscreen:
	rm -f ./screen/*

# Help
help:
	@echo "======================================"
	@echo "  Makefile - Equazione di Burgers 1D"
	@echo "======================================"
	@echo ""
	@echo "Setup iniziale:"
	@echo "  make setup        Crea venv e installa dipendenze"
	@echo ""
	@echo "Compilazione ed esecuzione:"
	@echo "  make all          Compila il programma"
	@echo "  make run          Simulazione standard + plot"
	@echo "  make shock        Analisi shock vs viscosita"
	@echo ""
	@echo "Plot:"
	@echo "  make plot         Plot simulazione standard"
	@echo "  make pshock       Plot analisi shock"
	@echo "  make bestfit      Trova best fit spettrale"
	@echo ""
	@echo "Paper:"
	@echo "  make paper        Compila paper LaTeX"
	@echo "  make paper-clean  Rimuove file ausiliari LaTeX"
	@echo ""
	@echo "Pulizia:"
	@echo "  make clean        Rimuove eseguibile e dati"
	@echo "  make cleanscreen  Rimuove screenshot"
	@echo ""
	@echo "  make help         Mostra questo messaggio"
