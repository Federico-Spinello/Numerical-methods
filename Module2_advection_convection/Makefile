CC = gcc
CFLAGS = -O2 -lm -g -Iinclude -lfftw3
TARGET = sim
SRC = ./src/main.c src/sim_adv_conv.c src/sim_shock_analysis.c src/params.c src/functions.c
PYTHON = ./venv/bin/python3
DATA = ./data
SCRIPTS = ./scripts/

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(SRC) -o $(TARGET) $(CFLAGS)

shock: $(TARGET)
	rm -f ./$(DATA)/shock_values.txt
	./$(TARGET) shock
	make pshock

plot:
	$(PYTHON) $(SCRIPTS)plot.py

bestfit:
	$(PYTHON) $(SCRIPTS)find_best_fit.py

run: $(TARGET)
	rm -f ./$(DATA)/data_*.dat ./$(DATA)/fft_*.dat
	./$(TARGET) sim
	make bestfit
	make plot


pshock:
	$(PYTHON) $(SCRIPTS)plot_shock.py

paper:
	cd paper && pdflatex -interaction=nonstopmode -file-line-error paper.tex || true
	cd paper && pdflatex -interaction=nonstopmode -file-line-error paper.tex || true

paper-clean:
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.toc

clean:
	rm -f $(TARGET) ./$(DATA)/data_*.dat ./$(DATA)/fft_*.dat ./$(DATA)/shock_values.txt
cleanscreen:
	rm -f $(TARGET) ./screen/*

.PHONY: all run shock plot plot_shock paper paper-clean clean
