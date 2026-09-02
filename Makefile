# Reproducibility Makefile.  Requires: `pip install -r requirements.txt`
# (and a TeX distribution with pdflatex for `make paper`).
# Usage:  make verify | make figures | make reproduce | make paper | make all | make clean

PYTHON ?= python
LATEX  ?= pdflatex

.PHONY: all verify figures reproduce paper clean help

help:
	@echo "make verify    - fast smoke test: the sum-rule identities + the interval battery catching a truncation the lone first moment misses (seconds)"
	@echo "make figures   - regenerate every native pgfplots data fragment (.dat) from the exact engine"
	@echo "make reproduce - run the master notebook end to end (every figure and number)"
	@echo "make paper     - compile paper/main.tex -> paper/main.pdf"
	@echo "make all       - figures + paper"
	@echo "make clean     - remove LaTeX aux files (keeps main.pdf)"

# ---- fast reproducibility check (numpy/scipy only, seconds) ----
# Recomputes the current-probe moments m0,m1,m2, the shot-budget interval, and shows the joint
# (m0,m1,m2)+Hankel battery catching (via m2) a truncation the lone first moment misses.
verify:
	$(PYTHON) src/verify.py

# ---- regenerate the data-driven figure fragments from the exact sector-Lanczos + Haydock engine ----
figures:
	cd src && $(PYTHON) export_akw_dat.py && $(PYTHON) export_sqw_dat.py && \
	          $(PYTHON) export_falsifier_dat.py && $(PYTHON) export_teeth_dat.py && \
	          $(PYTHON) export_christoffel_dat.py && $(PYTHON) export_bracketing_dat.py && \
	          $(PYTHON) export_momentcone_dat.py
	@echo "OK: native pgfplots .dat fragments regenerated (committed paper/figs/ remains authoritative)."

# ---- reproduce EVERYTHING in one coherent, narrated pass ----
reproduce:
	jupyter nbconvert --to notebook --execute --inplace notebooks/00_Reproduce_Everything.ipynb
	@echo "OK: full pipeline executed in notebooks/00_Reproduce_Everything.ipynb"

# ---- compile the preprint (uses the committed main.bbl; no bibtex needed) ----
paper:
	cd paper && $(LATEX) -interaction=nonstopmode main.tex >/dev/null && \
	            $(LATEX) -interaction=nonstopmode main.tex >/dev/null
	@echo "OK: paper/main.pdf built."

all: figures paper

clean:
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.blg paper/figs/*.aux
