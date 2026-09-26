# Reproducibility Makefile.  Requires: `pip install -r requirements.txt`
# (and a TeX distribution with pdflatex for `make paper`).
# Usage:  make verify | make cache | make figures | make reproduce | make paper | make all | make clean

PYTHON ?= python
LATEX  ?= pdflatex

.PHONY: all verify cache figures reproduce paper clean help

help:
	@echo "make verify    - fast smoke test: the sum-rule identities + the interval battery catching a truncation the lone first moment misses (seconds)"
	@echo "make cache     - rebuild the L=12 Lanczos caches in src/cache/ if missing (they are committed; make -B cache forces)"
	@echo "make figures   - regenerate the native pgfplots data fragments (.dat) that have a committed generator"
	@echo "make reproduce - run the master notebook end to end (needs nbconvert + ipykernel, ~30 min; see docs/REPRODUCE.md, Known gaps)"
	@echo "make paper     - compile the 2026-09-02 draft in paper/ (paper/main.tex -> paper/main.pdf)"
	@echo "make all       - figures + paper"
	@echo "make clean     - remove LaTeX aux files (keeps main.pdf)"

# ---- fast reproducibility check (numpy/scipy only, seconds) ----
# Recomputes the current-probe moments m0,m1,m2, the shot-budget interval, and shows the joint
# (m0,m1,m2)+Hankel battery catching (via m2) a truncation the lone first moment misses.
verify:
	$(PYTHON) src/verify.py

# ---- L=12 sector-Lanczos + Haydock caches (committed; rebuilt only if missing) ----
CACHE_FILES = src/cache/akw_L12.npz src/cache/sqw_L12.npz src/cache/teeth_L12.npz
cache: $(CACHE_FILES)
src/cache/akw_L12.npz:
	cd src && $(PYTHON) spectral_lanczos.py 12 akw
src/cache/sqw_L12.npz:
	cd src && $(PYTHON) spectral_lanczos.py 12 sqw
src/cache/teeth_L12.npz:
	cd src && $(PYTHON) -c "import spectral_lanczos as s; s.run_teeth()"

# ---- regenerate the data-driven figure fragments from the exact sector-Lanczos + Haydock engine ----
# (fig_gausslaw, fig_heron, fig_device and the sep_*.dat export have no committed generator yet:
#  see docs/REPRODUCE.md, Known gaps)
figures: cache
	cd src && $(PYTHON) export_akw_dat.py && $(PYTHON) export_sqw_dat.py && \
	          $(PYTHON) export_teeth_dat.py && \
	          $(PYTHON) export_christoffel_dat.py && $(PYTHON) export_bracketing_dat.py && \
	          $(PYTHON) export_momentcone_dat.py && $(PYTHON) separating_counts.py
	@echo "OK: native pgfplots .dat fragments regenerated; 'git status paper/figs' shows any change."

# ---- run the committed generators in one narrated pass (optional cells need qiskit-aer) ----
reproduce:
	$(PYTHON) -m nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 \
	          notebooks/00_Reproduce_Everything.ipynb
	@echo "OK: full pipeline executed in notebooks/00_Reproduce_Everything.ipynb"

# ---- compile the 2026-09-02 draft (uses the committed main.bbl; no bibtex needed) ----
paper:
	cd paper && $(LATEX) -interaction=nonstopmode main.tex >/dev/null && \
	            $(LATEX) -interaction=nonstopmode main.tex >/dev/null
	@echo "OK: paper/main.pdf built."

all: figures paper

clean:
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.blg paper/figs/*.aux
