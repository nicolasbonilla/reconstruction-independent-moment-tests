# Reproducibility Makefile.  Requires: `pip install -r requirements.txt`
# (qiskit for `make device`; nbconvert + ipykernel for `make reproduce`; pdflatex for `make paper`).
# Usage:  make verify | cache | figures | intervals | blind | theory | gauss | revision | device
#         | reproduce | paper | all | clean
# Every Python step runs in its own process (one recipe line = one shell). blind_addenda.py rebuilds the sealed
# truncations from a stored tie order; only its ARPACK diagnostic depends on process state.

PYTHON ?= python
LATEX  ?= pdflatex

.PHONY: all verify cache figures intervals blind theory gauss revision device reproduce paper clean help

help:
	@echo "make verify    - fast smoke test: the moment identities + the joint battery catching a truncation the lone first moment misses (seconds)"
	@echo "make cache     - rebuild the L=12 Lanczos caches in src/cache/ if missing (they are committed; make -B cache forces)"
	@echo "make figures   - regenerate the .dat fragments written by the export_* scripts and the figure-data scripts listed below"
	@echo "make intervals - interval-moment closure, joint battery (every-d scan, tie-break sensitivity) and Monte Carlo; deterministic order (~3 min)"
	@echo "make blind     - re-score the sealed blinded record, the post hoc addenda, and the fig_separating recount (~15 s)"
	@echo "make theory    - all keys of data/2026-09-27_theory_numerics.json (small_checks.py --all, ~15 min)"
	@echo "make gauss     - the fig_gausslaw data (export_gauss_state.py, ~25 s)"
	@echo "make revision  - intervals + blind + theory + gauss: every 2026-09-27 (phase A) step that needs only numpy/scipy/mpmath (~20 min)"
	@echo "make device    - reference Monte Carlo for the device coverage residual + fig_device data (needs qiskit; 7-20 min)"
	@echo "make reproduce - the master notebook end to end, phase-A scripts included (~45 min; needs nbconvert + ipykernel)"
	@echo "make paper     - compile the manuscript (paper/main.tex -> paper/main.pdf)"
	@echo "make all       - figures + paper"
	@echo "make clean     - remove LaTeX aux files (keeps main.pdf)"

# ---- fast reproducibility check (numpy/scipy only, seconds) ----
# Recomputes the current-probe moments m0,m1,m2, the shot-budget interval, and shows the joint
# (m0,m1,m2)+Hankel battery catching (via m2) a truncation the lone first moment passes (auto-scan: d=97).
# verify.py uses the deterministic lexsort key of the interval scripts (|dm1| = 0.186, |dm2| = 2.729 at d=97).
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

# ---- regenerate the data-driven figure fragments ----
# export_christoffel_dat.py also refreshes keys R9/m1/m2 of data/2026-09-27_theory_numerics.json;
# interval_moment_mc.py also writes data/2026-09-27_interval_moment_mc.json;
# export_gauss_state.py also writes data/2026-09-27_gauss_state.json;
# separating_counts.py also writes data/2026-09-27_separating_counts.json;
# export_device_dat.py only reformats the committed data/2026-09-27_delta0_reference_mc.json (see `make device`).
# The demonstrations that also write .dat files (run_collective_sumrule_falsifier.py, run_inverse_moment_falsifier.py,
# run_teeth_shared.py) run in the notebook (`make reproduce`); run_inverse_moment_falsifier.py without --r8 rewrites
# data/2026-08-28_inverse_moment_falsifier.json. export_heron_dat.py reformats data/heron_spectral.json (fig_heron).
figures: cache
	cd src && $(PYTHON) export_akw_dat.py
	cd src && $(PYTHON) export_sqw_dat.py
	cd src && $(PYTHON) check_sqw_extent.py
	cd src && $(PYTHON) export_teeth_dat.py
	cd src && $(PYTHON) export_christoffel_dat.py
	cd src && $(PYTHON) export_bracketing_dat.py
	cd src && $(PYTHON) export_momentcone_dat.py
	cd src && $(PYTHON) separating_counts.py
	cd src && $(PYTHON) interval_moment_mc.py
	cd src && $(PYTHON) export_gauss_state.py
	cd src && $(PYTHON) export_device_dat.py
	cd src && $(PYTHON) export_heron_dat.py
	@echo "OK: native pgfplots .dat fragments regenerated; 'git status paper/figs data' shows any change."

# ---- 2026-09-27 (phase A) recomputations: new dated outputs, no committed record overwritten ----
# Interval scripts: deterministic truncation order np.lexsort((index, -round(|psi0|^2, 12))).
intervals:
	cd src && $(PYTHON) interval_moment.py
	cd src && $(PYTHON) interval_battery.py
	cd src && $(PYTHON) interval_moment_mc.py

# Blinded test: blind_score.py re-writes data/2026-08-24_blind_harness_score.json from the sealed files
# (identical content); blind_addenda.py -> data/2026-09-27_blind_addenda.json (post hoc; part c rebuilds from
# the stored tie order and aborts without writing if the rebuild does not reproduce the sealed record);
# separating_counts.py -> paper/figs/sep_*.dat + data/2026-09-27_separating_counts.json.
blind:
	cd src && $(PYTHON) blind_score.py
	cd src && $(PYTHON) blind_addenda.py
	cd src && $(PYTHON) separating_counts.py

# Theory numerics: runs within_sector_lp.py, christoffel_tolerance_lp.py, necessary_sufficient_composition.py,
# estimator_form_offeigenstate.py, run_inverse_moment_falsifier.py --r8, export_christoffel_dat.py, then its
# own keys; each script can also be run alone and replaces only its own keys.
theory:
	cd src && $(PYTHON) small_checks.py --all

gauss:
	cd src && $(PYTHON) export_gauss_state.py

revision: intervals blind theory gauss
	@echo "OK: phase-A outputs regenerated (data/2026-09-27_*.json, paper/figs/*.dat); 'git status' shows any change."

# Reference Monte Carlo for Delta_0 (imports hardware_matched_job_L8.py; qiskit Statevector; no QPU).
device:
	cd src && $(PYTHON) delta0_reference_mc.py
	cd src && $(PYTHON) export_device_dat.py

# ---- run the committed generators, phase-A steps included, in one narrated pass ----
# The executed copy goes to out/ (gitignored), so the committed notebook keeps no outputs.
reproduce:
	$(PYTHON) -m nbconvert --to notebook --execute --ExecutePreprocessor.timeout=-1 \
	          --output-dir out --output 00_Reproduce_Everything.executed notebooks/00_Reproduce_Everything.ipynb
	@echo "OK: notebook executed -> out/00_Reproduce_Everything.executed.ipynb ('make device' is separate: it needs qiskit)."

# ---- compile the manuscript (uses the committed main.bbl; no bibtex needed) ----
paper:
	cd paper && $(LATEX) -interaction=nonstopmode main.tex >/dev/null && \
	            $(LATEX) -interaction=nonstopmode main.tex >/dev/null
	@echo "OK: paper/main.pdf built (text revised 2026-09-27; see paper/README.md for what is still open)."

all: figures paper

clean:
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.blg paper/figs/*.aux
