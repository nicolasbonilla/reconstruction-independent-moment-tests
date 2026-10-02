# Reproducibility Makefile.  Requires: `pip install -r requirements.txt`
# (qiskit for `make device` and `make hardware`; physics-tenpy for `make dmrg`; nbconvert + ipykernel for
# `make reproduce`; pdflatex + bibtex for `make paper`).
# Usage:  make verify | cache | figures | intervals | blind | theory | gauss | revision | device | hardware
#         | dmrg | nkrecord | nk-sealed-dir | reproduce | paper | all | clean
# Every Python step runs in its own process (one recipe line = one shell). blind_addenda.py rebuilds the sealed
# truncations from a stored tie order; only its ARPACK diagnostic depends on process state.
# After a run, `git diff --ignore-cr-at-eol --stat data paper/figs` shows which committed outputs moved
# (runtime_s, run_date and environment fields change on every run; numbers should not).

PYTHON ?= python
LATEX  ?= pdflatex
BIBTEX ?= bibtex

.PHONY: all verify cache figures intervals blind theory gauss revision device hardware dmrg nkrecord \
        nk-sealed-dir reproduce paper clean help

help:
	@echo "make verify        - smoke test + record check: the L=6 moments, thresholds, the d=97 catch, and the numbers main.tex prints from five records (seconds)"
	@echo "make cache         - rebuild the L=12 Lanczos caches in src/cache/ if missing (they are committed; make -B cache forces)"
	@echo "make figures       - regenerate the .dat fragments written by the export_* scripts and the figure-data scripts listed below"
	@echo "make intervals     - interval-moment closure, joint battery (every-d scan, tie-break sensitivity) and Monte Carlo; deterministic order (~3 min)"
	@echo "make blind         - re-score the sealed blinded record, the post hoc addenda, and the fig_separating recount (~15 s)"
	@echo "make theory        - all keys of data/2026-09-27_theory_numerics.json, incl. full_support_rr.py (small_checks.py --all, ~15 min)"
	@echo "make gauss         - the fig_gausslaw data (export_gauss_state.py, ~25 s)"
	@echo "make revision      - intervals + blind + theory + gauss: every step that needs only numpy/scipy/mpmath (~20 min)"
	@echo "make device        - reference Monte Carlo for the device coverage residual + fig_device data (needs qiskit; 7-30 min)"
	@echo "make hardware      - CZ counts of the deployed sampler circuits, the 295-CZ noise forecast and G-surface, the 2026-08-21 job reanalysis (qiskit, qiskit-aer, qiskit-ibm-runtime; no QPU)"
	@echo "make dmrg          - the L=24 DMRG code-scalability check (needs physics-tenpy; ~50 min)"
	@echo "make nkrecord      - the n_k device-attempt record, recomputed offline from its stored calibration tables (seconds)"
	@echo "make nk-sealed-dir - create 06_results/ at the repository root, where the five manifest-hashed nk_* scripts write"
	@echo "make reproduce     - the master notebook end to end (~45-60 min; needs nbconvert + ipykernel)"
	@echo "make paper         - compile the manuscript: pdflatex, bibtex, pdflatex, pdflatex (paper/main.tex -> paper/main.pdf)"
	@echo "make all           - figures + paper"
	@echo "make clean         - remove LaTeX aux files (keeps main.pdf and main.bbl)"

# ---- fast reproducibility check (numpy/scipy only, seconds) ----
# Recomputes the current-probe moments m0,m1,m2, the shot-budget interval, and shows the joint
# (m0,m1,m2)+Hankel battery catching (via m2) a truncation the lone first moment passes (auto-scan: d=97,
# deterministic lexsort order: |dm1| = 0.186, |dm2| = 2.729). Its record check matches the numbers main.tex
# prints from five records against their JSON keys; all other numbers: docs/REPRODUCE.md, Key numbers.
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
# data/2026-08-28_inverse_moment_falsifier.json. export_heron_dat.py reformats data/heron_spectral.json into
# heron_{exact,hw}.dat (read by figs/fig_heron.tex, which the manuscript no longer inputs).
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
	@echo "OK: native pgfplots .dat fragments regenerated; 'git diff --ignore-cr-at-eol --stat paper/figs data' shows any change."

# ---- recomputations with numpy/scipy/mpmath only: dated outputs, no sealed record overwritten ----
# Interval scripts: deterministic truncation order np.lexsort((index, -round(|psi0|^2, 12))).
intervals:
	cd src && $(PYTHON) interval_moment.py
	cd src && $(PYTHON) interval_battery.py
	cd src && $(PYTHON) interval_moment_mc.py

# Blinded test: blind_score.py re-writes data/2026-08-24_blind_harness_score.json from the sealed files
# (identical content); blind_addenda.py -> data/2026-09-27_blind_addenda.json (post hoc: frozen rule, sealed-text
# noise scale, variance inflation, rebuild of the truncations from the stored tie order -- it aborts without
# writing if the rebuild does not reproduce the sealed record); separating_counts.py -> paper/figs/sep_*.dat +
# data/2026-09-27_separating_counts.json.
blind:
	cd src && $(PYTHON) blind_score.py
	cd src && $(PYTHON) blind_addenda.py
	cd src && $(PYTHON) separating_counts.py

# Theory numerics: runs within_sector_lp.py, christoffel_tolerance_lp.py, necessary_sufficient_composition.py,
# estimator_form_offeigenstate.py, run_inverse_moment_falsifier.py --r8, export_christoffel_dat.py and
# full_support_rr.py (-> data/2026-10-02_full_support_rr.json), then its own keys; each script can also be run
# alone and replaces only its own keys.
theory:
	cd src && $(PYTHON) small_checks.py --all

gauss:
	cd src && $(PYTHON) export_gauss_state.py

revision: intervals blind theory gauss
	@echo "OK: outputs regenerated (data/2026-09-27_*.json, data/2026-10-02_full_support_rr.json, paper/figs/*.dat)."

# Reference Monte Carlo for Delta_0 (imports hardware_matched_job_L8.py; qiskit Statevector; no QPU).
device:
	cd src && $(PYTHON) delta0_reference_mc.py
	cd src && $(PYTHON) export_device_dat.py

# Device-side records, no QPU: sampler_cz_count.py re-counts the deposited ISA circuits offline (FakeFez
# transpilation needs qiskit-ibm-runtime); the forecast and the G-surface read its median (295 CZ) and need
# qiskit-aer; gate_tp_analysis.py re-analyses data/heron_gatetp_2026-08-21/ (9-12 min).
hardware:
	cd src && $(PYTHON) sampler_cz_count.py
	cd src && $(PYTHON) offdiag_gsurface.py
	cd src && $(PYTHON) offdiag_noise_forecast.py
	cd src && $(PYTHON) gate_tp_analysis.py

# L=24 DMRG code-scalability check (physics-tenpy); writes data/2026-10-02_beyond_ed_dmrg_catch.json.
dmrg:
	cd src && $(PYTHON) beyond_ed_dmrg_catch.py

# The n_k device-attempt record, recomputed from the calibration tables stored in it (no network).
nkrecord:
	$(PYTHON) src/nk_record/nk_attempt_record.py --offline

# The five manifest-hashed nk_* scripts (kept byte-for-byte as sealed) write to <repository root>/06_results/
# (src/../06_results, gitignored) and do not create it. Create it first, then run them, e.g. `python src/nk_falsifier.py`.
nk-sealed-dir:
	mkdir -p 06_results

# ---- run the committed generators in one narrated pass ----
# The executed copy goes to out/ (gitignored), so the committed notebook keeps no outputs.
reproduce:
	$(PYTHON) -m nbconvert --to notebook --execute --ExecutePreprocessor.timeout=-1 \
	          --output-dir out --output 00_Reproduce_Everything.executed notebooks/00_Reproduce_Everything.ipynb
	@echo "OK: notebook executed -> out/00_Reproduce_Everything.executed.ipynb ('make device', 'make hardware' and 'make dmrg' are separate)."

# ---- compile the manuscript: pdflatex, bibtex, pdflatex, pdflatex (refreshes main.bbl from refs.bib) ----
paper:
	cd paper && $(LATEX) -interaction=nonstopmode main.tex >/dev/null; \
	            $(BIBTEX) main >/dev/null && \
	            $(LATEX) -interaction=nonstopmode main.tex >/dev/null; \
	            $(LATEX) -interaction=nonstopmode main.tex >/dev/null; \
	            if grep -q "Rerun to get\|Label(s) may have changed" main.log; then \
	              $(LATEX) -interaction=nonstopmode main.tex >/dev/null; fi; \
	            if grep -E "^! |Citation .* undefined|Reference .* undefined" main.log; then exit 1; \
	            else echo "OK: paper/main.pdf built; no LaTeX error or undefined reference in main.log."; fi

all: figures paper

clean:
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.blg paper/*.toc paper/figs/*.aux
