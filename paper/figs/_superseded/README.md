# Orphan figure files (not read by any figure of the draft)

Moved here 2026-09-26; nothing was deleted. No `\input`/`table` in `paper/main.tex` or `paper/figs/*.tex`
reads these files.

- `falsifier_spectrum.dat`, `falsifier_sweep.dat` — data of the removed `fig_falsifier`; generator `src/_superseded/export_falsifier_dat.py`.
- `gauss.dat` — the bit-flip sweep of `⟨ΣG²⟩` (the content of `data/2026-08-18_gausslaw_falsifier.json`); `fig_gausslaw` plots `gauss_state.dat` instead, whose generator is not yet committed (known gap).
- `gauss_spectral.dat` — output of a circular Gauss-law "spectral screen" (clean-state moments against a reconstruction of the contaminated state) that the 2026-08-26 integrity pass dropped; its generator is not in this repository.
- `heron_sweep.dat` — a weight-misplacement sweep on the companion's Heron reconstruction (`f`, `r0`, `r1`), i.e. the reconstruction-vs-reconstruction lineage of `src/_superseded/run_heron_screen.py`; no committed generator; `fig_heron` does not read it.
- `sep_exemplar.dat` — an exemplar point set for `fig_separating` that the figure does not read.
- `fig_circuit.tex` — an earlier native-TikZ circuit; the draft includes `fig_circuit_qtk.pdf`, whose source is `paper/figs/src/fig_circuit_qtk.tex`.
