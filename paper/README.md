# Draft of 2026-09-02 (under revision)

This folder is the 2026-09-02 draft (30 pp). It is not on arXiv or in a journal and is being revised.
Known issues being corrected include (the list is not exhaustive):

- the sealed primary endpoint (TPR 53.7% [46.6, 60.7] at FPR 0.9%) is not reported;
- the hardware sampling-noise analysis is being redone (the analytic σ understates the spread);
- no spectral moment is estimated on hardware — only the sampled support is device-derived. The abstract
  and introduction say that each moment "is estimated from the same computational-basis samples that
  produced the spectrum" and that one moment is evaluated on hardware; neither holds in this work (the
  `m0_hat`/`m1_hat` values are exact classical ones, and the off-diagonal `m₁`, `m₂` would need
  rotated-basis circuits, not the computational-basis samples);
- "TREX" is a misnomer: the `ibm_fez` runs used SamplerV2 measurement twirling, Pauli gate twirling and
  dynamical decoupling, with no readout-error mitigation;
- the method text describes symmetry-restoring configuration recovery, but neither device run used it:
  the support is defined by post-selection alone. The companion's `L=6` run behind `fig_heron` was also
  post-selected in reversed qubit order, so its 300/300 coverage came from device errors (disclosed by the
  companion repository on 2026-09-26); the `L=8` run of this work reads the bits in the correct order;
- "conditioning is benign and detection is generous" (the within-sector section) overstates the power:
  `Δ₂^max = 7.1` is a worst case, and at the deployed order `m₀`–`m₂` a `K=1` redistribution moving 49% of
  the weight passes the battery (`data/2026-08-24_within_sector_control.json`), so there is no guaranteed
  power;
- the `fig_separating` caption: the recount (`src/separating_counts.py`) finds 26 instances in the shaded
  region, not 24 (`m₁` fires on 16, `m₂` alone on 10), including one determinant truncation (id 162) and
  the clean false positive (id 184); the exemplar id 222 is weight-preserving, not total-weight-changing;
- the `d=98` showcase (`|g₁| = 0.225 < τ₁`) does not reproduce from a fresh clone (`|g₁| = 0.254`, above
  the 2%-bias interval 0.229); `d=97` does (see the repository README, Known gaps);
- the n_k gate-negative outcome is not reported;
- the "necessary∘sufficient" composition wording and the m₂ ≤ 573.7 bound are being revised;
- the fig_sqw raster extent is misplaced.

The code and data in this repository are the reference.

`make paper` compiles this draft (`main.tex` with the committed `main.bbl`; no BibTeX). The source of
`figs/fig_circuit_qtk.pdf` is `figs/src/fig_circuit_qtk.tex` (`pdflatex fig_circuit_qtk.tex`). Files no
figure reads are in `figs/_superseded/`.
