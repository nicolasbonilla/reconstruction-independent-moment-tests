# Draft of 2026-09-02 (under revision)

This folder is the 2026-09-02 draft (30 pp). It is not on arXiv or in a journal and is being revised.
Known issues being corrected:

- the sealed primary endpoint (TPR 53.7% [46.6, 60.7] at FPR 0.9%) is not reported;
- the hardware sampling-noise analysis is being redone (the analytic σ understates the spread);
- no spectral moment is estimated on hardware — only the sampled support is device-derived;
- "TREX" is a misnomer (measurement twirling; no readout extinction applied);
- the n_k gate-negative outcome is not reported;
- the "necessary∘sufficient" composition wording and the m₂ ≤ 573.7 bound are being revised;
- the fig_sqw raster extent is misplaced.

The code and data in this repository are the reference.

`make paper` compiles this draft (`main.tex` with the committed `main.bbl`; no BibTeX). The source of
`figs/fig_circuit_qtk.pdf` is `figs/src/fig_circuit_qtk.tex` (`pdflatex fig_circuit_qtk.tex`). Files no
figure reads are in `figs/_superseded/`.
