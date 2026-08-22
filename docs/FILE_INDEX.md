# File index

Every file in this repository, described.

## `src/` — engine, screen, exporters
| File | Role |
|---|---|
| `verify.py` | Fast smoke test: moments + shot-budget interval + the joint battery closing `d=98`. |
| `spectral_lanczos.py` | Exact sector Lanczos + Haydock engine: ground state, `A(k,ω)`, `S(q,ω)`, `S^zz`, current probe, and the `run_teeth` truncation sweep. |
| `interval_moment.py` | Shot-budget interval-moment forecast at the real `ibm_fez` 50k-shot budget (local-estimator variance + swept bias). |
| `interval_battery.py` | Joint `(m₀,m₁,m₂)` + Hankel–Stieltjes feasibility battery; closes the single-moment blind spot. |
| `run_heron_screen.py` | Runs the screen on the real IBM Heron SQD spectral data. |
| `run_sumrule_falsifier.py` | The Drude weight-misplacement demonstration (`m₀` blind, `m₁` fires). |
| `run_gausslaw_falsifier.py` | The cross-domain U(1) lattice-Schwinger Gauss-law falsifier. |
| `run_falsifier_teeth.py` | The teeth: independent estimator vs the circular control. |
| `run_moment_bound_theorem.py` | The Christoffel miss-distance certificate and the Gauss–Radau bracketing. |
| `export_*.py` | One data-driven generator per figure fragment (writes `paper/figs/*.dat`). |
| `make_akw_figure.py`, `make_sqw_figure.py`, `make_paper_figs.py` | Raster-field + native-axes heatmap builders. |
| `paper_style.py` | Shared plotting identity for the raster fields. |

## `data/` — authoritative plotted values
`*.dat` — every plotted number for the native pgfplots figures. `heron_spectral.json` — the real
`ibm_fez` reconstructed spectra (companion study). `2026-08-22_interval_moment_closure.json` — the
committed interval-moment forecast output.

## `paper/` — the manuscript
`main.tex` (master), `figstyle.tex` (shared figure identity), `main.bbl` (frozen bibliography),
`refs.bib`, `main.pdf` (compiled, 18 pp), `arxiv-submission.tar.gz` (upload-ready bundle), and
`figs/` (native `.tex` fragments + `.dat` + raster `.png`).

## `notebooks/`
`00_Reproduce_Everything.ipynb` — the master notebook: every figure and number, narrated end to end.

## `docs/`
`REPRODUCE.md` (figure/number → script → command), `FILE_INDEX.md` (this file),
`ARXIV_SUBMISSION.md` (the verified upload bundle + metadata + steps), `img/` (README thumbnails).

## Root
`README.md`, `CITATION.cff`, `LICENSE` (MIT code + CC-BY-4.0 paper), `requirements.txt`, `Makefile`.
