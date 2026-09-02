# File index

Every file in this repository, described.

## `src/` — engine, moment tests, figure generators, blinded protocol
| File(s) | Role |
|---|---|
| `verify.py` | Fast smoke test: the current-probe moments, the shot-budget interval, and the joint battery catching (via `m₂`) a truncation the lone first moment misses. |
| `spectral_lanczos.py` | Exact sector Lanczos + Haydock engine: ground state, `A(k,ω)`, `S(q,ω)`, `S^zz`, current/density probes, truncation sweeps. |
| `interval_moment.py`, `interval_moment_mc.py`, `interval_battery.py` | Shot-budget interval-moment forecast (analytic + Monte-Carlo) and the joint `(m₀,m₁,m₂)`+Hankel–Stieltjes feasibility battery. |
| `necessary_sufficient_composition.py`, `markov_krein_window.py`, `within_sector_control.py` | The necessary∘sufficient composition, the two-sided Markov–Krein window-mass bound, and the within-sector detectability/completeness analysis. |
| `run_sumrule_falsifier.py`, `run_collective_sumrule_falsifier.py`, `run_inverse_moment_falsifier.py` | The current-response Drude, collective-channel, and negative-order (`m₋₁`/f-sum) demonstrations. |
| `run_falsifier_teeth.py`, `run_teeth_shared.py` | The teeth (independence vs a circular control) and the shared-state leak. |
| `run_gausslaw_falsifier.py` | The cross-domain U(1) lattice-Schwinger Gauss-law falsifier. |
| `run_moment_bound_theorem.py` | The Christoffel miss-distance certificate and the Gauss–Radau bracketing. |
| `nk_falsifier.py`, `joint_covariance_composition.py`, `lever1_flip_experiment.py` | The momentum-distribution `n_k` forecast, the same-sample covariance (Durbin–Wu–Hausman) analysis, and its pre-registered decision experiment. |
| `offdiag_noise_forecast.py`, `offdiag_gsurface.py`, `aer_noiseless_baseline.py`, `bootstrap_device_curve.py`, `analytic_coverage_leak.py` | The off-diagonal depolarizing-bias floor and the device-side coverage-leak forecasts. |
| `estimator_scaling.py`, `bond_moment_estimator.py` | Circuit-scaling and the reconstruction-independent bond-basis moment estimator. |
| `blind_preregister.py`, `blind_generate.py`, `blind_classify.py`, `blind_score.py` | The sealed pre-registered blinded protocol (seal → generate → classify → score). |
| `beyond_ed_dmrg_catch.py` | The `L=24` DMRG-oracle scalability check of the coverage falsifier. |
| `hardware_matched_job_L8.py`, `hardware_matched_job.py`, `hardware_blind_job.py`, `run_heron_screen.py` | Device-side job builders and the screen on real IBM Heron data (use your own credentials; never commit a token). |
| `export_*.py` | One data-driven generator per figure fragment (writes `paper/figs/*.dat`). |
| `make_akw_figure.py`, `make_sqw_figure.py`, `make_paper_figs.py`, `paper_style.py` | Raster-field + native-axes heatmap builders and the shared plotting identity. |

## `data/` — authoritative values and the sealed record
`*.json` — the committed outputs of the scripts above (interval forecasts, the blinded harness, the
`n_k`/covariance/off-diagonal/DMRG results, and the real `ibm_fez` `L=8` retained counts).
`prereg.json` + `prereg.sha256` and `blind_{instances_public,labels_sealed,verdicts}.json` +
`separating_demonstration.json` are the tamper-evident blinded pre-registration. `heron_spectral.json`
is the companion study's reconstructed spectra. The `paper/figs/*.dat` fragments hold every plotted
number for the native pgfplots figures.

## `paper/` — the manuscript
`main.tex` (master), `figstyle.tex` (shared figure identity), `main.bbl` (frozen bibliography),
`refs.bib`, `main.pdf` (compiled, 30 pp), `arxiv-submission.tar.gz` (upload-ready bundle), and
`figs/` (native `.tex` fragments + `.dat` + raster `.png`).

## `notebooks/`
`00_Reproduce_Everything.ipynb` — the master notebook: every figure and number, narrated end to end.

## `docs/`
`REPRODUCE.md` (figure/number → script → command), `FILE_INDEX.md` (this file),
`ARXIV_SUBMISSION.md` (the verified upload bundle + metadata + steps), `img/` (README thumbnails).

## Root
`README.md`, `CITATION.cff`, `LICENSE` (MIT code + CC-BY-4.0 paper), `requirements.txt`, `Makefile`.
