# Reproduce every number and figure

Every non-hardware result in the paper is an exact classical computation reproducible from this
repository with `numpy`/`scipy` alone (a few noise-model forecasts also use `qiskit-aer`). Run scripts
from `src/` (`cd src && python <script>.py`); the native pgfplots figures read the committed data
fragments in `paper/figs/*.dat`, which the `export_*`/`run_*` scripts regenerate.

## Fastest check (seconds)

```bash
python src/verify.py        # or: make verify
```
Recomputes the current-probe moments `m₀, m₁, m₂` of the doped `L=6, U/t=4` ring, the shot-budget
interval at the real `ibm_fez` 50k-shot budget, and shows the joint `(m₀,m₁,m₂)`+Hankel battery
**catching**, via `m₂`, a determinant truncation that the lone first moment misses (the truncation is
found by an auto-calibrated scan, so the demonstration tracks the current interval). Expected:
**ALL CHECKS PASS**.

## Figure → engine → data

| Figure (paper) | What it shows | Engine / script | Data fragment(s) |
|---|---|---|---|
| `fig_hero` | the reconstruction-independent moment loop (schematic) | native TikZ | — |
| `fig_momentcone` | a moment sequence refutable by geometry | `export_momentcone_dat.py` | `momentcone_point.dat` |
| `fig_akw` | `A(k,ω)` Mott map + per-`k` screen | `spectral_lanczos` → `export_akw_dat.py` | `akw_*.dat`, `akw_true/wrong.png` |
| `fig_sqw` | `S(q,ω)` gapped / `S^zz` gapless | `spectral_lanczos` → `export_sqw_dat.py` | `sqw_*.dat`, `dcp_*.dat` |
| `fig_collective_screen` | the screen transported to the collective channels | `run_collective_sumrule_falsifier.py` | `collective_screen_*.dat` |
| `fig_inverse_falsifier` | the negative-order (`m₋₁`, f-sum) test | `run_inverse_moment_falsifier.py` | `inverse_falsifier_*.dat` |
| `fig_teeth` | independence vs a circular control | `run_falsifier_teeth.py` → `export_teeth_dat.py` | `teeth.dat` |
| `fig_teeth_shared` | the shared-state leak that still fires the screen | `run_teeth_shared.py` | `teeth_shared.dat` |
| `fig_separating` | diagonal vs off-diagonal complementary reach (blinded) | sealed blind harness (`blind_*.py`) → `separating_demonstration.json` | `sep_*.dat` |
| `fig_device` | device-side coverage-leak curve (real `ibm_fez`, `L=8`) | `hardware_matched_job_L8.py`, `bootstrap_device_curve.py` | inline (retained counts) |
| `fig_heron` | executability precursor on real IBM Heron data | `run_heron_screen.py` | `heron_*.dat` |
| `fig_circuit` | the device-side moment-estimation protocol (schematic) | native TikZ / PDF | — |
| `fig_bracketing` + `fig_christoffel` | Gauss–Radau bracketing + certified `Wₙ(t)` | `markov_krein_window.py`, `run_moment_bound_theorem.py` → `export_bracketing_dat.py`, `export_christoffel_dat.py` | `bracket_*.dat`, `christoffel_*.dat` |
| `fig_momentmc` | Monte-Carlo interval-moment test at the `ibm_fez` budget | `interval_moment_mc.py` | `momentmc_*.dat` |
| `fig_gausslaw` | U(1) Gauss-law cross-domain falsifier | `run_gausslaw_falsifier.py` | `gauss*.dat` |

## Key numbers → script

| Number | Value | Where |
|---|---|---|
| current-probe moments `(m₀,m₁,m₂)` (L=6, U/t=4) | `0.5542, 4.2271, 33.3887` | `verify.py`, `interval_moment.py` |
| shot-budget 95% interval `z·δ₁` (ibm_fez, 2% bias) | `0.229` (5.4% of m₁) | `verify.py`, `interval_moment.py` |
| joint battery closes the lone-`m₁` blind spot | `|Δm₂| > z·δ₂` while `|Δm₁| ≤ z·δ₁` | `verify.py`, `interval_battery.py`, `within_sector_control.py` |
| blinded pre-registered rates | truncation `56/56`; spurious `45/61` (`25/27`+`20/34`); FP `1/112`; Krylov `0/71` | `blind_score.py`, `data/blind_verdicts.json` (seal `data/prereg.sha256`) |
| necessary∘sufficient composition (`m₂` bracket, L=6 U/t=8) | `123.3 ≤ m₂ ≤ 573.7` | `necessary_sufficient_composition.py`, `markov_krein_window.py` |
| off-diagonal device bias floor (ibm_fez depth) | `≈ 37–115 %` of the signal | `offdiag_noise_forecast.py`, `offdiag_gsurface.py` |
| `n_k` momentum-distribution forecast (FT-only bias) | `≈ 9 %` of the Fermi step | `nk_falsifier.py` |
| single-particle Mott edges / gap (L=12) | `μ± = ±2.484t`, `Δ ≈ 4.97t` | `spectral_lanczos` (`A(k,ω)`) |
| charge onset `Δc` (L=12) | `≈ 5.71t` | `spectral_lanczos` (structure factors) |

## Blinded pre-registration

The blinded test is sealed: `data/prereg.json` + its SHA-256 `data/prereg.sha256`
(`0e4c3368…`) fix the battery, thresholds, error catalog, and seed before any instance exists.
`blind_generate.py` produces the instances, `blind_classify.py`/`blind_score.py` apply the frozen
battery, and `data/blind_{instances_public,labels_sealed,verdicts}.json` + `separating_demonstration.json`
hold the record. The generator and classifier both re-verify the seal and refuse to run if it is broken.

## Hardware

Two device results use real `ibm_fez` (IBM Heron) data: the `L=6` full-sector executability precursor
(`fig_heron`, reprocessing the companion study [arXiv:2608.16436](https://arxiv.org/abs/2608.16436))
and a purpose-built `L=8` retained-count coverage-leak job (`fig_device`). Only the diagonal support
set is device-derived; the discriminating off-diagonal moments are evaluated in simulation and shown
presently infeasible on hardware. Re-acquiring device data needs an IBM Quantum account (use your own
credentials; never commit a token); the classical screen does not.
