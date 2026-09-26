# Reproduce the numbers and figures

Every non-hardware result in the paper is an exact classical computation. The core scripts need
`numpy`/`scipy`/`mpmath`; the noise-model forecasts and the re-analysis of the retained `ibm_fez` counts
also need `qiskit`/`qiskit-aer`, and the `L=24` DMRG check needs `physics-tenpy` (see
`requirements.txt`). Run scripts from the repository root (`python src/<script>.py`); their paths are
relative to `src/`, so they write their JSON into `data/` and their plotted data into `paper/figs/*.dat`,
and `git status` after a run shows whether a committed number moved. The items that cannot be
regenerated yet are listed under [Known gaps](#known-gaps).

## Fastest check (seconds)

```bash
python src/verify.py        # or: make verify
```
Recomputes the current-probe moments `m₀, m₁, m₂` of the doped `L=6, U/t=4` ring, the shot-budget
interval at the real `ibm_fez` 50k-shot budget, and shows the joint `(m₀,m₁,m₂)`+Hankel battery
**catching**, via `m₂`, a determinant truncation that the lone first moment misses (the truncation is
found by an auto-calibrated scan, currently `d=97`). Expected: **ALL CHECKS PASS**.

## Figure → engine → data

| Figure (paper) | What it shows | Engine / script | Data fragment(s) |
|---|---|---|---|
| `fig_hero` | the reconstruction-independent moment loop (schematic) | native TikZ | — |
| `fig_momentcone` | a moment sequence refutable by geometry | `export_momentcone_dat.py` (≈2 min) | `momentcone_point.dat` |
| `fig_akw` | `A(k,ω)` Mott map + per-`k` screen | `spectral_lanczos` → `export_akw_dat.py`; needs `src/cache/akw_L12.npz` (committed / `make cache`) | `akw_*.dat`, `akw_true/wrong.png` |
| `fig_sqw` | `S(q,ω)` gapped / `S^zz` gapless | `spectral_lanczos` → `export_sqw_dat.py`; needs `src/cache/sqw_L12.npz` (committed / `make cache`) | `sqw_*.dat`, `dcp_*.dat`, `sqw_charge/spin.png` |
| `fig_collective_screen` | the screen transported to the collective channels | `run_collective_sumrule_falsifier.py` | `collective_screen_*.dat` |
| `fig_inverse_falsifier` | the negative-order (`m₋₁`, f-sum) test | `run_inverse_moment_falsifier.py` | `inverse_falsifier_*.dat` |
| `fig_teeth` | independence vs a circular control | `spectral_lanczos.run_teeth` → `export_teeth_dat.py`; needs `src/cache/teeth_L12.npz` (committed / `make cache`) | `teeth.dat` |
| `fig_teeth_shared` | the shared-state leak that still fires the screen | `run_teeth_shared.py` | `teeth_shared.dat` |
| `fig_separating` | diagonal vs off-diagonal complementary reach (blinded) | `separating_counts.py` from the sealed `data/blind_*.json` (`.dat` export pending) | `sep_*.dat` |
| `fig_device` | device-side coverage-leak curve (real `ibm_fez`, `L=8`) | `hardware_matched_job_L8.py` (acquisition), `bootstrap_device_curve.py`, `aer_noiseless_baseline.py`; `retention_matched_control.py` → `data/2026-09-05_retention_matched_control.json` | coordinates typed inline in `fig_device.tex` (exporter pending) |
| `fig_heron` | executability precursor on real IBM Heron data | data: `paper/figs/heron_{exact,hw}.dat` from `data/heron_spectral.json`; exporter pending | `heron_*.dat` |
| `fig_circuit` | schematic of the simulated bond-basis estimator (not run on hardware) | `paper/figs/src/fig_circuit_qtk.tex` → `fig_circuit_qtk.pdf` (`pdflatex fig_circuit_qtk.tex`) | — |
| `fig_bracketing` + `fig_christoffel` | Gauss–Radau bracketing + the Christoffel width `Wₙ(t)` | `export_bracketing_dat.py`, `export_christoffel_dat.py` (≈2 min each); the two-sided window bound: `markov_krein_window.py` | `bracket_*.dat`, `christoffel_*.dat` |
| `fig_momentmc` | Monte-Carlo interval-moment test at the `ibm_fez` budget | `interval_moment_mc.py` | `momentmc_*.dat` |
| `fig_gausslaw` | U(1) Gauss-law cross-domain state-level check | data: `paper/figs/gauss_state.dat` — generator not yet committed (known gap); `run_gausslaw_falsifier.py` writes the bit-flip sweep JSON, which the figure does not plot | `gauss_state.dat` |

## Key numbers → script

| Number | Value | Where |
|---|---|---|
| current-probe moments `(m₀,m₁,m₂)` (L=6, U/t=4) | `0.5542, 4.2271, 33.3887` | `verify.py`, `interval_moment.py` |
| shot-budget 95% interval `z·δ₁` (ibm_fez, 2% bias) | `0.229` (5.4% of m₁) | `verify.py`, `interval_moment.py` |
| joint battery closes the lone-`m₁` blind spot | `|Δm₂| > z·δ₂` while `|Δm₁| ≤ z·δ₁` (auto-scan: `d=97`; the committed `d=98` case is under revision) | `verify.py`, `interval_battery.py`, `within_sector_control.py` |
| **sealed primary endpoint** | TPR `101/188 = 53.7%` [46.6, 60.7], FPR `1/112 = 0.9%` [0.2, 4.9] | `blind_score.py` → `data/2026-08-24_blind_harness_score.json` |
| blinded pre-registered per-class rates | truncation `56/56`; spurious `45/61` (`25/27`+`20/34`); FP `1/112`; Krylov `0/71` (one-node `0/15`) | `blind_score.py`, `data/blind_verdicts.json` (seal `data/prereg.sha256`) |
| `n_k` amendment gate | **DENIED (gate-negative)**: R1ro0.018 UCB `0.02867` > `0.02667` | `nk_stage0_gate_v4.py`, `nk_stage0_gate_v4_combine.py` → `data/2026-09-05_nk_v4_verdict.json` |
| retention-matched control (noiseless Δ₀ at the device's kept-shot count) | `0.203 / 0.248 / 0.329 / 0.418` vs device `0.0044 / 0.032 / 0.123 / 0.373` | `retention_matched_control.py` → `data/2026-09-05_retention_matched_control.json` |
| necessary∘sufficient composition (`m₂` bracket, L=6 U/t=8) (under revision) | `123.3 ≤ m₂ ≤ 573.7` | `necessary_sufficient_composition.py`, `markov_krein_window.py` |
| off-diagonal device bias floor (ibm_fez depth) | `≈ 37–115 %` of the signal | `offdiag_noise_forecast.py`, `offdiag_gsurface.py` |
| `n_k` momentum-distribution forecast (FT-only bias) | `≈ 9 %` of the Fermi step | `nk_falsifier.py` |
| single-particle Mott edges / gap (L=12) | `μ± = ±2.484t`, `Δ ≈ 4.97t` | `spectral_lanczos` (`A(k,ω)`) |
| charge onset `Δc` (L=12) | `≈ 5.71t` | `spectral_lanczos` (structure factors) |

## Blinded pre-registration

The blinded test is sealed: `data/prereg.json` + its SHA-256 `data/prereg.sha256`
(`0e4c3368…`) fix the battery, thresholds, error catalog, and seed before any instance exists.
`blind_generate.py` produces the instances, `blind_classify.py`/`blind_score.py` apply the frozen
battery, and `data/blind_{instances_public,labels_sealed,verdicts}.json` + `separating_demonstration.json`
hold the record. The generator and classifier both re-verify the seal and refuse to run if it is broken;
the sealed files are read-only, so the writers refuse to overwrite them (run in a scratch copy to
regenerate). `separating_counts.py` recounts `fig_separating` from the sealed files alone.

## Sealed records: how to check them

`.gitattributes` keeps these files byte-exact on every platform (no CRLF conversion).

```bash
sha256sum data/prereg.json                   # = data/prereg.sha256                (0e4c3368…2792)
sha256sum data/manifest_nk_device_v1.json    # = data/manifest_nk_device_v1.sha256 (496203db…ae80)
cd data && sha256sum -c 2026-09-05_nk_v4_verdict.json.sha256   # 784927f0…325a
```

The manifest also records the SHA-256 of the analysis code and gate records it was sealed against
(`analysis_code_sha256`, `gate_record_sha256`). Those eight files (`src/nk_stage0_gate{,_v2,_v3}.py`,
`src/nk_stage0_v3R2.py`, `src/nk_falsifier.py`, `src/spectral_lanczos.py`,
`data/2026-09-02_nk_stage0_gate_v3.json`, `data/2026-09-02_nk_stage0_v3R2.json`) are kept byte-for-byte
as sealed, including their original line endings, so every recorded hash verifies; for the same reason
the five `nk_*` scripts among them still write to `./06_results/` (gitignored), not to `data/`.

## Hardware

Two device results use real `ibm_fez` (IBM Heron) data: the `L=6` full-sector executability precursor
(`fig_heron`, reprocessing the companion study [arXiv:2608.16436](https://arxiv.org/abs/2608.16436))
and a purpose-built `L=8` retained-count coverage-leak job (`fig_device`). Only the diagonal support
set is device-derived: the `m0_hat`/`m1_hat` fields of `data/heron_counts_matched_L8_*.json` are exact
classical values, so no spectral moment is estimated on hardware. The retention-matched control
(`retention_matched_control.py`) re-runs the noiseless circuits at the device's kept-shot count
(≈23.6% of raw shots survive post-selection on the device, 100% in the noiseless simulation). The
discriminating off-diagonal moments are evaluated in simulation and shown presently infeasible on
hardware. The on-device `n_k` experiment was pre-registered (manifest v1 sealed 2026-09-02), its single
permitted amendment was denied by its sealed gate on 2026-09-05, and it was never run
(`nk_stage2_device_job.py` is the manifest-v1 job; `nk_stage2_device_job_v2_amendment.py` is the denied
amendment's job). Re-acquiring device data needs an IBM Quantum account (use your own credentials; never
commit a token).

## Known gaps

Pending generators and recomputations (being fixed with the revision; do not treat these as reproducible):

- `paper/figs/gauss_state.dat` (`fig_gausslaw`) has no committed generator.
- No exporter yet for `paper/figs/heron_{exact,hw}.dat` (`fig_heron`, from `data/heron_spectral.json`) or
  for the inline `fig_device` coordinates (from the retained counts).
- The `sep_*.dat` export from `separating_counts.py` (`fig_separating`) is pending.
- The hardware sampling-noise analysis (a reference Monte Carlo for Δ₀) is being redone; the draft's
  analytic σ understates the spread.
- The `m₂` bracket (`123.3 ≤ m₂ ≤ 573.7`) is being revised.
- The `d≈98` single-moment miss: a fresh run of `interval_moment.py` gives gap `0.2542` at `d=98` against
  `0.2625` committed (the truncation depends on tie-breaking); `verify.py` finds `d=97`;
  `interval_battery.py` prints a hard-coded "d=98" label.
- Figure fixes (including the `fig_sqw` raster extent) and regenerated README thumbnails.
