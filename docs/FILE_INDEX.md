# File index

Every file in this repository, described (regenerated from `git ls-files` on 2026-09-26).

## `src/` — engine, moment tests, figure generators, blinded protocol, device jobs

Scripts write their JSON into `data/` and plotted data into `paper/figs/` (paths relative to `src/`).

| File | Role |
|---|---|
| `verify.py` | Fast smoke test: the current-probe moments, the shot-budget interval, and the joint battery catching (via `m₂`) an auto-calibrated truncation the lone first moment misses. |
| `hubbard_ed.py` | Vendored exact-diagonalization utilities (Jordan–Wigner Hubbard chain; numpy/scipy only), imported by the exporters and demonstrations. |
| `spectral_lanczos.py` | Exact sector Lanczos + Haydock engine: ground state, `A(k,ω)`, `S(q,ω)`, `S^zz`, current probe, truncation sweeps; writes `src/cache/`. Hashed in the `n_k` manifest and kept byte-for-byte. |
| `interval_moment.py` | Shot-budget interval-moment forecast at the `ibm_fez` budget → `data/2026-08-22_interval_moment_closure.json`. |
| `interval_moment_mc.py` | Monte-Carlo interval-moment test (`fig_momentmc`) → `data/2026-08-24_interval_moment_mc.json`, `momentmc_*.dat`. |
| `interval_battery.py` | Joint `(m₀,m₁,m₂)`+Hankel–Stieltjes feasibility battery → `data/interval_battery_showcase.json`. |
| `within_sector_control.py` | Within-sector detectability: moment-preserving redistributions and Vandermonde conditioning → `data/2026-08-24_within_sector_control.json`. |
| `markov_krein_window.py` | Two-sided Markov–Krein window-mass bound → `data/2026-08-28_markov_krein_window.json`. |
| `necessary_sufficient_composition.py` | The `m₂` bracket from the Hankel (necessary) and support (sufficient) conditions → `data/necessary_sufficient_composition.json` (under revision). |
| `chigap_check.py` | Christoffel empty-region check on the `L=6` current measure → `data/2026-08-24_chigap_check.json`. |
| `run_sumrule_falsifier.py` | Current-response demonstration: the total-weight `m₀` sum rule is blind, `m₁` detects → `data/2026-08-18_sumrule_falsifier.json`; defines `current_operator`. |
| `run_collective_sumrule_falsifier.py` | The screen on the collective (charge/spin) channels → `data/2026-08-26_collective_sumrule_falsifier.json`, `collective_screen_*.dat`. |
| `run_inverse_moment_falsifier.py` | The negative-order `m₋₁` (optical f-sum) test → `data/2026-08-28_inverse_moment_falsifier.json`, `inverse_falsifier_*.dat`. |
| `run_falsifier_teeth.py` | The teeth at `L=6` (independent estimator vs a circular control) → `data/2026-08-18_falsifier_teeth.json`. |
| `run_teeth_shared.py` | The shared-state (device-realizable) teeth → `data/2026-08-26_teeth_shared_state.json`, `teeth_shared.dat`. |
| `run_gausslaw_falsifier.py` | U(1) quantum link model Gauss-law state-level check (bit-flip sweep) → `data/2026-08-18_gausslaw_falsifier.json`. |
| `export_akw_dat.py` | `fig_akw` fragments and rasters from `src/cache/akw_L12.npz`. |
| `export_sqw_dat.py` | `fig_sqw` fragments and rasters from `src/cache/sqw_L12.npz`. |
| `export_teeth_dat.py` | `fig_teeth` fragment from `src/cache/teeth_L12.npz`. |
| `export_bracketing_dat.py` | `fig_bracketing` fragments (Gauss–Radau bracketing). |
| `export_christoffel_dat.py` | `fig_christoffel` fragments (Christoffel width `Wₙ(t)`). |
| `export_momentcone_dat.py` | `fig_momentcone` point. |
| `blind_preregister.py` | Seals the blinded protocol (`data/prereg.json` + `.sha256`); refuses to overwrite the committed seal. |
| `blind_generate.py` | Generates the 300 blinded instances from the sealed prereg (refuses to overwrite the sealed outputs). |
| `blind_classify.py` | Applies the frozen battery blind → `data/blind_verdicts.json` (refuses to overwrite it). |
| `blind_score.py` | Unblinds and scores → `data/2026-08-24_blind_harness_score.json` (TPR/FPR with Wilson intervals, per class). |
| `separating_counts.py` | Recounts `fig_separating` from the sealed record alone (`.dat` export pending). |
| `bond_moment_estimator.py` | Bond-basis (Pauli-group) moment estimator, simulated → `data/2026-08-25_bond_moment_estimator_L{4,6}.json`. |
| `estimator_scaling.py` | Circuit-count scaling of that estimator → `data/2026-08-25_estimator_circuit_scaling.json`. |
| `verify_mm1_measurement_cost.py` | Checks that `m₋₁` needs no measurement settings beyond the `m₀`/`m₁` plan (prints only). |
| `offdiag_noise_forecast.py` | Off-diagonal depolarizing-bias floor under `ibm_fez`-anchored noise → `data/2026-08-30_offdiag_noise_forecast.json`. |
| `offdiag_gsurface.py` | That bias floor vs circuit depth → `data/2026-08-31_offdiag_gsurface.json`. |
| `joint_covariance_composition.py` | Same-sample covariance of the falsifier → `data/2026-09-01_joint_covariance.json`. |
| `lever1_flip_experiment.py` | Pre-registered decision experiment on that covariance → `data/2026-09-01_lever1_{prereg_seal,oracle}.json` (re-verifies the committed seal instead of overwriting it). |
| `beyond_ed_dmrg_catch.py` | `L=24` DMRG-oracle check of the coverage screen (needs `physics-tenpy`) → `data/2026-08-28_beyond_ed_dmrg_catch.json`. |
| `hardware_matched_job_L8.py` | The `ibm_fez` `L=8` coverage-leak job: Aer dry run (default) or QPU (`RUN=1`) → `data/heron_counts_matched_L8_*.json`. |
| `hardware_matched_job.py` | The `L=6` matched job on the companion's circuits → `data/heron_counts_matched_DRYRUN.json`. |
| `hardware_blind_job.py` | A blinded hardware job builder; never run on hardware; only its dry-run output is committed (`data/heron_counts_DRYRUN.json`). |
| `aer_noiseless_baseline.py` | Noiseless Aer run of the same `L=8` circuits → `data/aer_noiseless_baseline.json`. |
| `bootstrap_device_curve.py` | Bootstrap envelope of the device curve from the retained counts (upward-biased; not a confidence interval) → `data/bootstrap_device_curve.json`. |
| `retention_matched_control.py` | Noiseless re-run at the device's kept-shot count → `data/2026-09-05_retention_matched_control.json`. |
| `analytic_coverage_leak.py` | Superseded analytic shot-noise model (`fp_rate = 1.0`) → `data/_superseded/analytic_coverage_leak.json`. |
| `nk_falsifier.py` | `n_k` momentum-distribution forecast (committed output `data/2026-09-01_nk_falsifier.json`). Manifest-hashed; writes to `./06_results/`. |
| `nk_stage0_gate.py`, `nk_stage0_gate_v2.py`, `nk_stage0_gate_v3.py` | Stage-0 gate builds v1–v3 → `data/2026-09-02_nk_stage0_gate{,_v2,_v3}.json`. Manifest-hashed; write to `./06_results/`. |
| `nk_stage0_v3R2.py` | The sealed v3-R2 resolution → `data/2026-09-02_nk_stage0_v3R2.json`. Manifest-hashed; writes to `./06_results/`. |
| `nk_stage1_seal.py` | Seals `data/manifest_nk_device_v1.json` (+ `.sha256`); refuses to overwrite it. |
| `nk_stage2_device_job.py` | The manifest-v1 on-device `n_k` job; never run. |
| `nk_stage0_gate_v4.py` | The amendment-validation gate rows → `data/2026-09-04_nk_v4_row_*.json` (never overwrites a committed row). |
| `nk_stage0_gate_v4_combine.py` | Applies the sealed amendment rule → `data/2026-09-05_nk_v4_verdict.json` (+ `.sha256`): AMENDMENT-DENIED, gate-negative. |
| `nk_stage2_device_job_v2_amendment.py` | The denied amendment's job (`day_of_rule_v2`); manifest v2 never sealed; never run. |
| `paper_style.py` | Plotting identity used by the superseded raster scripts. |
| `cache/{akw,sqw,current,teeth,teeth_shared}_L12.npz` | Committed `L=12` Lanczos caches (≈0.22 MB) read by `export_akw/sqw/teeth_dat.py`; `make cache` rebuilds them if missing. |
| `_superseded/` | Quarantined prototypes (`run_moment_bound_theorem.py`, `run_heron_screen.py`, `make_akw_figure.py`, `make_sqw_figure.py`, `make_paper_figs.py`, `export_falsifier_dat.py`); stale or circular, see its README. |

## `data/` — committed results and the sealed records

- **Sealed records (read-only):** `prereg.json` + `prereg.sha256` and
  `blind_{instances_public,labels_sealed,verdicts}.json` (the blinded pre-registration);
  `manifest_nk_device_v1.json` + `.sha256` (the `n_k` device manifest, 2026-09-02);
  `2026-09-05_nk_v4_verdict.json` + `.sha256` (the amendment gate verdict: AMENDMENT-DENIED) with its
  row files `2026-09-04_nk_v4_row_{R1ro0.018,R2,R3}.json`. `.gitattributes` keeps every file whose
  SHA-256 is recorded in a seal byte-exact on every platform (see `docs/REPRODUCE.md`).
- `2026-08-24_blind_harness_score.json` (the sealed primary endpoint) and `separating_demonstration.json`.
- `heron_counts_matched_L8_*.json` — the real `ibm_fez` `L=8` retained counts (four shot budgets) and the
  Aer dry run; `m0_hat`/`m1_hat` in them are exact classical values. `heron_counts_*DRYRUN*.json` are dry runs.
- `heron_spectral.json` — the companion study's reconstructed spectra (source of `fig_heron`).
- `2026-09-05_retention_matched_control.json`, `aer_noiseless_baseline.json`, `bootstrap_device_curve.json` — the device-curve controls.
- Every other `*.json` is the committed output of the `src/` script named in its `_provenance` (see the table above).
- `_superseded/` — the superseded analytic coverage model, precursor-program data not used by the manuscript,
  and byte-identical duplicates of `paper/figs/*.dat` (see its README).

## `paper/` — 2026-09-02 draft under revision

`README.md` (status and known issues), `main.tex` (master), `figstyle.tex` (shared figure identity),
`main.bbl` (frozen bibliography), `refs.bib`, `main.pdf` (compiled draft, 30 pp), `figs/` (native `.tex`
fragments + `.dat` + raster `.png`), `figs/src/fig_circuit_qtk.tex` (source of `fig_circuit_qtk.pdf`),
`figs/_superseded/` (files no figure reads; see its README).

## `notebooks/`
`00_Reproduce_Everything.ipynb` — narrated run of the committed generators (optional cells need `qiskit-aer`).

## `docs/`
`REPRODUCE.md` (figure/number → script → command, sealed-record checks, known gaps), `FILE_INDEX.md`
(this file), `img/_superseded/` (the retired 2026-08-22 README thumbnails; see its README).

## `_superseded/`
The 2026-09-02 arXiv upload bundle and its upload guide; not submitted; kept as a record (see its README).

## Root
`README.md`, `CITATION.cff`, `LICENSE` (MIT code + CC-BY-4.0 paper), `requirements.txt`, `Makefile`,
`.gitattributes` (byte-exact files whose SHA-256 is recorded in a seal), `.gitignore`.
