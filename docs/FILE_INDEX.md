# File index

Every file in this repository, described. Updated 2026-09-27 from the working tree: the phase-A scripts
and their dated outputs are included (marked **2026-09-27**). Build by-products that `.gitignore` excludes
(`*.aux`, `*.log`, `*.out`, `*.blg`, `__pycache__/`) are not listed.

## `src/` — engine, moment tests, figure generators, blinded protocol, device jobs

Scripts write their JSON into `data/` and plotted data into `paper/figs/` (paths relative to `src/`).

| File | Role |
|---|---|
| `verify.py` | Fast smoke test: the current-probe moments, the shot-budget thresholds τ₀–τ₂, and the joint battery catching (via `m₂`) an auto-calibrated truncation (`d=97`) the lone first moment passes; deterministic `lexsort` order since **2026-09-27** (checks `|Δm₁| = 0.186`, `|Δm₂| = 2.729`). |
| `hubbard_ed.py` | Vendored exact-diagonalization utilities (Jordan–Wigner Hubbard chain; numpy/scipy only), imported by the exporters and demonstrations. |
| `spectral_lanczos.py` | Exact sector Lanczos + Haydock engine: ground state, `A(k,ω)`, `S(q,ω)`, `S^zz`, current probe, truncation sweeps; writes `src/cache/`. Hashed in the `n_k` manifest and kept byte-for-byte. |
| `interval_moment.py` | Shot-budget interval-moment closure at the `ibm_fez` budget, deterministic truncation order → `data/2026-09-27_interval_moment_closure.json` (**2026-09-27**; the argsort-era `data/2026-08-22_interval_moment_closure.json` is kept as the record). |
| `interval_moment_mc.py` | Monte-Carlo interval-moment test (`fig_momentmc`), power at every truncation depth → `data/2026-09-27_interval_moment_mc.json`, `momentmc_{hist,power,power_scan,scalars}.dat` (**2026-09-27**; `data/2026-08-24_interval_moment_mc.json` kept as the record). |
| `interval_battery.py` | Joint `(m₀,m₁,m₂)`+Hankel–Stieltjes battery: every-d scan at 0/2/4% bias, shell-boundary cuts, tie-break sensitivity → `data/2026-09-27_interval_battery.json` (**2026-09-27**; `data/interval_battery_showcase.json` kept as the record). |
| `within_sector_control.py` | Within-sector detectability: moment-preserving redistributions and Vandermonde conditioning → `data/2026-08-24_within_sector_control.json`. |
| `within_sector_lp.py` | **2026-09-27.** LP ranges of Δ_{K+1}, null-space dimensions, conditioning on raw and rescaled nodes, three explicit counterexamples, the shifted-power bound → keys `R3_within_sector_lp`, `M1_shifted_power_bound` of `data/2026-09-27_theory_numerics.json`. |
| `christoffel_tolerance_lp.py` | **2026-09-27.** Tolerance-inflated extremal atom (LP and closed form), U/t=4 and U/t=8 measures → key `R4_tolerance_inflated_atom`. |
| `estimator_form_offeigenstate.py` | **2026-09-27.** First-moment estimator forms on truncated (non-eigen) states → key `R6_estimator_forms_off_eigenstate`. |
| `small_checks.py` | **2026-09-27.** Single writer of `data/2026-09-27_theory_numerics.json`; `--all` runs the six theory scripts, then writes `R11_trace_over_dim`, `m4_psd_and_symmetry`, `M5_christoffel_claims`, `M6_inverse_moment_claim` (≈15 min). |
| `markov_krein_window.py` | Two-sided Markov–Krein window-mass bound → `data/2026-08-28_markov_krein_window.json`. |
| `necessary_sufficient_composition.py` | Rewritten **2026-09-27**: Hankel lower and Hausdorff upper bounds on `m₂` (two necessary conditions) for four support intervals → key `R5_hausdorff_weighted_support`. The earlier version and its output `data/necessary_sufficient_composition.json` are kept as the record. |
| `chigap_check.py` | Christoffel empty-region check on the `L=6` current measure → `data/2026-08-24_chigap_check.json`. |
| `run_sumrule_falsifier.py` | Current-response demonstration: the total-weight `m₀` sum rule is blind, `m₁` detects → `data/2026-08-18_sumrule_falsifier.json`; defines `current_operator`. |
| `run_collective_sumrule_falsifier.py` | The screen on the collective (charge/spin) channels → `data/2026-08-26_collective_sumrule_falsifier.json`, `collective_screen_*.dat`. |
| `run_inverse_moment_falsifier.py` | The negative-order `m₋₁` (optical f-sum) test → `data/2026-08-28_inverse_moment_falsifier.json` (a rerun reproduces it and the figure files bit for bit; checked 2026-09-27), `inverse_falsifier_*.dat`, `inverse_falsifier_table.tex`. Relative shifts only; the noise-normalized comparison is `--r8`. With `--r8` (**2026-09-27**) it writes only key `R8_inverse_moment_sensitivity` and no figure file. |
| `run_falsifier_teeth.py` | The teeth at `L=6` (independent estimator vs a circular control) → `data/2026-08-18_falsifier_teeth.json`. |
| `run_teeth_shared.py` | The shared-state teeth → `data/2026-08-26_teeth_shared_state.json`, `teeth_shared.dat`. |
| `run_gausslaw_falsifier.py` | U(1) quantum link model: Hamiltonian and Gauss-law operators (`build`) and the bit-flip sweep → `data/2026-08-18_gausslaw_falsifier.json` (not what `fig_gausslaw` plots). |
| `export_gauss_state.py` | **2026-09-27.** Generator of `paper/figs/gauss_state.dat` (`fig_gausslaw`) with an explicit admixture rule → also `data/2026-09-27_gauss_state.json`. |
| `export_akw_dat.py` | `fig_akw` fragments and rasters from `src/cache/akw_L12.npz`. |
| `export_sqw_dat.py` | `fig_sqw` fragments and rasters from `src/cache/sqw_L12.npz`. |
| `check_sqw_extent.py` | **2026-09-27.** What the `fig_sqw` rasters draw below the charge-onset line and what the display windows clip → `data/2026-09-27_sqw_extent_check.json`. |
| `export_teeth_dat.py` | `fig_teeth` fragment from `src/cache/teeth_L12.npz`. |
| `export_bracketing_dat.py` | `fig_bracketing` fragments (Gauss–Radau bracketing); also writes `bracket_maxW.dat`, which no figure reads. |
| `export_christoffel_dat.py` | `fig_christoffel` fragments; since **2026-09-27** exact `max_t Wₙ` by root-finding and weighted atoms only; also writes keys `R9_christoffel_max_Wn`, `m1_radau_outside_hull`, `m2_rescaled_frame`. |
| `export_momentcone_dat.py` | `fig_momentcone` point. |
| `blind_preregister.py` | Seals the blinded protocol (`data/prereg.json` + `.sha256`); refuses to overwrite the committed seal. |
| `blind_generate.py` | Generates the 300 blinded instances from the sealed prereg (refuses to overwrite the sealed outputs). |
| `blind_classify.py` | Applies the frozen battery blind → `data/blind_verdicts.json` (refuses to overwrite it). |
| `blind_score.py` | Unblinds and scores → `data/2026-08-24_blind_harness_score.json` (TPR/FPR with Wilson intervals, per class). |
| `blind_addenda.py` | **2026-09-27.** Post hoc addenda from the sealed record (frozen rule run verbatim, endpoint, FPR sensitivity, exact rebuild and post hoc shared-state rescoring of the truncations, Krylov residuals, calibrated scope) → `data/2026-09-27_blind_addenda.json`. The rebuild uses the stored tie order, re-validated against the sealed record on every run; it aborts and writes nothing otherwise. |
| `separating_counts.py` | Recounts `fig_separating` from the sealed record and, since **2026-09-27**, exports `paper/figs/sep_*.dat` and `data/2026-09-27_separating_counts.json`. |
| `bond_moment_estimator.py` | Bond-basis (Pauli-group) moment estimator, simulated → `data/2026-08-25_bond_moment_estimator_L{4,6}.json`. |
| `estimator_scaling.py` | Circuit-count scaling of that estimator → `data/2026-08-25_estimator_circuit_scaling.json`. |
| `verify_mm1_measurement_cost.py` | Checks that `m₋₁` needs no measurement settings beyond the `m₀`/`m₁` plan (prints only). |
| `offdiag_noise_forecast.py` | Off-diagonal depolarizing-bias floor under `ibm_fez`-anchored noise → `data/2026-08-30_offdiag_noise_forecast.json`. |
| `offdiag_gsurface.py` | That bias floor vs circuit depth → `data/2026-08-31_offdiag_gsurface.json`. |
| `joint_covariance_composition.py` | Same-shot covariance of the independent estimate and the reconstruction (simulation) → `data/2026-09-01_joint_covariance.json`. |
| `lever1_flip_experiment.py` | Pre-registered decision experiment on that covariance (simulation; no power gain) → `data/2026-09-01_lever1_{prereg_seal,oracle}.json` (re-verifies the committed seal instead of overwriting it). |
| `beyond_ed_dmrg_catch.py` | `L=24` DMRG-oracle check of the coverage screen on a sampled subspace (needs `physics-tenpy`) → `data/2026-08-28_beyond_ed_dmrg_catch.json`. |
| `hardware_matched_job_L8.py` | The `ibm_fez` `L=8` coverage job (determinant + R_xx·R_yy circuits; correct bit order; post-selection only): Aer dry run (default) or QPU (`RUN=1`; measurement twirling, Pauli gate twirling and dynamical decoupling, no readout-error mitigation) → `data/heron_counts_matched_L8_*.json`. |
| `delta0_reference_mc.py` | **2026-09-27.** Reference distributions for the device coverage residual (noiseless raw- and retained-matched replicas, uniform-noise mixture, uniform in-sector sampler; percentile ranks) → `data/2026-09-27_delta0_reference_mc.json`. Imports `hardware_matched_job_L8.py`; needs qiskit; no QPU. |
| `export_device_dat.py` | **2026-09-27.** `fig_device` data `device_{points,band,uniform}.dat` from `data/2026-09-27_delta0_reference_mc.json` (formatting only). |
| `export_heron_dat.py` | **2026-09-27.** `fig_heron` data `heron_{exact,hw}.dat` from `data/heron_spectral.json` (formatting only; byte-identical to the committed files). |
| `hardware_matched_job.py` | The `L=6` matched job on the companion's circuits; only its Aer dry run was executed → `data/heron_counts_matched_DRYRUN.json`. |
| `hardware_blind_job.py` | A blinded hardware job builder for the bond-basis estimator; never run on hardware; its state preparation is a placeholder, not a ground-state preparation; only a dry-run output is committed (`data/heron_counts_DRYRUN.json`, from an earlier version of the script). |
| `aer_noiseless_baseline.py` | One noiseless Aer seed of the same `L=8` circuits → `data/aer_noiseless_baseline.json` (superseded as a comparator by `delta0_reference_mc.py`). |
| `bootstrap_device_curve.py` | Bootstrap envelope of the device curve from the retained counts (upward-biased; not a test) → `data/bootstrap_device_curve.json`. |
| `retention_matched_control.py` | One noiseless draw per budget at the device's kept-shot count, split equally over circuits → `data/2026-09-05_retention_matched_control.json` (see `data/README.md` for its percentile ranks). |
| `analytic_coverage_leak.py` | Superseded analytic shot-noise model (`fp_rate = 1.0`) → `data/_superseded/analytic_coverage_leak.json`. |
| `nk_falsifier.py` | `n_k` momentum-distribution forecast (committed output `data/2026-09-01_nk_falsifier.json`). Manifest-hashed; writes to `./06_results/`. |
| `nk_stage0_gate.py`, `nk_stage0_gate_v2.py`, `nk_stage0_gate_v3.py` | Stage-0 gate builds v1–v3 → `data/2026-09-02_nk_stage0_gate{,_v2,_v3}.json`. Manifest-hashed; write to `./06_results/`. |
| `nk_stage0_v3R2.py` | The sealed v3-R2 resolution → `data/2026-09-02_nk_stage0_v3R2.json`. Manifest-hashed; writes to `./06_results/`. |
| `nk_stage1_seal.py` | Seals `data/manifest_nk_device_v1.json` (+ `.sha256`); refuses to overwrite it. |
| `nk_stage2_device_job.py` | The manifest-v1 on-device `n_k` job; never run. |
| `nk_stage0_gate_v4.py` | The amendment-validation gate rows → `data/2026-09-04_nk_v4_row_*.json` (never overwrites a committed row). |
| `nk_stage0_gate_v4_combine.py` | Applies the sealed amendment rule → `data/2026-09-05_nk_v4_verdict.json` (+ `.sha256`): AMENDMENT-DENIED, gate-negative. Emits the `finding` as corrected by hand on 2026-09-05 (see `docs/REPRODUCE.md`, Sealed records), so a recombination in a scratch copy reproduces the sealed file byte for byte. |
| `nk_stage2_device_job_v2_amendment.py` | The denied amendment's job (`day_of_rule_v2`); manifest v2 never sealed; never run. |
| `paper_style.py` | Plotting identity used by the superseded raster scripts. |
| `cache/{akw,sqw,current,teeth,teeth_shared}_L12.npz` | Committed `L=12` Lanczos caches (≈0.22 MB) read by `export_akw/sqw/teeth_dat.py`; `make cache` rebuilds them if missing. |
| `_superseded/` | Quarantined prototypes (`run_moment_bound_theorem.py`, `run_heron_screen.py`, `make_akw_figure.py`, `make_sqw_figure.py`, `make_paper_figs.py`, `export_falsifier_dat.py`) and, since 2026-09-27, `necessary_sufficient_composition_2026-08.py`; stale, circular or replaced; see its README. |

## `data/` — committed results and the sealed records

`data/README.md` lists which dated files supersede which. Nothing is deleted.

- **Sealed records (read-only):** `prereg.json` + `prereg.sha256` and
  `blind_{instances_public,labels_sealed,verdicts}.json` (the blinded pre-registration);
  `manifest_nk_device_v1.json` + `.sha256` (the `n_k` device manifest, 2026-09-02);
  `2026-09-05_nk_v4_verdict.json` + `.sha256` (the amendment gate verdict: AMENDMENT-DENIED) with its
  row files `2026-09-04_nk_v4_row_{R1ro0.018,R2,R3}.json`. `.gitattributes` keeps every file whose
  SHA-256 is recorded in a seal byte-exact on every platform (see `docs/REPRODUCE.md`).
- `2026-08-24_blind_harness_score.json` (the sealed primary endpoint) and `separating_demonstration.json`
  (the 2026-08-31 spurious-only count, 24 / 15 / 9; corrected by `2026-09-27_separating_counts.json`).
- `heron_counts_matched_L8_*.json` — the real `ibm_fez` `L=8` retained counts (four shot budgets) and the
  Aer dry run; `m0_hat`/`m1_hat` in them are exact classical values. `heron_counts_*DRYRUN*.json` are dry runs.
- `heron_spectral.json` — the companion study's reconstructed spectra (source of `fig_heron`); that `L=6` run
  was post-selected in reversed qubit order, so all its retained determinants came from device errors (see
  `docs/REPRODUCE.md`, Hardware).
- `2026-09-05_retention_matched_control.json`, `aer_noiseless_baseline.json`, `bootstrap_device_curve.json`
  — single-draw device-curve controls, kept as records; the reference distributions are in
  `2026-09-27_delta0_reference_mc.json`.
- **2026-09-27 (phase A):** `2026-09-27_delta0_reference_mc.json`, `2026-09-27_blind_addenda.json`,
  `2026-09-27_separating_counts.json`, `2026-09-27_interval_moment_closure.json`,
  `2026-09-27_interval_battery.json`, `2026-09-27_interval_moment_mc.json`, `2026-09-27_gauss_state.json`,
  `2026-09-27_sqw_extent_check.json`, `2026-09-27_theory_numerics.json` (one key per plan item, each with its
  own provenance). Generators in
  `docs/REPRODUCE.md`.
- Every other `*.json` is the committed output of the `src/` script named in its `_provenance` (see the table above).
- `_superseded/` — the superseded analytic coverage model, precursor-program data not used by the manuscript,
  copies of `.dat` files from 2026-08-22, and (`pre_wording_pass_2026-09-27/`) the earlier versions of seven
  records whose strings were rescoped on 2026-09-27 with every number unchanged (see its README and
  `data/README.md`).

## `paper/` — the manuscript, revised 2026-09-27

`README.md` (status, what is still open, and the defects of the 2026-09-05 text that the revision corrected), `main.tex` (master),
`figstyle.tex` (shared figure identity), `main.bbl` (frozen bibliography), `refs.bib`, `main.pdf` (38 pp,
built from the revised text), `figs/` (native `.tex` fragments + `.dat` + raster `.png`),
`figs/src/fig_circuit_qtk.tex` (source of `fig_circuit_qtk.pdf`), `figs/_superseded/` (files no figure
reads, and the versions replaced on 2026-09-27; see its README).

New or regenerated in `figs/` on 2026-09-27: `device_{points,band,uniform}.dat` (`fig_device`, now two
panels), `sep_{ac_pres_region,ac_chg_region,region_nonac}.dat` and the regenerated `sep_{clean,krylov,trunc,ac_other}.dat`
(`fig_separating`; the unused `sep_ac_sep.dat` moved to `figs/_superseded/sep_2026-08-31/`), `momentmc_power_scan.dat` and the regenerated `momentmc_power.dat` (`fig_momentmc`),
`christoffel_maxW.dat` and `christoffel_poles.dat` (exact maxima; weighted atoms), `gauss_state.dat`
(regenerated byte-identical); changed sources `fig_device.tex`, `fig_separating.tex`, `fig_momentmc.tex`,
`fig_sqw.tex`, `fig_christoffel.tex`, `fig_gausslaw.tex`.

## `notebooks/`
`00_Reproduce_Everything.ipynb` — narrated run of the committed generators, the 2026-09-27 (phase A) scripts
included, each in its own process from `src/`; committed without outputs (`make reproduce` writes the
executed copy to `out/`). The device reference Monte Carlo is behind a switch (`RUN_DEVICE_MC`, off by
default); the optional cells need `qiskit`/`qiskit-aer`.

## `docs/`
`REPRODUCE.md` (figure/number → script → command, phase-A outputs, sealed-record checks, known gaps),
`FILE_INDEX.md` (this file), `img/_superseded/` (the retired 2026-08-22 README thumbnails; see its README).

## `_superseded/`
The 2026-09-02 arXiv upload bundle and its upload guide; not submitted; kept as a record (see its README).

## Root
`README.md`, `CITATION.cff`, `LICENSE` (MIT code + CC-BY-4.0 paper), `requirements.txt`, `Makefile`,
`.gitattributes` (byte-exact files whose SHA-256 is recorded in a seal), `.gitignore`.
