# File index

Every file in this repository, described (state of 2026-10-02). Build by-products that `.gitignore` excludes
(`*.aux`, `*.log`, `*.out`, `*.blg`, `*.toc`, `__pycache__/`, `06_results/`, `out/`) are not listed. Dates in bold mark
the day a script or record was added or changed.

## `src/` — engine, moment tests, figure generators, blinded protocol, device jobs

Scripts write their JSON into `data/` and plotted data into `paper/figs/` (paths relative to `src/`).

| File | Role |
|---|---|
| `verify.py` | Fast smoke test: the current-probe moments, the shot-budget thresholds τ₀–τ₂, and the joint battery catching (via `m₂`) an auto-calibrated truncation (`d=97`) the lone first moment passes (deterministic `lexsort` order; checks `|Δm₁| = 0.186`, `|Δm₂| = 2.729`); since **2026-10-02** also a record check of every number `paper/main.tex` prints from five records. |
| `hubbard_ed.py` | Vendored exact-diagonalization utilities (Jordan–Wigner Hubbard chain; numpy/scipy only), imported by the exporters and demonstrations. |
| `spectral_lanczos.py` | Exact sector Lanczos + Haydock engine: ground state, `A(k,ω)`, `S(q,ω)`, `S^zz`, current probe, truncation sweeps; writes `src/cache/`. Hashed in the `n_k` manifest and kept byte-for-byte. |
| `interval_moment.py` | Shot-budget interval-moment closure at the `ibm_fez` budget, deterministic truncation order → `data/2026-09-27_interval_moment_closure.json` (**2026-09-27**; the argsort-era `data/2026-08-22_interval_moment_closure.json` is kept as the record). |
| `interval_moment_mc.py` | Monte-Carlo interval-moment test (`fig_momentmc`), power at every truncation depth → `data/2026-09-27_interval_moment_mc.json`, `momentmc_{hist,power,power_scan,scalars}.dat` (**2026-09-27**; `data/2026-08-24_interval_moment_mc.json` kept as the record). |
| `interval_battery.py` | Joint `(m₀,m₁,m₂)`+Hankel–Stieltjes battery: every-d scan at 0/2/4% bias, shell-boundary cuts, tie-break sensitivity → `data/2026-09-27_interval_battery.json` (**2026-09-27**; `data/interval_battery_showcase.json` kept as the record). |
| `within_sector_control.py` | Within-sector detectability: moment-preserving redistributions and Vandermonde conditioning → `data/2026-08-24_within_sector_control.json`. |
| `within_sector_lp.py` | **2026-09-27**, extended **2026-10-02.** LP ranges of Δ_{K+1} and the next-order shift at the optimum, null-space dimensions, conditioning on raw and rescaled nodes, three explicit counterexamples, the shifted-power bound with the sharp Chebyshev constant → keys `R3_within_sector_lp`, `M1_shifted_power_bound` of `data/2026-09-27_theory_numerics.json`. |
| `christoffel_tolerance_lp.py` | **2026-09-27.** Tolerance-inflated extremal atom (LP and closed form), U/t=4 and U/t=8 measures → key `R4_tolerance_inflated_atom`. |
| `estimator_form_offeigenstate.py` | **2026-09-27.** First-moment estimator forms on truncated (non-eigen) states → key `R6_estimator_forms_off_eigenstate`. |
| `small_checks.py` | **2026-09-27**, extended **2026-10-02.** Single writer of `data/2026-09-27_theory_numerics.json`; `--all` runs the seven theory scripts (incl. `full_support_rr.py`), then writes `R11_trace_over_dim`, `m4_psd_and_symmetry`, `M5_christoffel_claims`, `M6_inverse_moment_claim`, `MA6_estimator_shape_by_order`, `B2_one_node_m2_deficit`, `SOTA1_02_W1_transport` (≈15 min). |
| `markov_krein_window.py` | Two-sided Markov–Krein window-mass bound → `data/2026-08-28_markov_krein_window.json`. |
| `necessary_sufficient_composition.py` | Rewritten **2026-09-27**: Hankel lower and Hausdorff upper bounds on `m₂` (two necessary conditions) for four support intervals → key `R5_hausdorff_weighted_support`. The earlier version and its output `data/necessary_sufficient_composition.json` are kept as the record. |
| `chigap_check.py` | Christoffel empty-region check on the `L=6` current measure → `data/2026-08-24_chigap_check.json`. |
| `run_sumrule_falsifier.py` | Current-response demonstration: the total-weight `m₀` sum rule is blind, `m₁` detects → `data/2026-08-18_sumrule_falsifier.json`; defines `current_operator`. |
| `run_collective_sumrule_falsifier.py` | The screen on the collective (charge/spin) channels → `data/2026-08-26_collective_sumrule_falsifier.json`, `collective_screen_*.dat`. |
| `run_inverse_moment_falsifier.py` | Misplaced weight at fixed total and the negative-order `m₋₁` (optical f-sum) test → `data/2026-08-28_inverse_moment_falsifier.json` (a rerun reproduces every number and the figure files bit for bit; claim string corrected **2026-10-02**), `inverse_falsifier_*.dat`, `inverse_falsifier_table.tex`. Relative shifts only; the noise-normalized comparison is `--r8`. With `--r8` (**2026-09-27**) it writes only key `R8_inverse_moment_sensitivity` and no figure file. |
| `run_falsifier_teeth.py` | The truncation sweep at `L=6` (independent estimator vs a circular control) → `data/2026-08-18_falsifier_teeth.json`. |
| `run_teeth_shared.py` | The shared-state truncation sweep at `L=12` (`fig_teeth_shared`) → `data/2026-08-26_teeth_shared_state.json`, `teeth_shared.dat`. |
| `run_gausslaw_falsifier.py` | U(1) quantum link model: Hamiltonian and Gauss-law operators (`build`) and the bit-flip sweep → `data/2026-08-18_gausslaw_falsifier.json` (not what `fig_gausslaw` plots). |
| `export_gauss_state.py` | **2026-09-27.** Generator of `paper/figs/gauss_state.dat` (`fig_gausslaw`) with an explicit admixture rule → also `data/2026-09-27_gauss_state.json`. |
| `export_akw_dat.py` | `fig_akw` fragments and rasters from `src/cache/akw_L12.npz`. |
| `export_sqw_dat.py` | `fig_sqw` fragments and rasters from `src/cache/sqw_L12.npz`. |
| `check_sqw_extent.py` | **2026-09-27.** What the `fig_sqw` rasters draw below the charge-onset line and what the display windows clip → `data/2026-09-27_sqw_extent_check.json`. |
| `export_teeth_dat.py` | `teeth.dat` (idealized curve and circular control, plotted inside `fig_teeth_shared`) from `src/cache/teeth_L12.npz`. |
| `export_bracketing_dat.py` | `fig_bracketing` fragments (Gauss–Radau bracketing); since **2026-10-02** the band includes the exact `Wₙ` peaks; also writes `bracket_maxW.dat` (equal to `christoffel_maxW.dat`; checked by `verify.py`). |
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
| `offdiag_noise_forecast.py` | Off-diagonal depolarizing-bias floor under `ibm_fez`-anchored noise at the reference depth read from `data/sampler_cz_count.json` (295 CZ since **2026-10-02**) → `data/2026-08-30_offdiag_noise_forecast.json`. |
| `offdiag_gsurface.py` | That bias floor vs circuit depth (local, global, coherent), reference depths 295 and 213 CZ from `data/sampler_cz_count.json` (**2026-10-02**) → `data/2026-08-31_offdiag_gsurface.json`. |
| `sampler_cz_count.py` | **2026-10-02.** Transpiled CZ counts of the sampler circuits on `ibm_fez`: `FakeFez` transpilation of both job scripts' circuits and the ISA circuits actually submitted in the six device jobs (stored encoded, so the re-count runs offline) → `data/sampler_cz_count.json`. Needs qiskit, qiskit-ibm-runtime. |
| `joint_covariance_composition.py` | Same-shot covariance of the independent estimate and the reconstruction (simulation) → `data/2026-09-01_joint_covariance.json`. |
| `lever1_flip_experiment.py` | Pre-registered decision experiment on that covariance (simulation; no power gain) → `data/2026-09-01_lever1_{prereg_seal,oracle}.json` (re-verifies the committed seal instead of overwriting it). |
| `beyond_ed_dmrg_catch.py` | Code-scalability check of the shared-state test on sampled subspaces at `L = 8` (ED), 20 and 24 (DMRG references; needs `physics-tenpy`); since **2026-10-02** a degenerate subspace ground level is reported as a range → `data/2026-10-02_beyond_ed_dmrg_catch.json`. |
| `hardware_matched_job_L8.py` | The `ibm_fez` `L=8` coverage job (determinant + R_xx·R_yy circuits; correct bit order; post-selection only): Aer dry run (default) or QPU (`RUN=1`; measurement twirling, Pauli gate twirling and dynamical decoupling, no readout-error mitigation) → `data/heron_counts_matched_L8_*.json` (jobs created 2026-08-29 UTC; the records' `date` field is the script constant `DATE`). |
| `full_support_rr.py` | **2026-10-02.** Rayleigh–Ritz on the full ground-state support of the calibration ring (U/t = 3–8) under the default battery, the one-node Krylov deficit, and the deterministic leak scores of every truncation at U/t = 4 against the L1 error → `data/2026-10-02_full_support_rr.json`. |
| `gate_tp_analysis.py` | **2026-10-02.** Reanalysis of the two 2026-08-21 `L=6` `ibm_fez` jobs from their raw records (`data/heron_gatetp_2026-08-21/`): decode, the submitting notebook's pipeline in both bit orders, a noiseless bit-order check, what a support-based test sees when the support is complete → `data/2026-10-02_gatetp_reanalysis.json`. Needs qiskit. |
| `delta0_reference_mc.py` | **2026-09-27.** Reference distributions for the device coverage residual (noiseless raw- and retained-matched replicas, uniform-noise mixture, uniform in-sector sampler; percentile ranks) → `data/2026-09-27_delta0_reference_mc.json`. Imports `hardware_matched_job_L8.py`; needs qiskit; no QPU. |
| `export_device_dat.py` | **2026-09-27.** `fig_device` data `device_{points,band,uniform}.dat` from `data/2026-09-27_delta0_reference_mc.json` (formatting only). |
| `export_heron_dat.py` | **2026-09-27.** `fig_heron` data `heron_{exact,hw}.dat` from `data/heron_spectral.json` (formatting only; byte-identical to the committed files). |
| `hardware_matched_job.py` | The `L=6` matched job on the companion's circuits; only its Aer dry run was executed from it → `data/heron_counts_matched_DRYRUN.json`. Its `build_circuit` is the circuit family of the 2026-08-21 jobs (`sampler_cz_count.py`). |
| `hardware_blind_job.py` | A blinded hardware job builder for the bond-basis estimator; never run on hardware; its state preparation is a placeholder, not a ground-state preparation; only a dry-run output is committed (`data/heron_counts_DRYRUN.json`, from an earlier version of the script). |
| `aer_noiseless_baseline.py` | One noiseless Aer seed of the same `L=8` circuits → `data/aer_noiseless_baseline.json` (superseded as a comparator by `delta0_reference_mc.py`). |
| `bootstrap_device_curve.py` | Bootstrap envelope of the device curve from the retained counts (upward-biased; not a test) → `data/bootstrap_device_curve.json`. |
| `retention_matched_control.py` | One noiseless draw per budget at the device's kept-shot count, split equally over circuits → `data/2026-09-05_retention_matched_control.json` (see `data/README.md` for its percentile ranks). |
| `analytic_coverage_leak.py` | Superseded analytic shot-noise model (`fp_rate = 1.0`) → `data/_superseded/analytic_coverage_leak.json`. |
| `nk_falsifier.py` | `n_k` momentum-distribution forecast (committed output `data/2026-09-01_nk_falsifier.json`). Manifest-hashed; writes to `06_results/` at the repository root (create it first). |
| `nk_stage0_gate.py`, `nk_stage0_gate_v2.py`, `nk_stage0_gate_v3.py` | Stage-0 gate builds v1–v3 → `data/2026-09-02_nk_stage0_gate{,_v2,_v3}.json`. Manifest-hashed; write to `06_results/` at the repository root (create it first). |
| `nk_stage0_v3R2.py` | The sealed v3-R2 resolution → `data/2026-09-02_nk_stage0_v3R2.json`. Manifest-hashed; writes to `06_results/` at the repository root (create it first). |
| `nk_stage1_seal.py` | Seals `data/manifest_nk_device_v1.json` (+ `.sha256`); refuses to overwrite it. |
| `nk_stage2_device_job.py` | The manifest-v1 on-device `n_k` job; it submitted no job (the versions launched on 2026-09-03 are in `nk_record/`). |
| `nk_stage0_gate_v4.py` | The amendment-validation gate rows → `data/2026-09-04_nk_v4_row_*.json` (never overwrites a committed row). |
| `nk_stage0_gate_v4_combine.py` | Applies the sealed amendment rule → `data/2026-09-05_nk_v4_verdict.json` (+ `.sha256`): AMENDMENT-DENIED, gate-negative. Emits the `finding` as corrected by hand on 2026-09-05 (see `docs/REPRODUCE.md`, Sealed records), so a recombination in a scratch copy reproduces the sealed file byte for byte. |
| `nk_stage2_device_job_v2_amendment.py` | The denied amendment's job (`day_of_rule_v2`); manifest v2 never sealed; never run. |
| `nk_record/` | **2026-10-02.** The scripts of the `n_k` device attempt as run on 2026-09-03 to 2026-09-05 (two job launches, calibration watchers, the v1 chain checks), the logged outputs, and `nk_attempt_record.py`, the generator of `data/2026-09-05_nk_device_attempt_record.json` (`--offline` recomputes it from the stored calibration tables); see its README. |
| `paper_style.py` | Plotting identity used by the superseded raster scripts. |
| `cache/{akw,sqw,current,teeth,teeth_shared}_L12.npz` | Committed `L=12` Lanczos caches (≈0.22 MB) read by `export_akw/sqw/teeth_dat.py`; `make cache` rebuilds them if missing. |
| `_superseded/` | Quarantined prototypes (`run_moment_bound_theorem.py`, `run_heron_screen.py`, `make_akw_figure.py`, `make_sqw_figure.py`, `make_paper_figs.py`, `export_falsifier_dat.py`) and, since 2026-09-27, `necessary_sufficient_composition_2026-08.py`; stale, circular or replaced; see its README. |

## `data/` — committed results, sealed records, device records

`data/README.md` says which dated file supersedes which; nothing is deleted. Every record and its
generator (in `src/` unless stated):

| Record | Generator | Note |
|---|---|---|
| `prereg.json`, `prereg.sha256` | `blind_preregister.py` | **sealed** blinded pre-registration |
| `blind_instances_public.json`, `blind_labels_sealed.json` | `blind_generate.py` | **sealed** |
| `blind_verdicts.json` | `blind_classify.py` | **sealed** |
| `2026-08-24_blind_harness_score.json` | `blind_score.py` | primary endpoint as implemented |
| `separating_demonstration.json` | an earlier version of `separating_counts.py` (2026-08-31) | spurious-only count 24 / 15 / 9; corrected by `2026-09-27_separating_counts.json` |
| `2026-08-18_sumrule_falsifier.json` | `run_sumrule_falsifier.py` | |
| `2026-08-18_falsifier_teeth.json` | `run_falsifier_teeth.py` | |
| `2026-08-18_gausslaw_falsifier.json` | `run_gausslaw_falsifier.py` | bit-flip sweep (not what `fig_gausslaw` plots) |
| `2026-08-22_interval_moment_closure.json` | `interval_moment.py` (argsort era) | superseded by `2026-09-27_interval_moment_closure.json` |
| `interval_battery_showcase.json` | `interval_battery.py` (argsort era) | superseded by `2026-09-27_interval_battery.json` |
| `2026-08-24_interval_moment_mc.json` | `interval_moment_mc.py` (grid only) | superseded by `2026-09-27_interval_moment_mc.json` |
| `2026-08-24_chigap_check.json` | `chigap_check.py` | |
| `2026-08-24_within_sector_control.json` | `within_sector_control.py` | `weight_moved_frac` is ‖δ‖₁/m₀ |
| `2026-08-25_bond_moment_estimator_L{4,6}.json` | `bond_moment_estimator.py` (`BME_L=4`; default L=6) | noiseless simulation |
| `2026-08-25_estimator_circuit_scaling.json` | `estimator_scaling.py` | |
| `2026-08-26_collective_sumrule_falsifier.json` | `run_collective_sumrule_falsifier.py` | |
| `2026-08-26_teeth_shared_state.json` | `run_teeth_shared.py` | plotted record (see `docs/REPRODUCE.md`, Known gaps) |
| `2026-08-28_beyond_ed_dmrg_catch.json` | `beyond_ed_dmrg_catch.py` (earlier version) | superseded by `2026-10-02_beyond_ed_dmrg_catch.json` |
| `2026-08-28_inverse_moment_falsifier.json` | `run_inverse_moment_falsifier.py` | |
| `2026-08-28_markov_krein_window.json` | `markov_krein_window.py` | |
| `2026-08-30_offdiag_noise_forecast.json` | `offdiag_noise_forecast.py` | fixed name; regenerated 2026-10-02 at 295 CZ |
| `2026-08-31_offdiag_gsurface.json` | `offdiag_gsurface.py` | fixed name; regenerated 2026-10-02 |
| `2026-09-01_joint_covariance.json` | `joint_covariance_composition.py` | simulation |
| `2026-09-01_lever1_prereg_seal.json`, `2026-09-01_lever1_oracle.json` | `lever1_flip_experiment.py` | seal and oracle stage of the same-sample study |
| `2026-09-01_nk_falsifier.json` | `nk_falsifier.py` | manifest-hashed script |
| `2026-09-02_nk_stage0_gate{,_v2,_v3}.json` | `nk_stage0_gate{,_v2,_v3}.py` | `_v3` is manifest-hashed |
| `2026-09-02_nk_stage0_v3R2.json` | `nk_stage0_v3R2.py` | manifest-hashed |
| `manifest_nk_device_v1.json`, `manifest_nk_device_v1.sha256` | `nk_stage1_seal.py` | **sealed** `n_k` manifest |
| `2026-09-04_nk_v4_row_{R1ro0.018,R2,R3}.json` | `nk_stage0_gate_v4.py` | amendment gate rows |
| `2026-09-05_nk_v4_verdict.json`, `.sha256` | `nk_stage0_gate_v4_combine.py` | **sealed** amendment verdict (DENIED) |
| `2026-09-05_nk_device_attempt_record.json` | `nk_record/nk_attempt_record.py` | **2026-10-02**; 0 jobs, 0 shots |
| `2026-09-05_retention_matched_control.json` | `retention_matched_control.py` | single draw; ranked in `2026-09-27_delta0_reference_mc.json` |
| `aer_noiseless_baseline.json` | `aer_noiseless_baseline.py` | one Aer seed; superseded as a comparator |
| `bootstrap_device_curve.json` | `bootstrap_device_curve.py` | upward-biased; not a test |
| `heron_counts_matched_L8_<job id>.json` (4) | `hardware_matched_job_L8.py` (`RUN=1`) | **raw** `ibm_fez` counts, 7 circuits × shots per job; jobs created 2026-08-29 UTC |
| `heron_counts_matched_L8_DRYRUN.json` | `hardware_matched_job_L8.py` (default) | Aer dry run |
| `heron_counts_matched_DRYRUN.json` | `hardware_matched_job.py` | Aer dry run |
| `heron_counts_DRYRUN.json`, `heron_counts_DRYRUN_m0m1.json` | earlier versions of `hardware_blind_job.py` (`NMOM` = 3 and 2) | dry runs; never run on hardware |
| `heron_spectral.json` | the companion's pipeline (copied from its deposit) | spectra of the companion's `L=6` run; read by `export_heron_dat.py` and `verify.py` |
| `heron_gatetp_2026-08-21/` | IBM Quantum REST API (`/jobs/{id}`, `/results`, `/metrics`), copied byte for byte | **2026-10-02**; raw records of the two 2026-08-21 jobs; SHA-256 in its README |
| `sampler_cz_count.json` | `sampler_cz_count.py` | **2026-10-02**; includes the encoded ISA circuits of all six device jobs |
| `necessary_sufficient_composition.json` | `_superseded/necessary_sufficient_composition_2026-08.py` | superseded by key `R5_hausdorff_weighted_support` |
| `2026-09-27_delta0_reference_mc.json` | `delta0_reference_mc.py` | regenerated 2026-10-02 (values bit-identical; monotonicity flags corrected) |
| `2026-09-27_blind_addenda.json` | `blind_addenda.py` | parts (f), (g) added 2026-10-02 |
| `2026-09-27_separating_counts.json` | `separating_counts.py` | |
| `2026-09-27_interval_{moment_closure,battery,moment_mc}.json` | `interval_moment.py`, `interval_battery.py`, `interval_moment_mc.py` | |
| `2026-09-27_gauss_state.json` | `export_gauss_state.py` | |
| `2026-09-27_sqw_extent_check.json` | `check_sqw_extent.py` | |
| `2026-09-27_theory_numerics.json` | `small_checks.py --all` (each key names its own script) | keys added 2026-10-02 |
| `2026-10-02_full_support_rr.json` | `full_support_rr.py` | |
| `2026-10-02_beyond_ed_dmrg_catch.json` | `beyond_ed_dmrg_catch.py` | |
| `2026-10-02_gatetp_reanalysis.json` | `gate_tp_analysis.py` | |
| `_superseded/` | see its README | superseded analytic coverage model, precursor-program data, 2026-08-22 `.dat` copies, pre-wording-pass versions of seven records |

`.gitattributes` keeps every file whose SHA-256 is recorded in a seal, and the raw IBM job records,
byte-exact on every platform, and checks out the inputs whose LF-normalised hashes the generators record
with LF line endings (see `docs/REPRODUCE.md`).

## `paper/` — the manuscript

`README.md` (structure, files, build), `main.tex` (master), `refs.bib`, `main.bbl` (generated by BibTeX),
`figstyle.tex` (shared figure identity), `main.pdf`, `SciPost.cls` + `SciPost_bibstyle.bst` (SciPost's
submission class and bibliography style, under their own terms), `figs/` (native `.tex` fragments + `.dat` +
raster `.png`; `fig_heron.tex` and `fig_teeth.tex` are kept but not input by `main.tex`),
`figs/src/fig_circuit_qtk.tex` (source of `fig_circuit_qtk.pdf`), `figs/_superseded/` (files no figure reads,
and earlier versions of replaced figure files; see its README).

## `notebooks/`
`00_Reproduce_Everything.ipynb` — narrated run of the committed generators, each in its own process from
`src/`; committed without outputs (`make reproduce` writes the executed copy to `out/`). The device reference
Monte Carlo is behind a switch (`RUN_DEVICE_MC`, off by default); the optional cells need
`qiskit`/`qiskit-aer`/`physics-tenpy`.

## `docs/`
`REPRODUCE.md` (figure/number → script → command, dated recomputations, sealed-record checks, hardware,
known gaps), `FILE_INDEX.md` (this file), `PROCESS_RECORD.md` (design history of the sealed tests, the
estimator form, the tie order, and the manuscript's earlier versions), `img/_superseded/` (retired
prototype thumbnails; see its README).

## `_superseded/`
An arXiv upload bundle of the 2026-09-02 draft and its upload guide; not submitted; kept as a record (see
its README).

## Root
`README.md`, `CITATION.cff`, `.zenodo.json` (archive metadata), `LICENSE` (MIT code; CC BY 4.0 paper, data
and documentation; SciPost template files excluded), `requirements.txt`, `Makefile`, `.gitattributes`
(byte-exact sealed and raw records; LF for hashed inputs), `.gitignore`.
