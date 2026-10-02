# `data/` — committed results, sealed records, device records

A correction is a **new dated file** (or, for the two forecast records with fixed names, a regeneration
recorded in the file's `run_date`); the earlier file stays as the record, and nothing is deleted. This page
says which file supersedes which. The generator of every record is listed in
[`../docs/FILE_INDEX.md`](../docs/FILE_INDEX.md); commands and key numbers are in
[`../docs/REPRODUCE.md`](../docs/REPRODUCE.md). Licence: CC BY 4.0 (see `../LICENSE`).

## Sealed records (read-only; their writers refuse to overwrite them)

`prereg.json` + `prereg.sha256`, `blind_{instances_public,labels_sealed,verdicts}.json`,
`manifest_nk_device_v1.json` + `.sha256`, `2026-09-05_nk_v4_verdict.json` + `.sha256` and its rows
`2026-09-04_nk_v4_row_{R1ro0.018,R2,R3}.json`. How to check them: `docs/REPRODUCE.md`, Sealed records.

## Recomputations of 2026-09-27

| File | Generator (`src/`) |
|---|---|
| `2026-09-27_delta0_reference_mc.json` | `delta0_reference_mc.py` (needs qiskit); regenerated 2026-10-02 with every value bit-identical, the monotonicity flags computed over all rows sorted by p (`exact_E_monotone_decreasing_in_p` is `false` at every budget) and LF-normalised input hashes |
| `2026-09-27_blind_addenda.json` | `blind_addenda.py` (reads the sealed record only); parts (f) noise scale of the sealed text and variance inflation, and (g) shared-state identity on the faithful controls, added 2026-10-02 |
| `2026-09-27_separating_counts.json` | `separating_counts.py` (reads the sealed record only) |
| `2026-09-27_interval_moment_closure.json` | `interval_moment.py` |
| `2026-09-27_interval_battery.json` | `interval_battery.py` |
| `2026-09-27_interval_moment_mc.json` | `interval_moment_mc.py` |
| `2026-09-27_gauss_state.json` | `export_gauss_state.py` |
| `2026-09-27_sqw_extent_check.json` | `check_sqw_extent.py` (`fig_sqw`: below-onset display values, window clipping) |
| `2026-09-27_theory_numerics.json` | `small_checks.py --all` (one key per item; each key names its own script); keys `MA6_estimator_shape_by_order`, `B2_one_node_m2_deficit`, `SOTA1_02_W1_transport` and the next-order fields of `R3_within_sector_lp` added 2026-10-02 |

## Records of 2026-10-02

| File | Generator | Content |
|---|---|---|
| `2026-10-02_full_support_rr.json` | `src/full_support_rr.py` | full-support Rayleigh–Ritz on the calibration ring; leak scores of every truncation at U/t = 4 |
| `sampler_cz_count.json` | `src/sampler_cz_count.py` | CZ counts of the sampler circuits on `ibm_fez`, with the encoded ISA circuits of all six device jobs and their creation times |
| `2026-10-02_beyond_ed_dmrg_catch.json` | `src/beyond_ed_dmrg_catch.py` | `L = 8, 20, 24` shared-state check against ED/DMRG, with the degenerate 500-shot row at `L=24` reported as a range |
| `2026-10-02_gatetp_reanalysis.json` | `src/gate_tp_analysis.py` | reanalysis of the two 2026-08-21 `L=6` jobs |
| `heron_gatetp_2026-08-21/` | IBM Quantum REST API, copied byte for byte | raw job, result and metrics records of those jobs (see its README) |
| `2026-09-05_nk_device_attempt_record.json` | `src/nk_record/nk_attempt_record.py` | the sealed `n_k` device attempt of 2026-09-03 to 2026-09-05: 0 jobs, 0 shots, with the logged outputs and re-queried calibrations |

Regenerated in place on 2026-10-02 (fixed file names): `2026-08-30_offdiag_noise_forecast.json` and
`2026-08-31_offdiag_gsurface.json` now take their reference depth from `sampler_cz_count.json` (295 CZ, the
median over the deployed `L=8` circuits; the earlier 292 is one deployed count, kept as a grid point and
reproduced to 10⁻¹³). `2026-08-28_inverse_moment_falsifier.json` changed only in its claim and provenance
strings. `heron_spectral.json` was replaced by the companion's corrected deposit: the arrays are identical,
and the provenance now says that `A_hw` and `A_exact` agree to 7.9×10⁻¹⁵ (not byte for byte) and that the
support was retained by post-selection, not recovered.

## Earlier records kept in place, and what supersedes them

| Record (kept) | Superseded by | Why |
|---|---|---|
| `separating_demonstration.json` | `2026-09-27_separating_counts.json` | It counted spurious-feature instances only (24 / 15 / 9). The shaded region holds 26: 20 weight-preserving + 4 weight-changing spurious features, truncation id 162 and clean false positive id 184; m₁ fires on 16, m₂ alone on 10. |
| `2026-08-22_interval_moment_closure.json`, `interval_battery_showcase.json`, `2026-08-24_interval_moment_mc.json` | `2026-09-27_interval_moment_closure.json`, `2026-09-27_interval_battery.json`, `2026-09-27_interval_moment_mc.json` | Truncation order by `argsort` broke ties in degenerate \|ψ₀\|² shells by floating-point noise; their d=98 values are tie-break draws. The 2026-09-27 files use a deterministic order and measure the tie sensitivity. |
| `2026-08-28_beyond_ed_dmrg_catch.json` | `2026-10-02_beyond_ed_dmrg_catch.json` | At `L=24` and 500 shots the subspace ground level is two-fold degenerate (splitting ~2×10⁻¹⁶; gap 0.0722 to the next level), so the earlier starting value Δ₀ = 44.2 depends on the eigenvector returned; over the ground manifold Δ₀ spans 40.1–44.5. Every other number reproduces to about 10⁻¹³. |
| `necessary_sufficient_composition.json` | key `R5_hausdorff_weighted_support` of `2026-09-27_theory_numerics.json` | Its "sufficient" upper bound `m₂ ≤ 573.68` is the Hausdorff bound on the a-priori range [0.148, 52.21]t, a necessary condition; the weighted support gives 147.6, but that support is oracle knowledge. |
| `aer_noiseless_baseline.json` | `2026-09-27_delta0_reference_mc.json` | One Aer seed; it sits at the 19.5 / 17.7 / 81.4 / 98.3th percentiles of 2000 noiseless replicas at 50k / 30k / 16k / 4k. |
| `2026-09-05_retention_matched_control.json` | `2026-09-27_delta0_reference_mc.json` | One draw per budget, split equally over circuits; it sits at the 95.4 / 93.6 / 99.65 / 78.9th percentiles of the equal-split distribution. |
| `bootstrap_device_curve.json` | (none; not a test) | Resampling of the retained counts; upward-biased (its 50k 95% interval excludes the point estimate). |
| `_superseded/analytic_coverage_leak.json` | `2026-09-27_delta0_reference_mc.json` | Analytic shot-noise model; modelled false-positive rate 1.0; its σ is 34–195× smaller than the Monte Carlo spread at the device's raw shots. |

### Wording pass of 2026-09-27 (numbers unchanged)

The generators of the records below were rescoped on 2026-09-27 (no internal workflow identifiers, no
"device-measurable", "same-sample" or "confirmed" claims; correct script paths and fillings; relative shift
distinguished from detection power). Each record was regenerated in a scratch copy and replaced only
because every number reproduced bit for bit; only provenance, claim and verdict strings changed. The
versions from before the pass are kept in `_superseded/pre_wording_pass_2026-09-27/`.

| Record | Generator (`src/`) | What changed |
|---|---|---|
| `2026-08-18_sumrule_falsifier.json` | `run_sumrule_falsifier.py` | script path; filling (hole-doped, N = 4 on 6 sites, not half filling); claim and verdict; workflow id removed; `wording_note` added |
| `2026-08-18_falsifier_teeth.json` | `run_falsifier_teeth.py` | script path; claim and verdict (idealised: independent branch on the exact ground state); workflow id removed; `wording_note` added |
| `2026-08-18_gausslaw_falsifier.json` | `run_gausslaw_falsifier.py` | script path; claim and verdict (no same-sample or cross-domain claim; evaluated exactly); workflow id removed; `wording_note` added |
| `2026-08-25_bond_moment_estimator_L4.json`, `2026-08-25_bond_moment_estimator_L6.json` | `bond_moment_estimator.py` (`BME_L=4`; default L=6) | claim: a noiseless-simulation check, not run on hardware |
| `2026-08-25_estimator_circuit_scaling.json` | `estimator_scaling.py` | method: `<ρ_π> = 0` by the reflection symmetry of the open chain (any filling), not half filling |
| `2026-08-28_inverse_moment_falsifier.json` | `run_inverse_moment_falsifier.py` (no `--r8`) | script path; claim (relative shift is not detection power; see key `R8_inverse_moment_sensitivity`); `sensitivity_meaning` and `wording_note` added |

Not regenerated: `2026-08-26_teeth_shared_state.json` (with `paper/figs/teeth_shared.dat`) keeps its
2026-08-26 wording ("device-realizable"), because a rerun moves the three smallest-d sweep points
(`argsort` tie order; see `docs/REPRODUCE.md`, Known gaps); the script now writes rescoped strings.

Two committed files need a reading note:
- `2026-08-24_within_sector_control.json`: `cases[0].weight_moved_frac = 0.494` is ‖δ‖₁/m₀, i.e. 24.7% of
  the weight transferred. The LP ranges and explicit counterexamples are in key `R3_within_sector_lp`.
- `2026-08-18_gausslaw_falsifier.json` is the bit-flip sweep; `fig_gausslaw` plots the ε-admixture curve
  of `2026-09-27_gauss_state.json`.

## Device records

- `heron_counts_matched_L8_<job id>.json` — the **raw** `ibm_fez` counts of this work's four `L=8` jobs
  (`per_time_counts`: 7 circuits × the job's shots per circuit, every bitstring, before post-selection),
  plus `heron_counts_matched_L8_DRYRUN.json` (Aer). The jobs were created on 2026-08-29 UTC
  (`sampler_cz_count.json`, `deployed.jobs.<id>.created`); the records' `date` field (2026-08-28) is the
  script constant `DATE` of `hardware_matched_job_L8.py`, not the run date. `m0_hat`/`m1_hat` in them are
  exact classical values; only the support set is device-derived.
- `heron_gatetp_2026-08-21/` — the raw records of the two `L=6` jobs of 2026-08-21, which enter no result:
  their retained support is the whole sector under either bit order.
- `heron_spectral.json` — the companion's `L=6` reconstructed spectra (copied from its deposit); that run
  was post-selected in reversed bit order, so all its retained determinants came from device errors; its
  raw counts are not deposited.
- The IBM job records carry the submitting account's `user_id` (an IBMid account identifier) as returned by
  the API; no token, key or instance identifier.

## `_superseded/`

The analytic coverage model (read by `delta0_reference_mc.py` only as the comparator for its spread
ratio), precursor-program data, the 2026-08-22 `.dat` copies and the pre-wording-pass versions of seven
records (above). No number of the manuscript rests on these files. See
[`_superseded/README.md`](_superseded/README.md).
