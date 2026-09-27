# `data/` — committed results, sealed records, device counts

Nothing in this folder is deleted or overwritten when a result is corrected: a correction is a **new dated
file**, and the earlier file stays as the record. This page says which is which. Generators and commands are
in [`../docs/REPRODUCE.md`](../docs/REPRODUCE.md); every file is described in
[`../docs/FILE_INDEX.md`](../docs/FILE_INDEX.md).

## Sealed records (read-only; their writers refuse to overwrite them)

`prereg.json` + `prereg.sha256`, `blind_{instances_public,labels_sealed,verdicts}.json`,
`manifest_nk_device_v1.json` + `.sha256`, `2026-09-05_nk_v4_verdict.json` + `.sha256` and its rows
`2026-09-04_nk_v4_row_{R1ro0.018,R2,R3}.json`. How to check them: `docs/REPRODUCE.md`, Sealed records.

## Revision phase A (2026-09-27)

| File | Generator (`src/`) |
|---|---|
| `2026-09-27_delta0_reference_mc.json` | `delta0_reference_mc.py` (needs qiskit) |
| `2026-09-27_blind_addenda.json` | `blind_addenda.py` (reads the sealed record only) |
| `2026-09-27_separating_counts.json` | `separating_counts.py` (reads the sealed record only) |
| `2026-09-27_interval_moment_closure.json` | `interval_moment.py` |
| `2026-09-27_interval_battery.json` | `interval_battery.py` |
| `2026-09-27_interval_moment_mc.json` | `interval_moment_mc.py` |
| `2026-09-27_gauss_state.json` | `export_gauss_state.py` |
| `2026-09-27_theory_numerics.json` | `small_checks.py --all` (one key per item; each key names its own script) |

## Earlier records kept in place, and what supersedes them

| Record (kept) | Superseded by | Why |
|---|---|---|
| `separating_demonstration.json` | `2026-09-27_separating_counts.json` | It counted spurious-feature instances only (24 / 15 / 9). The shaded region holds 26: 20 weight-preserving + 4 weight-changing spurious features, truncation id 162 and clean false positive id 184; m₁ fires on 16, m₂ alone on 10. |
| `2026-08-22_interval_moment_closure.json`, `interval_battery_showcase.json`, `2026-08-24_interval_moment_mc.json` | `2026-09-27_interval_moment_closure.json`, `2026-09-27_interval_battery.json`, `2026-09-27_interval_moment_mc.json` | Truncation order by `argsort` broke ties in degenerate \|ψ₀\|² shells by floating-point noise; their d=98 values are tie-break draws. The 2026-09-27 files use a deterministic order and measure the tie sensitivity. |
| `necessary_sufficient_composition.json` | key `R5_hausdorff_weighted_support` of `2026-09-27_theory_numerics.json` | Its "sufficient" upper bound `m₂ ≤ 573.68` is the Hausdorff bound on the a-priori range [0.148, 52.21]t, a necessary condition; the weighted support gives 147.6, but that support is oracle knowledge (it needs the exact spectrum). |
| `aer_noiseless_baseline.json` | `2026-09-27_delta0_reference_mc.json` | One Aer seed; it sits at the 19.5 / 17.7 / 81.4 / 98.3th percentiles of 2000 noiseless replicas at 50k / 30k / 16k / 4k. |
| `2026-09-05_retention_matched_control.json` | `2026-09-27_delta0_reference_mc.json` | One draw per budget, split equally over circuits; it sits at the 95.4 / 93.6 / 99.65 / 78.9th percentiles of the equal-split distribution. |
| `bootstrap_device_curve.json` | (none; not a test) | Resampling of the retained counts; upward-biased (its 50k 95% interval excludes the point estimate). |
| `_superseded/analytic_coverage_leak.json` | `2026-09-27_delta0_reference_mc.json` | Analytic shot-noise model; modelled false-positive rate 1.0; its σ is 34–195× smaller than the Monte Carlo spread at the device's raw shots. |

Two committed files are unchanged but need a reading note:
- `2026-08-24_within_sector_control.json`: `cases[0].weight_moved_frac = 0.494` is ‖δ‖₁/m₀, i.e. 24.7% of
  the weight transferred. The LP ranges and explicit counterexamples are in key `R3_within_sector_lp`.
- `2026-08-18_gausslaw_falsifier.json` is the bit-flip sweep; `fig_gausslaw` plots the ε-admixture curve
  of `2026-09-27_gauss_state.json`.

## Device data

`heron_counts_matched_L8_*.json` are the retained `ibm_fez` counts of this work (four budgets, 2026-08-29)
plus the Aer dry run; `m0_hat`/`m1_hat` in them are exact classical values, and only the support set is
device-derived. `heron_spectral.json` holds the companion's `L=6` reconstructed spectra; that run was
post-selected in reversed bit order, so all its retained determinants came from device errors, and its raw
counts are not deposited.

## `_superseded/`

The analytic coverage model (read by `delta0_reference_mc.py` only as the comparator for its spread
ratio), precursor-program data and the 2026-08-22 `.dat` copies. The 2026-09-05 text still quotes the
analytic model's σ (9.3×10⁻⁵ at 50k to 9.2×10⁻⁴ at 4k) as "shot noise is negligible"; the revision
replaces it with the Monte Carlo spread. No other manuscript number rests on these files. See
[`_superseded/README.md`](_superseded/README.md).
