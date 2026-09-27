# Reproduce the numbers and figures

Every non-hardware result is an exact classical computation. The core scripts need
`numpy`/`scipy`/`mpmath`; the re-analysis of the retained `ibm_fez` counts (including the reference Monte
Carlo for Δ₀) and the noise-model forecasts also need `qiskit`/`qiskit-aer`, and the `L=24` DMRG check
needs `physics-tenpy` (see `requirements.txt`). Scripts resolve their paths relative to `src/`, so
`python src/<script>.py` from the root and `cd src && python <script>.py` behave the same: JSON goes to
`data/`, plotted data to `paper/figs/*.dat`, and `git status` after a run shows whether a committed number
moved. The exception is the five manifest-hashed `nk_*` scripts, which write to `./06_results/` relative to
the working directory (see Sealed records). What cannot be regenerated yet is listed under
[Known gaps](#known-gaps).

The manuscript in `paper/` was **revised on 2026-09-27** against the phase-A recomputations (section
below); where the text and a data file named here disagree, the data file is the reference.

## Fastest check (seconds)

```bash
python src/verify.py        # or: make verify
```
Recomputes the current-probe moments `m₀, m₁, m₂` of the doped `L=6, U/t=4` ring, the shot-budget
thresholds `τ_k` at `N_s = 5×10⁴` (the per-circuit budget of the companion's `ibm_fez` run; assumed 2% bias
budget), and shows the joint `(m₀,m₁,m₂)`+Hankel battery **catching**, via `m₂`, a determinant truncation
that the lone first moment passes (found by an auto-calibrated scan: `d=97`). Since 2026-09-27 it orders the
truncation with the same deterministic `lexsort` key as the interval scripts and checks the quoted values:
`|Δm₁| = 0.186 < τ₁ = 0.229`, `|Δm₂| = 2.729 > τ₂ = 1.848`. Expected: **ALL CHECKS PASS**. (The earlier
`argsort` order printed 0.155 and 2.560 on the authors' machine; the cut keeps 7 of a 24-fold degenerate
`|ψ₀|²` shell, and the m₁-pass / m₂-reject verdict holds for all 300 tie choices tested in
`interval_battery.py`.)

## Make targets

| Target | Runs | Needs | Time |
|---|---|---|---|
| `make verify` | `verify.py` | numpy/scipy | seconds |
| `make cache` | rebuilds `src/cache/*_L12.npz` if missing (committed) | numpy/scipy | minutes |
| `make figures` | the `export_*` scripts, `separating_counts.py`, `interval_moment_mc.py` (list in the Makefile); the `run_*` demonstrations that also write `.dat` files run in the notebook | numpy/scipy/mpmath, matplotlib | ≈10 min |
| `make intervals` | `interval_moment.py`, `interval_battery.py`, `interval_moment_mc.py` | numpy/scipy | ≈3 min |
| `make blind` | `blind_score.py`, `blind_addenda.py`, `separating_counts.py` | numpy/scipy | ≈25 s |
| `make theory` | `small_checks.py --all` (six scripts, then its own keys) | numpy/scipy/mpmath | ≈15 min |
| `make gauss` | `export_gauss_state.py` | numpy/scipy | ≈25 s |
| `make revision` | `intervals` + `blind` + `theory` + `gauss` | numpy/scipy/mpmath | ≈20 min |
| `make device` | `delta0_reference_mc.py`, `export_device_dat.py` | qiskit (Statevector; no QPU) | 7–20 min |
| `make reproduce` | the notebook (phase-A scripts included; device Monte Carlo off by default) | + nbconvert, ipykernel | ≈45 min |
| `make paper` | `pdflatex` twice on `paper/main.tex` (committed `main.bbl`) | TeX | ≈1 min |

Times are from the phase-A runs on a laptop and vary with load (for example, `delta0_reference_mc.py` took
7.2 min on 2 BLAS threads and 19.3 min under concurrent load).

## Figure → engine → data

| Figure (paper) | What it shows | Engine / script | Data fragment(s) |
|---|---|---|---|
| `fig_hero` | the reconstruction-independent moment loop (schematic) | native TikZ | — |
| `fig_momentcone` | a moment sequence outside the moment cone | `export_momentcone_dat.py` (≈2 min) | `momentcone_point.dat` |
| `fig_akw` | `A(k,ω)` Mott map + per-`k` screen | `spectral_lanczos` → `export_akw_dat.py`; needs `src/cache/akw_L12.npz` | `akw_*.dat`, `akw_true/wrong.png` |
| `fig_sqw` | `S(q,ω)` gapped / `S^zz` gapless | `spectral_lanczos` → `export_sqw_dat.py`; needs `src/cache/sqw_L12.npz`. Since 2026-09-27 the rasters are placed at their true extents (`sqw_extent.dat`: q 0.1667–1.8333, charge 0–11t, spin 0–2.6t) and clipped by the 10t and 2.4t axis windows | `sqw_*.dat`, `dcp_*.dat`, `sqw_charge/spin.png` |
| `fig_collective_screen` | the screen on the collective channels | `run_collective_sumrule_falsifier.py` | `collective_screen_*.dat` |
| `fig_inverse_falsifier` | the negative-order (`m₋₁`, f-sum) test (relative shifts; the noise-normalized comparison is `--r8`) | `run_inverse_moment_falsifier.py` (without `--r8`; rewrites `data/2026-08-28_inverse_moment_falsifier.json`, reproduced bit for bit on 2026-09-27; ≈9 min) | `inverse_falsifier_*.dat` |
| `fig_teeth` | independence vs a circular control | `spectral_lanczos.run_teeth` → `export_teeth_dat.py`; needs `src/cache/teeth_L12.npz` | `teeth.dat` |
| `fig_teeth_shared` | the shared-state leak that still fires the screen | `run_teeth_shared.py` | `teeth_shared.dat` |
| `fig_separating` | diagonal vs off-diagonal reach (blinded record) | `separating_counts.py` from the sealed `data/blind_*.json` → also `data/2026-09-27_separating_counts.json`. The legend counts for the shaded region count in-region points only | `sep_{clean,krylov,trunc,ac_other,ac_pres_region,ac_chg_region,region_nonac}.dat` |
| `fig_device` | device coverage residual Δ₀ (real `ibm_fez`, `L=8`) against reference samplers; two panels: (a) linear, device vs the noiseless raw-matched band (mean ± 1 sd of 2000 replicas; the distribution is right-skewed, so this is not a 68% interval); (b) log scale, adding the uniform in-sector sampler | `delta0_reference_mc.py` → `data/2026-09-27_delta0_reference_mc.json` → `export_device_dat.py` | `device_points.dat`, `device_band.dat`, `device_uniform.dat` |
| `fig_heron` | the companion's `L=6` run on real IBM Heron data (post-selected in reversed bit order; see Hardware) | `export_heron_dat.py` (**2026-09-27**) from `data/heron_spectral.json`; reproduces the committed files byte for byte | `heron_{exact,hw}.dat` |
| `fig_circuit` | schematic of the simulated bond-basis estimator (not run on hardware) | `paper/figs/src/fig_circuit_qtk.tex` → `fig_circuit_qtk.pdf` (`pdflatex fig_circuit_qtk.tex`) | — |
| `fig_bracketing` + `fig_christoffel` | Gauss–Radau bracketing + the Christoffel width `Wₙ(t)` | `export_bracketing_dat.py`, `export_christoffel_dat.py` (≈3 min; exact `max_t Wₙ` by root-finding since 2026-09-27); the two-sided window bound: `markov_krein_window.py` | `bracket_*.dat`, `christoffel_*.dat` (`christoffel_maxW.dat` is read by `fig_christoffel` (b) and by the `fig_bracketing` inset) |
| `fig_momentmc` | Monte-Carlo interval-moment test at the `ibm_fez` budget; panel (b) now plots every truncation depth | `interval_moment_mc.py` | `momentmc_{hist,power,power_scan,scalars}.dat` |
| `fig_gausslaw` | U(1) Gauss-law state-level check | `export_gauss_state.py` → also `data/2026-09-27_gauss_state.json`. ε is an amplitude: the unphysical weight is ε²/((1−ε)²+ε²), 0.155 at ε = 0.3 | `gauss_state.dat` |

## Revision phase A (2026-09-27): scripts and outputs

Every phase-A script writes its JSON to a **new dated file** in `data/`; no sealed or committed data record
is overwritten, and the records they supersede stay in `data/` (see `data/README.md`). The figure `.dat`
files they regenerate replace the previous ones. Earlier versions of `fig_device.tex`, `fig_separating.tex`,
the `sep_*.dat` files and `christoffel_{maxW,poles}.dat` are kept in `paper/figs/_superseded/` (see its
README); the other changed figure files (`momentmc_power.dat`, `fig_momentmc.tex`, `fig_sqw.tex`,
`fig_christoffel.tex`, `fig_gausslaw.tex`) are in git history at commit `c4f0d77`, the last commit before
phase A. Item labels in the scripts (R1, B2, …) refer to the internal revision plan.

| Script (new or changed) | Command (from `src/`) | Output | Content | Deps, time |
|---|---|---|---|---|
| `delta0_reference_mc.py` (new) | `python delta0_reference_mc.py` (seed 20260927; env overrides in its docstring) | `data/2026-09-27_delta0_reference_mc.json` | Δ₀, \|S\|, Δ₁ and relative L1 of the device against: noiseless multinomial replicas at the device's raw shots (2000) and at its retained counts (per circuit and equal split); a uniform-noise mixture; a uniform in-sector sampler; percentile ranks; bit-order and Lanczos checks | qiskit, 7–20 min |
| `export_device_dat.py` (new) | `python export_device_dat.py` | `paper/figs/device_{points,band,uniform}.dat` | reformats the JSON above; no computation | seconds |
| `export_heron_dat.py` (new) | `python export_heron_dat.py` | `paper/figs/heron_{exact,hw}.dat` | reformats `data/heron_spectral.json`; byte-identical to the committed files | < 1 s |
| `check_sqw_extent.py` (new) | `python check_sqw_extent.py` | `data/2026-09-27_sqw_extent_check.json` | `fig_sqw`: displayed value below the charge-onset line for the current and the pre-2026-09-27 raster placement, weight clipped by the 10t / 2.4t windows, spinon-peak offsets (the phase-A scratch check, committed) | < 1 s |
| `verify.py` (changed) | `python verify.py` | prints only | deterministic `lexsort` truncation order; checks τ₀–τ₂ and the d=97 values quoted in the text | ≈10 s |
| `blind_addenda.py` (new) | `python blind_addenda.py` | `data/2026-09-27_blind_addenda.json` | post hoc addenda: (0) frozen rule rerun verbatim, 300/300; (a) primary endpoint, per-class rates, the two failed sealed predictions; (b) modelled FPR with no bias and ±2% bias; (c) exact rebuild of the 56 sealed truncations (stored tie order, re-validated on every run; abort otherwise) and a post hoc shared-state rescoring; (d) Krylov residuals; (e) calibrated scope | numpy/scipy, ≈11 s |
| `separating_counts.py` (extended) | `python separating_counts.py` | `paper/figs/sep_*.dat` (7 files), `data/2026-09-27_separating_counts.json` | the `fig_separating` recount from the sealed record, now exported | < 1 s |
| `interval_moment.py` (deterministic order) | `python interval_moment.py` | `data/2026-09-27_interval_moment_closure.json` | shot-budget closure on the 22-point grid; shell at each cut; comparison with the 2026-08-22 argsort run | ≈1 s |
| `interval_battery.py` (deterministic order, every-d scan) | `python interval_battery.py` | `data/2026-09-27_interval_battery.json` | joint battery at 0/2/4% bias for every d = 1…185; shell-boundary cuts; tie-break sensitivity (300 seeded draws per cut); comparison with `interval_battery_showcase.json` | ≈2 min |
| `interval_moment_mc.py` (deterministic order, every d) | `python interval_moment_mc.py` (seed 1, 4000 replicas) | `data/2026-09-27_interval_moment_mc.json`, `paper/figs/momentmc_{hist,power,power_scan,scalars}.dat` | power at every d and on the grid; ROC at d=98 and at the lowest-power cuts | ≈40 s |
| `export_gauss_state.py` (new) | `python export_gauss_state.py` | `data/2026-09-27_gauss_state.json`, `paper/figs/gauss_state.dat` (byte-identical to the committed file) | the ε-admixture curve ⟨ΣG²⟩ = g_u·w(ε) with an explicit selection rule (g_u = 2); diagnostic of the old rule | ≈6–25 s |
| `small_checks.py` (new) | `python small_checks.py --all` | `data/2026-09-27_theory_numerics.json` (all keys) | runs the six scripts below, then writes `R11_trace_over_dim`, `m4_psd_and_symmetry`, `M5_christoffel_claims`, `M6_inverse_moment_claim` | ≈15 min |
| `within_sector_lp.py` (new) | `python within_sector_lp.py` | keys `R3_within_sector_lp`, `M1_shifted_power_bound` | LP ranges of Δ_{K+1}, null-space dimensions, Vandermonde conditioning on raw and rescaled nodes, three explicit counterexamples, shifted-power bound checks | ≈16 s |
| `christoffel_tolerance_lp.py` (new) | `python christoffel_tolerance_lp.py` | key `R4_tolerance_inflated_atom` | extremal atom of positive measures within the tolerance box, U/t=4 and U/t=8 measures | ≈3 min |
| `necessary_sufficient_composition.py` (rewritten) | `python necessary_sufficient_composition.py` | key `R5_hausdorff_weighted_support` | Hankel lower and Hausdorff upper bounds on m₂ for four support intervals, a priori to oracle; no longer writes `data/necessary_sufficient_composition.json` (kept as the record) | ≈2.5 min |
| `estimator_form_offeigenstate.py` (new) | `python estimator_form_offeigenstate.py` | key `R6_estimator_forms_off_eigenstate` | first-moment estimator forms on truncated states, lexsort and argsort orders | ≈1 s |
| `run_inverse_moment_falsifier.py` (`--r8` added) | `python run_inverse_moment_falsifier.py --r8` | key `R8_inverse_moment_sensitivity` | noise-normalized sensitivity of m₋₁ against m₁, m₂ on the ring and the open chain; writes no figure file | ≈2.5 min |
| `export_christoffel_dat.py` (changed) | `python export_christoffel_dat.py` | `paper/figs/christoffel_*.dat`; keys `R9_christoffel_max_Wn`, `m1_radau_outside_hull`, `m2_rescaled_frame` | exact `max_t Wₙ` by root-finding; `christoffel_poles.dat` lists weighted atoms only | ≈3 min |

`small_checks.py --all` took 14–15 min in the phase-A runs (its docstring now says so).

Wording pass (2026-09-27, numbers unchanged). Docstrings, printed messages and provenance strings were
rescoped in `interval_battery.py`, `interval_moment.py`, `interval_moment_mc.py`, `hardware_blind_job.py`,
`bond_moment_estimator.py`, `estimator_scaling.py`, `run_sumrule_falsifier.py`, `run_falsifier_teeth.py`,
`run_gausslaw_falsifier.py`, `run_teeth_shared.py`, `run_inverse_moment_falsifier.py`, `blind_score.py`,
`delta0_reference_mc.py`, `small_checks.py` and `nk_stage2_device_job_v2_amendment.py` (no "remedy", no
"device-realizable", "device-measurable" or "same-sample" without qualification, relative shift kept apart
from detection power, no internal workflow identifiers, no "TREX"). The seven JSON records whose strings
changed were regenerated in a scratch copy and replaced only because every number reproduced bit for bit;
their earlier versions are in `data/_superseded/pre_wording_pass_2026-09-27/` (list in `data/README.md`).
`run_teeth_shared.py`'s record was not regenerated (Known gaps).

## Key numbers → script

| Number | Value | Where |
|---|---|---|
| current-probe moments `(m₀,m₁,m₂)` (L=6, U/t=4) | `0.5542, 4.2271, 33.3887` | `verify.py`, `interval_moment.py` |
| shot-budget 95% interval `z·δ₁` (50k shots, 2% bias) | `0.229` (5.4% of m₁) | `verify.py`, `data/2026-09-27_interval_moment_closure.json` |
| **sealed primary endpoint** | TPR `101/188 = 53.7%` [46.6, 60.7], FPR `1/112 = 0.9%` [0.2, 4.9]; nominal joint target 14.3% | `blind_score.py` → `data/2026-08-24_blind_harness_score.json`; `data/2026-09-27_blind_addenda.json` `a_primary_endpoint` |
| per-class rates | truncation `56/56`; spurious `45/61` (`25/27`+`20/34`); Krylov `0/71` (n_l = 1…4: 0/15, 0/20, 0/15, 0/21) | same |
| failed sealed secondary predictions | one-node Krylov `0/15`; m₀-preserving spurious `20/34` = 59% [42, 74] | `a_primary_endpoint.sealed_secondary_predictions_that_failed` |
| modelled FPR under the frozen rule (112 clean instances) | 1.3% unbiased; 8.0% with +2% bias; 11.5% with −2% | `b_fpr_sensitivity` |
| shared-state rescoring of the 56 sealed truncations (post hoc) | 56/56 rejected; every one fires on m₀ (min r₀ = 5.94) | `c_truncations_rebuild_and_shared_state` |
| one-node Krylov residual | m₂ off by −4.19% to −2.35%; max g₂/δ₂ = 0.986; 6/15 would fire at shot-only thresholds | `d_krylov` |
| calibrated scope | truncations keep 29–138 of 186 support configurations (75.1–98.7% of the \|ψ₀\|² weight); spurious w_s/m₀ = 2.1–24.9%; U/t ∈ {3,4,5,6} | `e_calibrated_scope` |
| `fig_separating` recount | 26 in the shaded region = 20 weight-preserving + 4 weight-changing spurious + truncation id 162 + clean id 184; m₁ fires on 16, m₂ alone on 10 | `separating_counts.py` → `data/2026-09-27_separating_counts.json` |
| same-sample simulation | ρ ≈ 0.43 at 84% coverage; no power gain at matched FPR | `lever1_flip_experiment.py` → `data/2026-09-01_lever1_oracle.json` `summary["84"]` |
| device Δ₀ (50k/30k/16k/4k) | 0.0044 / 0.032 / 0.123 / 0.373; coverage 94 / 85 / 71 / 35% | `data/2026-09-27_delta0_reference_mc.json` `device` |
| noiseless Δ₀ at the device's raw shots (2000 replicas) | 0.0361±0.0181 / 0.0660±0.0267 / 0.1208±0.0350 / 0.2691±0.0318; device percentile rank 0 / 1.7 / 57.9 / 100 | `budgets.<N>.noiseless_raw_matched` |
| noiseless Δ₀ at the device's retained counts | per circuit 0.224 / 0.279 / 0.343 / 0.443; equal split 0.150 / 0.205 / 0.273 / 0.404 | `noiseless_retained_matched`, `noiseless_retained_matched_equal_split` |
| committed single-draw retention-matched control | 0.203 / 0.248 / 0.329 / 0.418, at the 95.4 / 93.6 / 99.65 / 78.9th percentiles of the equal-split distribution | `data/2026-09-05_retention_matched_control.json`; rank in `retention_matched_control_seed1_rank` |
| uniform in-sector sampler at the device's kept counts | E[Δ₀] = 3.4×10⁻¹⁰ / 1.5×10⁻⁶ / 5.2×10⁻⁴ / 0.093; relative L1 < 4×10⁻⁹ (Lanczos floor) / < 4×10⁻⁹ / 0.014 / 0.90 | `uniform_in_sector_at_device_kept` |
| relative L1 of the reconstructed A(ω) | device 0.35 / 0.56 / 0.90 / 1.12; noiseless raw-matched at 50k 0.566 ± 0.032 (device below all 80 replicas) | `data/heron_counts_matched_L8_*.json` `falsifier.rel_L1`; `noiseless_raw_matched.rel_L1` |
| analytic σ vs Monte Carlo spread | MC sd / analytic σ = 195 / 177 / 125 / 34 (raw-matched comparator) | `sd_over_analytic_sigma` |
| `n_k` amendment gate | **DENIED (gate-negative)**: R1ro0.018 (ε = 3.5×10⁻³) UCB `0.02867` > `0.02667`; R2 (ε = 2.85×10⁻³) passes | `nk_stage0_gate_v4.py`, `nk_stage0_gate_v4_combine.py` → `data/2026-09-05_nk_v4_verdict.json` (see Sealed records) |
| truncation at d=98 (deterministic order, 2% bias) | \|g₁\| = 0.284 > 0.229 (m₁ rejects); over 300 tie choices m₁ rejects on 57%, the battery on 100% | `data/2026-09-27_interval_battery.json` `rows`, `tie_break_sensitivity` |
| tie-robust lone-m₁ miss caught by m₂ | d=97: \|g₁\| = 0.186 < 0.229, \|g₂\| = 2.729 > 1.848 (lexsort) | `full_scan_lexsort`, `full_scan_summary_lexsort.bias2.verify_py_rule_pick_d` |
| every-d scan (2% bias) | 17 lone-m₁ misses; m₂ catches 3 (d = 97, 96, 85); the whole battery passes 12 (d = 185, 184, 183, 165–160, 95, 87, 86). 0% bias: 13 and 6; 4% bias: 27 and 22. Over 300 tie choices the battery rejects d = 86 on 1.3%, d = 95 on 10.7%, d = 87 on 44.3%, d = 160 on 39.7%, d = 165 on 0.3% and d = 161–164, 183–185 on 0% | `full_scan_summary_lexsort`, `tie_break_sensitivity` |
| Monte-Carlo power of lone m₁ (shot-only threshold) | < 50% at 13 depths; minima 0.050 (d=162), 0.055 (d=86), 0.061 (d=95); FP rate 0.045; grid d=98 power 0.94 | `data/2026-09-27_interval_moment_mc.json` |
| no guaranteed power at m₀–m₂ | null-space dimension 2 (5 poles); a redistribution keeping m₀, m₁, m₂ exactly moves 71.7% of the weight; one keeping m₀, m₁ with \|Δ₂\| ≤ τ₂ moves 85.9%; both pass | `data/2026-09-27_theory_numerics.json` `R3_within_sector_lp` |
| LP range of Δ₂ over m₀,m₁-preserving redistributions | [−0.92, 7.09] (max 3.84 τ₂; contains 0) | `R3_within_sector_lp.lp_ranges` |
| Vandermonde conditioning | κ = 50.2 on nodes rescaled to [−1,1] (σ_min 0.056); 1.25×10⁶ on raw nodes (unit-dependent) | `R3_within_sector_lp.conditioning` |
| committed random K=1 case | ‖δ‖₁/m₀ = 0.494 (24.7% of the weight transferred); Δ₂ = 0.082 = 0.044 τ₂ | `within_sector_control.py` → `data/2026-08-24_within_sector_control.json` |
| shifted-power bound \|Δ_{s−1}\| ≤ ρ(D/2)^{s−1} | 0 violations in 20000 random cases; tight at s = 2 | `M1_shifted_power_bound` |
| tolerance-inflated extremal atom | 3.6–5.7× W₁ at t = 12–20 (U/t=4); 3.2–7.7× at t = 14–20 and 10–11× at the figure's t* (U/t=8) | `R4_tolerance_inflated_atom` |
| exact `max_t Wₙ`, n = 1…4 | 1 (at the mean) / 0.969 / 0.886 / 0.882 | `R9_christoffel_max_Wn`; `paper/figs/christoffel_maxW.dat` |
| Hankel / Hausdorff bounds on m₂ (L=6, U/t=8) | 123.3 ≤ m₂ (true 125.4); m₂ ≤ 147.6 on the weighted support [2.49, 13.92]t (oracle); a priori 579.8 (full Fock) or 243.3 (symmetry sector); 573.7 for [0.148, 52.21]t | `R5_hausdorff_weighted_support` |
| estimator forms off an eigenstate (44.4% coverage) | ⟨J(H−e₀)J⟩ = 4.558, ½⟨[J,[H,J]]⟩ = 3.717, m̄₁ = 3.743, exact 4.227 | `R6_estimator_forms_off_eigenstate` |
| m₋₁ sensitivity (L=12, U/t=8; 2% bias, spurious pole at 2.5t) | ring: f* = 0.197 (m₋₁), 0.051 (m₁), 0.042 (m₂); open chain: 0.168, 0.055, 0.042; ring cancellation D/(½⟨−T⟩) = 0.94 | `R8_inverse_moment_sensitivity`, `M6_inverse_moment_claim` |
| Tr(O)/dim for the L=8 addition moment | 6.80 (closed form U(1/4+(L−1)/8) − E₀/2, E₀ = −4.6035); O is not PSD (λ_min = −0.818) | `R11_trace_over_dim`, `m4_psd_and_symmetry` |
| Gauss-law curve | ⟨ΣG²⟩ = 2·w(ε), w(0.3) = 0.155 | `data/2026-09-27_gauss_state.json` |
| off-diagonal device bias floor (ibm_fez depth, 292 CZ) | 37.3 % (exact local, lower bound) to 115.6 % (global) of the signal, `≈ 37–116 %` | `offdiag_noise_forecast.py`, `offdiag_gsurface.py` → `data/2026-08-31_offdiag_gsurface.json` (`G = 292`) |
| `n_k` momentum-distribution forecast (FT-only bias) | `≈ 9 %` of the Fermi step | `nk_falsifier.py` |
| single-particle Mott edges / gap (L=12) | `μ± = ±2.484t`, `Δ ≈ 4.97t` | `spectral_lanczos` (`A(k,ω)`) |
| charge onset `Δc` (L=12) | `≈ 5.71t` | `spectral_lanczos` (structure factors) |

## Blinded pre-registration

The blinded test is sealed: `data/prereg.json` + its SHA-256 `data/prereg.sha256` (`0e4c3368…`) fix the
battery, thresholds, error catalogue and seed before any instance exists. `blind_generate.py` produces the
instances, `blind_classify.py`/`blind_score.py` apply the frozen battery, and
`data/blind_{instances_public,labels_sealed,verdicts}.json` hold the record. The generator and classifier
re-verify the seal and refuse to run if it is broken; the sealed files are read-only, so the writers refuse
to overwrite them (run in a scratch copy to regenerate).

What the sealed test did and did not do:
- The independent estimate `m̂` was simulated as the exact moment plus unbiased Gaussian noise
  (`blind_generate.py`); no same-sample covariance, no realized bias.
- The spurious feature is a point atom in the generator; the prereg text calls it a "Gaussian atom".
- The seal is self-attested (no third-party timestamp); the simulation is deterministic given the sealed
  seed, so the seal alone cannot exclude pre-seal iteration.
- Regenerating the truncation instances from the seed depends on the platform: the cut often falls inside
  a degenerate |ψ₀|² shell (51 of 56). `blind_addenda.py` rebuilds all 56 exactly from the stored recovered
  tie order, which every run re-validates against the sealed m̄, var_loc and m̂ (it aborts and writes
  nothing otherwise). The sealed files are the canonical record.

`blind_addenda.py` and `separating_counts.py` read only the sealed files and the frozen pipeline.
`data/separating_demonstration.json` (24 / 15 / 9) is kept unchanged; its correction is
`data/2026-09-27_separating_counts.json` (26 / 16 / 10).

## Sealed records: how to check them

`.gitattributes` keeps these files byte-exact on every platform (no CRLF conversion).

```bash
sha256sum data/prereg.json                   # = data/prereg.sha256                (0e4c3368…2792)
sha256sum data/manifest_nk_device_v1.json    # = data/manifest_nk_device_v1.sha256 (496203db…ae80)
cd data && sha256sum -c 2026-09-05_nk_v4_verdict.json.sha256   # 784927f0…325a
```

The `n_k` amendment verdict is not the raw output of the first combiner. On 2026-09-05 its `finding`
said that tightening readout "did not move the canary"; the ladder in the same file refutes that (the
canary fell 14.2% of the 18.4% needed). The same day the `finding` was rewritten by hand, the caveats on
the axis attribution and a `corrections` record (with the superseded digest) were added, and the file was
re-sealed. `nk_stage0_gate_v4_combine.py` now emits the corrected text; run in a scratch copy (with the
sealed verdict and its `.sha256` removed), it rewrites the sealed file byte for byte (`784927f0…325a`).

The manifest also records the SHA-256 of the analysis code and gate records it was sealed against
(`analysis_code_sha256`, `gate_record_sha256`). Those eight files (`src/nk_stage0_gate{,_v2,_v3}.py`,
`src/nk_stage0_v3R2.py`, `src/nk_falsifier.py`, `src/spectral_lanczos.py`,
`data/2026-09-02_nk_stage0_gate_v3.json`, `data/2026-09-02_nk_stage0_v3R2.json`) are kept byte-for-byte
as sealed, including their original line endings, so every recorded hash verifies; for the same reason
the five `nk_*` scripts among them still write to `./06_results/` (gitignored), not to `data/`.

## Hardware

**The `L=8` jobs of this work** (`ibm_fez`, 2026-08-29; `src/hardware_matched_job_L8.py`;
`data/heron_counts_matched_L8_*.json` at 50k/30k/16k/4k shots per circuit):
- Seven circuits prepare a product determinant with X gates (5 up electrons on even qubits, 4 down on odd
  qubits) and apply three layers of RZZ on-site terms and R_xx·R_yy nearest-neighbour hopping (periodic
  wrap), without Jordan–Wigner strings; circuit k evolves for time proportional to k, so k=0 is the bare
  determinant. They are number- and S_z-conserving samplers of the addition sector, not the Hubbard
  propagator of c†|ψ₀⟩.
- Bitstrings are read in the correct qubit order (`int(bs, 2)`) and post-selected to (N↑, N↓) = (5, 4);
  no configuration recovery. SamplerV2 ran with measurement twirling, Pauli gate twirling and XpXm dynamical
  decoupling; no readout-error mitigation was applied (SamplerV2 returns twirled raw counts).
- Only the support set S is device-derived. `m0_hat`/`m1_hat` in the records are exact classical values;
  Δ₀, Δ₁ and the relative L1 error are classical functions of S. No moment is estimated from counts.
- The reference distributions come from `delta0_reference_mc.py`, which draws from the exact noiseless
  output distributions of the same seven circuits (qiskit Statevector) and applies the same post-selection
  and falsifier arithmetic. Use its percentile ranks, not Gaussian z-scores. The single noiseless Aer seed
  (`aer_noiseless_baseline.json`) sits at the 19.5 / 17.7 / 81.4 / 98.3th percentiles of the raw-matched
  distribution, and `bootstrap_device_curve.py` resamples the retained counts (upward-biased; not a test).

**The companion's `L=6` run** (behind `fig_heron`) was post-selected in reversed qubit order, so all its
retained determinants came from device errors; its raw counts are not deposited.

The discriminating off-diagonal moments are evaluated in simulation only. The on-device `n_k` experiment
was pre-registered (manifest v1 sealed 2026-09-02), its single permitted amendment was denied by its sealed
gate on 2026-09-05, and it was never run (`nk_stage2_device_job.py` is the manifest-v1 job;
`nk_stage2_device_job_v2_amendment.py` is the denied amendment's job). Re-acquiring device data needs an
IBM Quantum account (use your own credentials; never commit a token).

## Known gaps

- **Manuscript.** `paper/main.tex` was revised on 2026-09-27 against the phase-A data; the length pass, the
  SciPost port and the companion-v3 dependence are open (`paper/README.md`).
- **`blind_addenda.py`, part (c):** the ARPACK default-start rebuild that recovered the sealed tie order
  was verified only with scipy 1.13.1 and numpy 1.26.4; it is now a diagnostic (reported under
  `arpack_default_start_rebuild_this_process`), and the rebuild uses the stored order. The shared-state
  thresholds reuse the sealed `var_loc` (centred on the exact E₀), and 4 of the 56 truncations have zero
  reconstructed weight (m̄ ≈ 10⁻³⁰), so their match is absolute, not relative; both are recorded in the JSON.
- **`run_teeth_shared.py` (L=12, `argsort` order):** a rerun (2026-09-27, scratch copy) reproduces the two 5%
  crossings (d ≈ 15977 and 43554, relative change < 2×10⁻⁶) but moves the three smallest-d points (m̄₁ by
  1.4%, 0.4% and 6.2% at d = 337, 246, 180) against the committed 2026-08-26 record, which is kept.
- **`fig_sqw` windows:** the display windows clip 6.3% of the broadened charge-grid weight (above 10t; 8.5% at
  the worst q) and 0.25% of the spin weight (above 2.4t); `src/check_sqw_extent.py` →
  `data/2026-09-27_sqw_extent_check.json` measures this and what the raster draws below the onset line
  (displayed value ≤ 0.19, i.e. ≤ 7% of vmax in linear scale, just below the onset; ≤ 0.05 below 5t).
- **`paper/figs/bracket_maxW.dat`** (from `export_bracketing_dat.py`, read by no figure) holds grid maxima
  (0.998659 / 0.947611 / 0.883640 / 0.868666) that differ from the exact `christoffel_maxW.dat`.
- **Provenance hash:** `delta0_reference_mc.py` hashes `data/_superseded/analytic_coverage_leak.json` as
  checked out; a CRLF checkout gives a different hash for the same content.
- **Notebook:** `notebooks/00_Reproduce_Everything.ipynb` runs the phase-A scripts too. It also runs
  `run_teeth_shared.py`, which moves `teeth_shared.dat` (above).
- **Unquoted extensions:** higher-order within-sector LPs (shifts of m₃, m₄ over battery-passing
  redistributions) were computed only outside this repository, during verification; add them to
  `within_sector_lp.py` before quoting them.
- **Legacy wording kept:** `data/2026-08-26_teeth_shared_state.json` (not regenerated, see above) still
  says "device-realizable"; `beyond_ed_dmrg_catch.py` (L=24 DMRG, not rerun) still calls its screen
  "device-realizable" and its record names an internal script path, as does
  `run_collective_sumrule_falsifier.py`'s; `lever1_flip_experiment.py` prints "make-or-break". The `nk_*`
  gate and seal scripts and their records keep their internal workflow identifiers, because the records are
  sealed and five of the scripts are hashed in the manifest.
- No Zenodo DOI yet; README thumbnails not regenerated.
