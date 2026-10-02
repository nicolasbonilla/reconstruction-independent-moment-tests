# Reproduce the numbers and figures

Every result other than the device sampling is an exact classical computation. The core scripts need
`numpy`/`scipy`/`mpmath`; the re-analysis of the `ibm_fez` records (including the reference Monte Carlo
for Δ₀) and the noise-model forecasts also need `qiskit`/`qiskit-aer` (the CZ count also
`qiskit-ibm-runtime`, for the `FakeFez` target), and the `L=24` DMRG check needs `physics-tenpy` (see
`requirements.txt`). Scripts resolve their paths relative to `src/`, so `python src/<script>.py` from the
root and `cd src && python <script>.py` behave the same: JSON goes to `data/`, plotted data to
`paper/figs/*.dat`. The exception is the five manifest-hashed `nk_*` scripts, which write to
`06_results/` at the repository root (`src/../06_results`, gitignored) and do not create it: run
`make nk-sealed-dir` (`mkdir -p 06_results`) first (see Sealed records).

After a run, `git diff --ignore-cr-at-eol --stat data paper/figs` lists the outputs that changed. Run time,
run date and environment fields change on every run; numbers should not. What does not regenerate exactly
is listed under [Known gaps](#known-gaps).

## Fastest check (seconds)

```bash
python src/verify.py        # or: make verify
```
Recomputes the current-probe moments `m₀, m₁, m₂` of the doped `L=6, U/t=4` ring, the shot-budget
thresholds `τ_k` at `N_s = 5×10⁴` (assumed 2% bias budget), and shows the joint `(m₀,m₁,m₂)`+Hankel battery
**catching**, via `m₂`, a determinant truncation that the lone first moment passes (auto-calibrated scan,
deterministic `lexsort` order: `d=97`, `|Δm₁| = 0.186 < τ₁ = 0.229`, `|Δm₂| = 2.729 > τ₂ = 1.848`). Its
record check then matches every number `paper/main.tex` prints from five records
(`data/2026-08-28_inverse_moment_falsifier.json`, `data/heron_spectral.json`, keys `R8`, `R9` and
`m2_rescaled_frame` of `data/2026-09-27_theory_numerics.json`, the forecast percentages of
`data/2026-08-31_offdiag_gsurface.json`, and `paper/figs/bracket_maxW.dat`) against
its key, at the printed precision. Expected: **ALL CHECKS PASS**. All other numbers: [Key numbers](#key-numbers--script).

## Make targets

| Target | Runs | Needs | Time |
|---|---|---|---|
| `make verify` | `verify.py` | numpy/scipy | seconds |
| `make cache` | rebuilds `src/cache/*_L12.npz` if missing (committed) | numpy/scipy | minutes |
| `make figures` | the `export_*` scripts, `separating_counts.py`, `interval_moment_mc.py`, `check_sqw_extent.py` (list in the Makefile); the `run_*` demonstrations that also write `.dat` files run in the notebook | numpy/scipy/mpmath, matplotlib | ≈10 min |
| `make intervals` | `interval_moment.py`, `interval_battery.py`, `interval_moment_mc.py` | numpy/scipy | ≈3 min |
| `make blind` | `blind_score.py`, `blind_addenda.py`, `separating_counts.py` | numpy/scipy | ≈25 s |
| `make theory` | `small_checks.py --all` (seven scripts incl. `full_support_rr.py`, then its own keys) | numpy/scipy/mpmath | ≈15 min |
| `make gauss` | `export_gauss_state.py` | numpy/scipy | ≈25 s |
| `make revision` | `intervals` + `blind` + `theory` + `gauss` | numpy/scipy/mpmath | ≈20 min |
| `make device` | `delta0_reference_mc.py`, `export_device_dat.py` | qiskit (Statevector; no QPU) | 7–30 min |
| `make hardware` | `sampler_cz_count.py`, `offdiag_gsurface.py`, `offdiag_noise_forecast.py`, `gate_tp_analysis.py` | qiskit, qiskit-aer, qiskit-ibm-runtime (no QPU) | ≈25–35 min |
| `make dmrg` | `beyond_ed_dmrg_catch.py` | physics-tenpy | ≈50 min |
| `make nkrecord` | `nk_record/nk_attempt_record.py --offline` | numpy | seconds |
| `make nk-sealed-dir` | `mkdir -p 06_results` (output folder of the five manifest-hashed `nk_*` scripts) | — | — |
| `make reproduce` | the notebook (device Monte Carlo off by default) | + nbconvert, ipykernel | ≈45–60 min |
| `make paper` | `pdflatex`, `bibtex`, `pdflatex` ×2 on `paper/main.tex` (a third pass if LaTeX asks for a rerun); fails on a LaTeX error or an undefined citation or reference | TeX | ≈1 min |

Times are from runs on a laptop and vary with load: `delta0_reference_mc.py` took 7.2 min on 2 BLAS threads
and 29.2 min (its 2026-10-02 regeneration) under concurrent load; `gate_tp_analysis.py` 9–12 min;
`beyond_ed_dmrg_catch.py` 48 min.

## Figure → engine → data

| Figure (paper) | What it shows | Engine / script | Data fragment(s) |
|---|---|---|---|
| `fig_hero` (Sec. 1) | one subspace state, two operator evaluations (schematic) | native TikZ | — |
| `fig_teeth_shared` (Sec. 3.1) | the shared-state leak under truncation (`L=12`), with the idealized curve and the circular control of the former `fig_teeth` | `run_teeth_shared.py`; `spectral_lanczos.run_teeth` → `export_teeth_dat.py` (needs `src/cache/teeth_L12.npz`) | `teeth_shared.dat`, `teeth.dat` |
| `fig_inverse_falsifier` (Sec. 3.2) | misplaced weight at fixed total and the negative-order (`m₋₁`, f-sum) moment (relative shifts; the noise-normalized comparison is `--r8`) | `run_inverse_moment_falsifier.py` (without `--r8`; rewrites `data/2026-08-28_inverse_moment_falsifier.json` and the figure files, reproduced bit for bit; ≈9 min) | `inverse_falsifier_*.dat`, `inverse_falsifier_table.tex` |
| `fig_separating` (Sec. 4) | diagonal vs off-diagonal reach (blinded record) | `separating_counts.py` from the sealed `data/blind_*.json` → also `data/2026-09-27_separating_counts.json`. The legend counts for the shaded region count in-region points only | `sep_{clean,krylov,trunc,ac_other,ac_pres_region,ac_chg_region,region_nonac}.dat` |
| `fig_device` (Sec. 6) | device coverage residual Δ₀ (`ibm_fez`, `L=8`) against reference samplers; (a) linear, device vs the noiseless raw-matched band (mean ± 1 sd of 2000 replicas; right-skewed, not a 68% interval); (b) log scale, adding the uniform in-sector sampler | `delta0_reference_mc.py` → `data/2026-09-27_delta0_reference_mc.json` → `export_device_dat.py` | `device_points.dat`, `device_band.dat`, `device_uniform.dat` |
| `fig_momentmc` (App. B) | Monte-Carlo interval-moment test at `N_s = 5×10⁴`; panel (b) plots every truncation depth | `interval_moment_mc.py` | `momentmc_{hist,power,power_scan,scalars}.dat` |
| `fig_akw` (App. D) | `A(k,ω)` Mott map + per-`k` screen | `spectral_lanczos` → `export_akw_dat.py`; needs `src/cache/akw_L12.npz` | `akw_*.dat`, `akw_true/wrong.png` |
| `fig_sqw` (App. D) | `S(q,ω)` gapped / `S^zz` gapless; rasters at their true extents (`sqw_extent.dat`: q 0.1667–1.8333, charge 0–11t, spin 0–2.6t), clipped by the 10t and 2.4t axis windows | `spectral_lanczos` → `export_sqw_dat.py`; needs `src/cache/sqw_L12.npz` | `sqw_*.dat`, `dcp_*.dat`, `sqw_charge/spin.png` |
| `fig_collective_screen` (App. D) | the test on the collective channels | `run_collective_sumrule_falsifier.py` | `collective_screen_*.dat` |
| `fig_momentcone` (App. E) | a moment sequence outside the moment cone | `export_momentcone_dat.py` (≈2 min) | `momentcone_point.dat` |
| `fig_bracketing` + `fig_christoffel` (App. F) | Gauss–Radau bracketing + the Christoffel width `Wₙ(t)` | `export_bracketing_dat.py`, `export_christoffel_dat.py` (≈3 min; exact `max_t Wₙ` by root-finding); the two-sided window bound: `markov_krein_window.py` | `bracket_*.dat`, `christoffel_*.dat` (`christoffel_maxW.dat` is read by `fig_christoffel` (b) and the `fig_bracketing` inset; the bracketing band includes the exact `Wₙ` peaks, so `bracket_maxW.dat` equals it) |
| `fig_gausslaw` (App. H) | U(1) Gauss-law state-level check | `export_gauss_state.py` → also `data/2026-09-27_gauss_state.json`. ε is an amplitude: the unphysical weight is ε²/((1−ε)²+ε²), 0.155 at ε = 0.3 | `gauss_state.dat` |
| `fig_circuit` (App. J) | schematic of the simulated bond-basis estimator (not run on hardware) | `paper/figs/src/fig_circuit_qtk.tex` → `fig_circuit_qtk.pdf` (`pdflatex fig_circuit_qtk.tex`) | — |
| `fig_heron` (not input by `main.tex`) | the companion's `L=6` run, kept as a source file | `export_heron_dat.py` from `data/heron_spectral.json`; reproduces the committed files byte for byte | `heron_{exact,hw}.dat` |

## Dated recomputations

Each recomputation writes its JSON to a **dated file** in `data/`; sealed records are never overwritten, and
the records a new file supersedes stay in `data/` (see `data/README.md`). Regenerated figure `.dat` files
replace the previous ones; earlier figure files are in `paper/figs/_superseded/` or in git history. Item
labels inside the scripts and records (R1, B2, H8, …) refer to internal review plans.

### 2026-09-27

| Script (new or changed) | Command (from `src/`) | Output | Content | Deps, time |
|---|---|---|---|---|
| `delta0_reference_mc.py` (new) | `python delta0_reference_mc.py` (seed 20260927; env overrides in its docstring) | `data/2026-09-27_delta0_reference_mc.json` | Δ₀, \|S\|, Δ₁ and relative L1 of the device against: noiseless multinomial replicas at the device's raw shots (2000) and at its retained counts (per circuit and equal split); a uniform-noise mixture with the exact E[Δ₀](p) and its shape (`monotonicity_detail`); a uniform in-sector sampler; percentile ranks; bit-order and Lanczos checks. Input hashes over LF-normalised bytes. Regenerated 2026-10-02 (every Monte Carlo and analytic value bit-identical; monotonicity flags now over all rows sorted by p) | qiskit, 7–30 min |
| `export_device_dat.py` (new) | `python export_device_dat.py` | `paper/figs/device_{points,band,uniform}.dat` | reformats the JSON above | seconds |
| `export_heron_dat.py` (new) | `python export_heron_dat.py` | `paper/figs/heron_{exact,hw}.dat` | reformats `data/heron_spectral.json`; byte-identical to the committed files | < 1 s |
| `check_sqw_extent.py` (new) | `python check_sqw_extent.py` | `data/2026-09-27_sqw_extent_check.json` | `fig_sqw`: displayed value below the charge-onset line, weight clipped by the 10t / 2.4t windows, spinon-peak offsets | < 1 s |
| `blind_addenda.py` (new; extended 2026-10-02) | `python blind_addenda.py` | `data/2026-09-27_blind_addenda.json` | post hoc addenda: (0) frozen rule rerun verbatim, 300/300; (a) primary endpoint, per-class rates, the two failed sealed predictions; (b) modelled FPR with no bias and ±2% bias; (c) exact rebuild of the 56 sealed truncations (stored tie order, re-validated on every run; abort otherwise) and the shared-state rescoring; (d) Krylov residuals; (e) calibrated scope; (f) the noise scale of the sealed text and variance inflation f = 1, 2, 4, 10; (g) shared-state identity on the 112 faithful controls. Hashes over LF-normalised bytes; the tie order is hashed as an array | numpy/scipy, ≈11 s |
| `separating_counts.py` (extended) | `python separating_counts.py` | `paper/figs/sep_*.dat` (7 files), `data/2026-09-27_separating_counts.json` | the `fig_separating` recount from the sealed record | < 1 s |
| `interval_moment.py` (deterministic order) | `python interval_moment.py` | `data/2026-09-27_interval_moment_closure.json` | shot-budget closure on the 22-point grid; shell at each cut; comparison with the 2026-08-22 argsort run | ≈1 s |
| `interval_battery.py` (deterministic order, every-d scan) | `python interval_battery.py` | `data/2026-09-27_interval_battery.json` | joint battery at 0/2/4% bias for every d = 1…185; shell-boundary cuts; tie-break sensitivity (300 seeded draws per cut) | ≈2 min |
| `interval_moment_mc.py` (deterministic order, every d) | `python interval_moment_mc.py` (seed 1, 4000 replicas) | `data/2026-09-27_interval_moment_mc.json`, `paper/figs/momentmc_*.dat` | power at every d and on the grid; ROC at d=98 and at the lowest-power cuts | ≈40 s |
| `export_gauss_state.py` (new) | `python export_gauss_state.py` | `data/2026-09-27_gauss_state.json`, `paper/figs/gauss_state.dat` | the ε-admixture curve ⟨ΣG²⟩ = g_u·w(ε) with an explicit selection rule (g_u = 2) | ≈6–25 s |
| `small_checks.py` (new; extended 2026-10-02) | `python small_checks.py --all` | `data/2026-09-27_theory_numerics.json` (all keys) | runs the seven scripts below, then writes `R11_trace_over_dim`, `m4_psd_and_symmetry`, `M5_christoffel_claims`, `M6_inverse_moment_claim`, `MA6_estimator_shape_by_order`, `B2_one_node_m2_deficit`, `SOTA1_02_W1_transport` | ≈15 min |
| `within_sector_lp.py` (new) | `python within_sector_lp.py` | keys `R3_within_sector_lp`, `M1_shifted_power_bound` | LP ranges of Δ_{K+1} (and Δ_{K+2} at the optimum), null-space dimensions, Vandermonde conditioning, three explicit counterexamples, shifted-power bound with the sharp Chebyshev constant | ≈16 s |
| `christoffel_tolerance_lp.py` (new) | `python christoffel_tolerance_lp.py` | key `R4_tolerance_inflated_atom` | extremal atom of positive measures within the tolerance box, U/t=4 and U/t=8 measures | ≈3 min |
| `necessary_sufficient_composition.py` (rewritten) | `python necessary_sufficient_composition.py` | key `R5_hausdorff_weighted_support` | Hankel lower and Hausdorff upper bounds on m₂ for four support intervals, a priori to oracle | ≈2.5 min |
| `estimator_form_offeigenstate.py` (new) | `python estimator_form_offeigenstate.py` | key `R6_estimator_forms_off_eigenstate` | first-moment estimator forms on truncated states | ≈1 s |
| `run_inverse_moment_falsifier.py` (`--r8` added) | `python run_inverse_moment_falsifier.py --r8` | key `R8_inverse_moment_sensitivity` | noise-normalized sensitivity of m₋₁ against m₁, m₂ on the ring and the open chain; writes no figure file | ≈2.5 min |
| `export_christoffel_dat.py` (changed) | `python export_christoffel_dat.py` | `paper/figs/christoffel_*.dat`; keys `R9_christoffel_max_Wn`, `m1_radau_outside_hull`, `m2_rescaled_frame` | exact `max_t Wₙ` by root-finding; `christoffel_poles.dat` lists weighted atoms only | ≈3 min |
| `verify.py` (changed) | `python verify.py` | prints only | deterministic `lexsort` order; record check of the numbers main.tex prints from five records | ≈10 s |

Wording pass (2026-09-27, numbers unchanged): docstrings, printed messages and provenance strings were
rescoped in fifteen scripts; the seven JSON records whose strings changed were regenerated in a scratch copy
and replaced only because every number reproduced bit for bit. Their earlier versions are in
`data/_superseded/pre_wording_pass_2026-09-27/` (list in `data/README.md`).

### 2026-10-02

| Script (new or changed) | Command | Output | Content | Deps, time |
|---|---|---|---|---|
| `full_support_rr.py` (new) | `cd src && python full_support_rr.py` | `data/2026-10-02_full_support_rr.json` | Rayleigh–Ritz on the full ground-state support of the calibration ring at U/t = 3, 4, 5, 6, 8: moments k = 0…4 of both spectra, the leak form, Δ_k/τ_k and the default verdict, relative L1 at η = 0.15t on three windows, poles, Wasserstein-1; the one-node Krylov deficit; and the deterministic leak scores of every truncation of the U/t=4 support with Spearman correlations against L1 | numpy/scipy, ≈10 s |
| `sampler_cz_count.py` (new) | `cd src && python sampler_cz_count.py` | `data/sampler_cz_count.json` | transpiled CZ counts of the sampler circuits: (A) `build_circuit(k)` of both job scripts transpiled to `FakeFez`, seeds 0–4; (B) the ISA circuits actually submitted in the six `ibm_fez` jobs (encoded pubs stored in the JSON, so the re-count runs offline), identified by simulation | qiskit, qiskit-ibm-runtime, 2–4 min |
| `offdiag_gsurface.py`, `offdiag_noise_forecast.py` (reference depth read from `sampler_cz_count.json`) | `cd src && python offdiag_gsurface.py && python offdiag_noise_forecast.py` | `data/2026-08-31_offdiag_gsurface.json`, `data/2026-08-30_offdiag_noise_forecast.json` (fixed names; `run_date` inside) | first-moment bias floor of the surrogate density probe at G = 295 CZ, with 213 and 292 kept as grid points | qiskit, qiskit-aer; the forecast ≈9 min |
| `beyond_ed_dmrg_catch.py` (degeneracy handling) | `cd src && python beyond_ed_dmrg_catch.py` | `data/2026-10-02_beyond_ed_dmrg_catch.json` (supersedes `data/2026-08-28_beyond_ed_dmrg_catch.json`) | the shared-state test on sampled subspaces at L = 8 (ED), 20 and 24 (DMRG); a degenerate subspace ground level is reported as a range over the ground manifold and left out of the summaries | physics-tenpy, ≈50 min |
| `gate_tp_analysis.py` (new) | `python src/gate_tp_analysis.py` | `data/2026-10-02_gatetp_reanalysis.json` | the two 2026-08-21 `L=6` jobs: decode of the raw records, the submitting notebook's pipeline in both bit orders, a noiseless bit-order check, and what a support-based moment test can see when the support is complete | qiskit (qiskit-ibm-runtime optional), 9–12 min |
| `nk_record/nk_attempt_record.py` (new) | `python src/nk_record/nk_attempt_record.py --offline` | `data/2026-09-05_nk_device_attempt_record.json` | the sealed `n_k` device attempt: manifest hash, the deposited scripts as run, the logged outputs with UTC timestamps, the historical calibrations re-queried at the logged instants with the sealed chain rule applied, the account's jobs after the seal | numpy (live mode: qiskit-ibm-runtime, read-only), seconds |
| `export_bracketing_dat.py` (changed) | `cd src && python export_bracketing_dat.py` | `paper/figs/bracket_{env_n2,env_n4,maxW}.dat` | the bracketing band now includes the exact `Wₙ` peaks (0.968792 for n = 2, 0.882243 for n = 4) | minutes |
| `run_inverse_moment_falsifier.py` (claim string) | as above | `data/2026-08-28_inverse_moment_falsifier.json` | claim and provenance strings only (on the ring m₋₁ moves most in relative terms; on the open chain least); every number and figure file reproduced byte for byte | ≈9 min |
| `data/heron_spectral.json` | copied from the companion's deposit | — | arrays unchanged; the provenance corrected by the companion (max \|A_hw − A_exact\| = 7.9×10⁻¹⁵, not a byte-for-byte identity; retained, not recovered, support); `export_heron_dat.py` reproduces the `.dat` files byte for byte | — |

## Key numbers → script

| Number | Value | Where |
|---|---|---|
| current-probe moments `(m₀,m₁,m₂)` (L=6, U/t=4) | `0.5542, 4.2271, 33.3887` | `verify.py`, `interval_moment.py` |
| shot-budget 95% interval `z·δ₁` (50k shots, 2% bias) | `0.229` (5.4% of m₁) | `verify.py`, `data/2026-09-27_interval_moment_closure.json` |
| estimator forms off an eigenstate (44.4% coverage, d=100) | ⟨J(H−e₀)J⟩ = 4.558, ½⟨[J,[H,J]]⟩ = 3.717, m̄₁ = 3.743, exact 4.227 | key `R6_estimator_forms_off_eigenstate` |
| **blinded test, as implemented (sealed primary endpoint)** | TPR `101/188 = 53.7%` [46.6, 60.7], FPR `1/112 = 0.9%` [0.2, 4.9]; nominal joint target 14.3% | `blind_score.py` → `data/2026-08-24_blind_harness_score.json`; `data/2026-09-27_blind_addenda.json` `a_primary_endpoint` |
| per-class rates | truncation `56/56`; spurious `45/61` (`25/27` weight-changing + `20/34` m₀-preserving); Krylov `0/71` (n_l = 1…4: 0/15, 0/20, 0/15, 0/21) | same |
| failed sealed secondary predictions | one-node Krylov `0/15`; m₀-preserving spurious `20/34` = 59% [42, 74] | `a_primary_endpoint.sealed_secondary_predictions_that_failed` |
| **noise scale of the sealed text** (same draws × 1.96, threshold unchanged; post hoc) | FPR `40/112 = 35.7%` [27.4, 44.9]; TPR `126/188 = 67.0%`, including 17/56 Krylov n_l ≥ 2 rejected by noise alone; expected FPR over fresh noise 34.4% | `f_noise_model_sensitivity.literal_sealed_model` |
| variance inflation f (noise and threshold) | f = 4: m₀-preserving 17/34, weight-changing 22/27, FPR 15/112 (expected 7.1%); f = 10: 15/34, 19/27, 19/112 (10.7%); truncations 56/56 throughout | `f_noise_model_sensitivity.variance_inflation` |
| modelled FPR under the frozen rule (112 clean instances) | 1.3% unbiased; 8.0% with +2% bias; 11.5% with −2% | `b_fpr_sensitivity` |
| shared-state rescoring of the 56 sealed truncations (post hoc) | leak score above 10⁻⁶ on m₀, m₁, m₂ for all 56 (smallest Δ₀/m̂₀ = 0.25); with the sealed draws and τ_k as a convention 56/56 rejected, every one on m₀ (min r₀ = 5.94) | `c_truncations_rebuild_and_shared_state` |
| shared-state estimate on the 112 faithful controls | equal to the oracle-centred one (max relative difference 2.6×10⁻¹⁵); FPR 1/112 (40/112 at the sealed-text noise scale) | `g_shared_state_faithful_controls` |
| one-node Krylov residual | m₂ off by 2.35–4.19%; passes because τ₂ includes the 2% bias term (τ₂ ≥ 1.96 × 2% of m₂); largest residual 0.986 τ₂; 6/15 would fire at shot-only thresholds | `d_krylov`; `data/2026-10-02_full_support_rr.json` `summary.one_node_deficit_rel_m2_by_U` |
| calibrated scope | truncations keep 29–138 of 186 support configurations (75.1–98.7% of the \|ψ₀\|² weight); spurious w_s/m₀ = 2.1–24.9%; U/t ∈ {3,4,5,6} | `e_calibrated_scope` |
| `fig_separating` recount | 26 in the shaded region = 20 weight-preserving + 4 weight-changing spurious + truncation id 162 + clean id 184; m₁ fires on 16, m₂ alone on 10 | `separating_counts.py` → `data/2026-09-27_separating_counts.json` |
| same-sample simulation | ρ ≈ 0.43 at 84% coverage; no power gain at matched FPR | `lever1_flip_experiment.py` → `data/2026-09-01_lever1_oracle.json` `summary["84"]` |
| full-support Rayleigh–Ritz (U/t=4) | Δ₀, Δ₁ < 10⁻¹⁴ relative; Δ₂ = 0.62 (1.9% of m₂ = 33.39) = 0.34 τ₂, passes; relative L1 at η = 0.15t 0.93 (0.930 wide window, 0.936 on [−2, 25]t, 0.929 window-free) | `data/2026-10-02_full_support_rr.json` `summary` |
| full-support RR over U/t = 3–8 | Δ₂ = 0.20–0.36 τ₂, relative L1 0.68–0.99; passes at every U/t | same |
| leak score vs L1 error (U/t=4 truncations, d = 7…185, non-degenerate) | Spearman 0.44 / 0.43 / 0.31 for k = 0 / 1 / 2 | `truncation_scan_U4.spearman_truncations_only` |
| device Δ₀ (50k/30k/16k/4k) | 0.0044 / 0.032 / 0.123 / 0.373; coverage 94 / 85 / 71 / 35% | `data/2026-09-27_delta0_reference_mc.json` `device` |
| noiseless Δ₀ at the device's raw shots (2000 replicas) | 0.0361±0.0181 / 0.0660±0.0267 / 0.1208±0.0350 / 0.2691±0.0318; device percentile rank 0 / 1.7 / 57.9 / 100 | `budgets.<N>.noiseless_raw_matched` |
| uniform-noise mixture, exact E[Δ₀](p) | not monotone at any budget: one minimum, then a rise; at 50k 0.036 at p = 0, minimum 0.0020 near p ≈ 0.93, 0.0024 at p = 1 (minima: p ≈ 0.88 at 30k, 0.77 at 16k, 0.09 at 4k) | `budgets.<N>.uniform_noise_mixture_raw_matched.monotonicity_detail` |
| noiseless Δ₀ at the device's retained counts | per circuit 0.224 / 0.279 / 0.343 / 0.443; equal split 0.150 / 0.205 / 0.273 / 0.404 | `noiseless_retained_matched`, `noiseless_retained_matched_equal_split` |
| committed single-draw retention-matched control | 0.203 / 0.248 / 0.329 / 0.418, at the 95.4 / 93.6 / 99.65 / 78.9th percentiles of the equal-split distribution | `data/2026-09-05_retention_matched_control.json`; rank in `retention_matched_control_seed1_rank` |
| uniform in-sector sampler at the device's kept counts | E[Δ₀] = 3.4×10⁻¹⁰ / 1.5×10⁻⁶ / 5.2×10⁻⁴ / 0.093; relative L1 < 4×10⁻⁹ (Lanczos floor) / < 4×10⁻⁹ / 0.014 / 0.90 | `uniform_in_sector_at_device_kept` |
| relative L1 of the reconstructed A(ω) | device 0.35 / 0.56 / 0.90 / 1.12; noiseless raw-matched at 50k 0.566 ± 0.032 (device below all 80 replicas) | `data/heron_counts_matched_L8_*.json` `falsifier.rel_L1`; `noiseless_raw_matched.rel_L1` |
| analytic σ vs Monte Carlo spread (superseded model) | MC sd / analytic σ = 195 / 177 / 125 / 34 | `sd_over_analytic_sigma` |
| CZ count of the sampler circuits on `ibm_fez` | median 295 over the 24 deployed `L=8` circuits with k ≥ 1 (range 278–313); 213 for the 12-qubit (`L=6`) circuit family | `sampler_cz_count.py` → `data/sampler_cz_count.json` `headline` |
| off-diagonal bias floor of a surrogate density probe (2.5×10⁻³ two-qubit error, no mitigation) | at 295 CZ: 37.6% (exact static local depolarizing) to 116.4% (global) of the signal; coherent ZZ 46.9%; at 100 CZ: 13.5% / 49.4% | `offdiag_gsurface.py` → `data/2026-08-31_offdiag_gsurface.json`; `data/2026-08-30_offdiag_noise_forecast.json` `at_ibm_fez_anchor` |
| `n_k` momentum-distribution forecast (FT-only bias) | `≈ 9 %` of the Fermi step | `nk_falsifier.py` |
| `n_k` amendment gate | **DENIED**: R1ro0.018 (ε = 3.5×10⁻³) UCB `0.02867` > `0.02667`; R2 (ε = 2.85×10⁻³) passes | `nk_stage0_gate_v4.py`, `nk_stage0_gate_v4_combine.py` → `data/2026-09-05_nk_v4_verdict.json` |
| `n_k` device attempt | 0 jobs, 0 shots; every re-queried calibration snapshot forfeits under the sealed v1 chain rule (T1 ≥ 150 µs, T2 ≥ 100 µs on a 12-qubit chain) | `data/2026-09-05_nk_device_attempt_record.json` `summary`, `requery` |
| 2026-08-21 `L=6` jobs | retained support 300/300 (the whole sector) in both arms under either bit order; retention in the correct order 29.9% / 30.0% (73 364 and 73 475 of 245 000) | `data/2026-10-02_gatetp_reanalysis.json` `summary`, `bitorder` |
| `L=24` DMRG check (U/t=8, open chain) | from 1000 to 15000 shots Δ₀ falls from 40.6 to 24.5 (DMRG m₀ = 55.0); relative error of m̄₀ 0.89 → 0.55; m̂₀ 10.1–15.5% below m₀; at 500 shots the subspace ground level is 2-fold degenerate and Δ₀ spans 40.1–44.5 over it | `data/2026-10-02_beyond_ed_dmrg_catch.json` `summary`, `runs` |
| truncation at d=98 (deterministic order, 2% bias) | \|g₁\| = 0.284 > 0.229 (m₁ rejects); over 300 tie choices m₁ rejects on 57%, the battery on 100% | `data/2026-09-27_interval_battery.json` `rows`, `tie_break_sensitivity` |
| tie-robust lone-m₁ miss caught by m₂ | d=97: \|g₁\| = 0.186 < 0.229, \|g₂\| = 2.729 > 1.848 (lexsort) | `full_scan_lexsort`, `full_scan_summary_lexsort.bias2.verify_py_rule_pick_d` |
| every-d scan (2% bias) | 17 lone-m₁ misses; m₂ catches 3 (d = 97, 96, 85); the whole battery passes 12 (d = 185, 184, 183, 165–160, 95, 87, 86). 0% bias: 13 and 6; 4% bias: 27 and 22. Over 300 tie choices the battery rejects d = 86 on 1.3%, d = 95 on 10.7%, d = 87 on 44.3%, d = 160 on 39.7%, d = 165 on 0.3% and d = 161–164, 183–185 on 0% | `full_scan_summary_lexsort`, `tie_break_sensitivity` |
| Monte-Carlo power of lone m₁ (shot-only threshold) | < 50% at 13 depths; minima 0.050 (d=162), 0.055 (d=86), 0.061 (d=95); FP rate 0.045; grid d=98 power 0.94 | `data/2026-09-27_interval_moment_mc.json` |
| no guaranteed power at m₀–m₂ | null-space dimension 2 (5 poles); a redistribution keeping m₀, m₁, m₂ exactly moves 71.7% of the weight; one keeping m₀, m₁ with \|Δ₂\| ≤ τ₂ moves 85.9%; both pass | `data/2026-09-27_theory_numerics.json` `R3_within_sector_lp` |
| LP range of Δ₂ over m₀,m₁-preserving redistributions | [−0.92, 7.09] (max 3.84 τ₂; contains 0); over m₀–m₂-preserving ones, the redistribution that maximizes Δ₃ moves m₄ by 1.186 τ₄ (unique optimum) | `R3_within_sector_lp.lp_ranges` (`next_order_at_optimum`) |
| Vandermonde conditioning | κ = 50.2 on nodes rescaled to [−1,1] (σ_min 0.056); 1.25×10⁶ on raw nodes (unit-dependent) | `R3_within_sector_lp.conditioning` |
| committed random K=1 case | ‖δ‖₁/m₀ = 0.494 (24.7% of the weight transferred); Δ₂ = 0.082 = 0.044 τ₂ | `within_sector_control.py` → `data/2026-08-24_within_sector_control.json` |
| shifted-power bound with the sharp Chebyshev constant, \|Δ_{s−1}\| ≤ 2ρ(D/4)^{s−1} | 0 violations in 20000 random cases; attained at s = 2 (largest ratio 1.000), 0.998 at s = 3 | `M1_shifted_power_bound` (`all_cases_hold_sharp`, `random_null_space_mc`) |
| tolerance-inflated extremal atom | 3.6–5.7× W₁ at t = 12–20 (U/t=4); 3.2–7.7× at t = 14–20 and 10–11× at the figure's t* (U/t=8) | `R4_tolerance_inflated_atom` |
| exact `max_t Wₙ`, n = 1…4 | 1 (at the mean) / 0.969 / 0.886 / 0.882 | `R9_christoffel_max_Wn`; `paper/figs/christoffel_maxW.dat` |
| Hankel / Hausdorff bounds on m₂ (L=6, U/t=8) | 123.3 ≤ m₂ (true 125.4); m₂ ≤ 147.6 on the weighted support [2.49, 13.92]t (oracle); a priori 579.8 (full Fock) or 243.3 (symmetry sector); 573.7 for [0.148, 52.21]t | `R5_hausdorff_weighted_support` |
| m₋₁ sensitivity (L=12, U/t=8; 2% bias, spurious pole at 2.5t) | default rule (shot term + 2% bias): ring f* = 0.198 (m₋₁), 0.057 (m₁), 0.046 (m₂); open chain 0.169, 0.059, 0.045 (bias-only: ring 0.197 / 0.051 / 0.042, open chain 0.168 / 0.055 / 0.042); ring cancellation D/(½⟨−T⟩) = 0.94 | `R8_inverse_moment_sensitivity.*.f_star_full_rule`, `.f_star_bias_only`; `M6_inverse_moment_claim` |
| local-estimator shape (L=6, U/t=4) | skewness and excess kurtosis of O_k^loc for k = 0…4 and of its 5×10⁴-shot mean (m₁: −0.011 and 3.5×10⁻⁴; the Monte Carlo values of Fig. momentmc, −0.04 and −0.05, agree within their sampling errors) | `MA6_estimator_shape_by_order` |
| Tr(O)/dim for the L=8 addition moment | 6.80 (closed form U(1/4+(L−1)/8) − E₀/2, E₀ = −4.6035); O is not PSD (λ_min = −0.818) | `R11_trace_over_dim`, `m4_psd_and_symmetry` |
| Gauss-law curve | ⟨ΣG²⟩ = 2·w(ε), w(0.3) = 0.155 | `data/2026-09-27_gauss_state.json` |
| single-particle Mott edges / gap (L=12) | `μ± = ±2.484t`, `Δ ≈ 4.97t` | `spectral_lanczos` (`A(k,ω)`) |
| charge onset `Δc` (L=12) | `≈ 5.71t` | `spectral_lanczos` (structure factors) |

## Blinded pre-registration

The blinded test is sealed: `data/prereg.json` + its SHA-256 `data/prereg.sha256` (`0e4c3368…`) fix the
battery, the decision rule and its threshold inputs (`N_s`, `z`, `b`), the error catalogue, the primary
endpoint, the secondary predictions and the seed before any instance exists. The formula for `τ_k` is fixed
in the classifier code, not in the hashed text. `blind_generate.py` produces the instances,
`blind_classify.py`/`blind_score.py` apply the frozen battery, and
`data/blind_{instances_public,labels_sealed,verdicts}.json` hold the record. The generator and classifier
re-verify the seal and refuse to run if it is broken; the sealed files are read-only, so the writers refuse
to overwrite them (run in a scratch copy to regenerate).

- The independent estimate `m̂` was simulated as the exact moment plus unbiased Gaussian noise
  (`blind_generate.py`), with the standard error the threshold's shot term assumes; the sealed text
  specifies a 1.96-fold larger noise scale (addenda part f).
- The spurious feature is a point atom in the generator; the prereg text calls it a "Gaussian atom".
- The seal is self-attested (no third-party timestamp). Sealing, generation, classification and scoring
  ran in a single invocation, so the code predates the seal but is not covered by it. The scripts in `src/`
  differ from the 2026-08-24 working copies that wrote the sealed files only in the output folder
  (`06_results/` then, `data/` now), the refuse-to-overwrite guards and, in `blind_score.py`, the printed
  summary. SHA-256 of those working copies (not deposited): `blind_preregister.py`
  `ff29c3371b11829f9246cfbbdfa05499bb0c4d7da19b4e7ff831fee1fe1b8d79`, `blind_generate.py`
  `de7f299aa46f616f384fa326d7b2c67cc9a46e37478278a193ddee83836b02a7`, `blind_classify.py`
  `0c14b1c9c87669877bcf8fdfaa0f34fbccfa39818cab8442da29453628d540f2`, `blind_score.py`
  `b98868d24cd62b2ce8fea5a2e6a66c1b174bfcdfb0f04c9ff5401d5eaedf8f52`.
- Regenerating the truncation instances from the seed depends on tie-breaking inside degenerate |ψ₀|²
  shells (the cut falls inside one in 51 of 56). `blind_addenda.py` rebuilds all 56 exactly from the stored
  tie order, which every run re-validates against the sealed m̄, var_loc and m̂ (it aborts and writes nothing
  otherwise). The sealed files are the canonical record.

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

The amendment verdict file was corrected and re-sealed on the day of the verdict; the record is in
[`PROCESS_RECORD.md`](PROCESS_RECORD.md), section 3. `nk_stage0_gate_v4_combine.py` emits the corrected
text; run in a scratch copy (with the sealed verdict and its `.sha256` removed), it rewrites the sealed file
byte for byte (`784927f0…325a`).

The manifest also records the SHA-256 of the analysis code and gate records it was sealed against
(`analysis_code_sha256`, `gate_record_sha256`). Those eight files (`src/nk_stage0_gate{,_v2,_v3}.py`,
`src/nk_stage0_v3R2.py`, `src/nk_falsifier.py`, `src/spectral_lanczos.py`,
`data/2026-09-02_nk_stage0_gate_v3.json`, `data/2026-09-02_nk_stage0_v3R2.json`) are kept byte-for-byte
as sealed, including their original line endings, so every recorded hash verifies; for the same reason the
five `nk_*` scripts among them still write to `06_results/` at the repository root (create it first:
`make nk-sealed-dir`), not to `data/`.

## Hardware

**The `L=8` jobs** (`ibm_fez`, created 2026-08-29 UTC; `src/hardware_matched_job_L8.py`;
`data/heron_counts_matched_L8_<job id>.json` at 50k/30k/16k/4k shots per circuit). The records' `date`
field is the script constant `DATE` of `hardware_matched_job_L8.py` (2026-08-28), not the run date; the
creation times are in `data/sampler_cz_count.json` (`deployed.jobs.<id>.created`).
- Seven circuits prepare a product determinant with X gates (5 up electrons on even qubits, 4 down on odd
  qubits) and apply three layers of RZZ on-site terms and R_xx·R_yy nearest-neighbour hopping (periodic
  wrap), without Jordan–Wigner strings; circuit k evolves for time proportional to k, so k=0 is the bare
  determinant. They are number- and S_z-conserving samplers of the addition sector, not the Hubbard
  propagator of c†|ψ₀⟩.
- `per_time_counts` are the raw counts (7 × shots per job, all bitstrings). The analysis reads them in the
  correct qubit order (`int(bs, 2)`) and post-selects to (N↑, N↓) = (5, 4); no configuration recovery.
  SamplerV2 ran with measurement twirling, Pauli gate twirling and XpXm dynamical decoupling; no
  readout-error mitigation was applied.
- Only the support set S is device-derived. `m0_hat`/`m1_hat` in the records are exact classical values;
  Δ₀, Δ₁ and the relative L1 error are classical functions of S. No moment is estimated from counts.
- The reference distributions come from `delta0_reference_mc.py`, which draws from the exact noiseless
  output distributions of the same seven circuits (qiskit Statevector) and applies the same post-selection
  and arithmetic. Use its percentile ranks, not Gaussian z-scores. The single noiseless Aer seed
  (`aer_noiseless_baseline.json`) sits at the 19.5 / 17.7 / 81.4 / 98.3th percentiles of the raw-matched
  distribution, and `bootstrap_device_curve.py` resamples the retained counts (upward-biased; not a test).

**The 2026-08-21 `L=6` jobs** (`da44t8u1vhnc73fkie30`, U = 4; `da459juaa69c739jcb80`, U = 3.5; 7 circuits ×
35 000 shots each). Their raw REST records are in `data/heron_gatetp_2026-08-21/` (byte for byte; SHA-256
in its README). The submitting notebook read bitstrings in reversed order (`int(bs[::-1], 2)`);
`gate_tp_analysis.py` reproduces its outputs in that order and recomputes them in the correct one. In both
orders the retained support is the whole 300-configuration sector, so a reconstruction on it is exact and a
support-based moment test returns a verdict set by the assumed Hamiltonian, not by the device data. They
enter no result of the paper.

**The companion's `L=6` run** was post-selected in reversed qubit order, so all its retained determinants
came from device errors; its spectra are in `data/heron_spectral.json`, its raw counts are not deposited.

**The `n_k` test.** The on-device `n_k` experiment was pre-registered (manifest v1 sealed 2026-09-02). Its
job was launched twice on 2026-09-03 against `ibm_fez` and submitted nothing: no 12-qubit chain met the
sealed coherence floors. Its single permitted amendment was denied by its sealed gate on 2026-09-05. The
record is `data/2026-09-05_nk_device_attempt_record.json` with the scripts as run in `src/nk_record/`;
`nk_stage2_device_job.py` is the manifest-v1 job and `nk_stage2_device_job_v2_amendment.py` the denied
amendment's job. Re-acquiring device data needs an IBM Quantum account of your own; never commit a token.

## Known gaps

- **`blind_addenda.py`, part (c):** the ARPACK default-start rebuild that first recovered the sealed tie
  order was verified only with scipy 1.13.1 and numpy 1.26.4; it is now a diagnostic (reported under
  `arpack_default_start_rebuild_this_process`), and the rebuild uses the stored order. The shared-state
  thresholds reuse the sealed `var_loc` (centred on the exact E₀), and 4 of the 56 truncations have zero
  reconstructed weight (m̄ ≈ 10⁻³⁰), so their match is absolute, not relative; both are recorded in the JSON.
- **`run_teeth_shared.py` (L=12, `argsort` order):** a rerun reproduces the two 5% crossings (d ≈ 15977 and
  43554, relative change < 2×10⁻⁶) but moves the three smallest-d points (m̄₁ by 1.4%, 0.4% and 6.2% at
  d = 337, 246, 180) against the committed 2026-08-26 record, which is kept and plotted.
- **`fig_sqw` windows:** the display windows clip 6.3% of the broadened charge-grid weight (above 10t; 8.5% at
  the worst q) and 0.25% of the spin weight (above 2.4t); `src/check_sqw_extent.py` →
  `data/2026-09-27_sqw_extent_check.json` measures this and what the raster draws below the onset line
  (displayed value ≤ 0.19, i.e. ≤ 7% of vmax in linear scale, just below the onset; ≤ 0.05 below 5t).
- **Notebook:** `notebooks/00_Reproduce_Everything.ipynb` also runs `run_teeth_shared.py`, which moves
  `teeth_shared.dat` (above); the device reference Monte Carlo is behind a switch that is off by default.
- **`L=24` DMRG check:** at 500 shots the subspace ground level is degenerate, so the point values of that
  row refer to the vector the eigensolver returns; the record gives the range over the ground manifold.
- **Not deposited:** the companion's raw `L=6` counts; the notebook that submitted the 2026-08-21 jobs (its
  pipeline is re-implemented in `gate_tp_analysis.py`, which compares its outputs with the values that
  notebook printed); the 2026-08-24 working copies of the blind scripts (hashes above);
  `src/nk_record/nk_stage2_device_job_launch2.py` is rebuilt from the session log, not the original file.
- **Legacy wording kept:** `data/2026-08-26_teeth_shared_state.json` (not regenerated, see above) still
  says "device-realizable"; `run_collective_sumrule_falsifier.py`'s record names an internal script path;
  `lever1_flip_experiment.py` prints "make-or-break". The `nk_*` gate and seal scripts and their records
  keep their internal workflow identifiers, because the records are sealed and five of the scripts are
  hashed in the manifest; the logged outputs in `src/nk_record/` keep the local paths of the machine that
  ran them.
- **Archive:** the Zenodo DOI is added at submission.
