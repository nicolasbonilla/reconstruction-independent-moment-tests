# Process record

This file keeps the design and analysis history behind the paper. The paper states results and the true account of
each test; how the sealed designs, records and figures reached their final form is recorded here. Every statement
points to a committed file. Licence: CC BY 4.0, like the other documentation in `docs/`.

Paper pointer: `paper/main.tex`, Data availability and Appendix "Cost of a device-native estimator" (sealed `n_k` test).

---

## 1. The sealed on-device `n_k` test: Stage-0 design gate (2026-09-02)

The pre-registered `n_k` test (doped L=6 ring, U/t=8, 12 qubits) was sealed on 2026-09-02 (manifest
`data/manifest_nk_device_v1.json`, SHA-256 `4962...ae80`). Before the seal, a simulated Stage-0 gate (a noisy
simulation of the sealed circuits and estimators, run before any device time is spent) checked the design. Stage-0 is
a design check, not itself a pre-registered test; the gate records call its outcome a design result. It reached a pass
only at its fourth step, after two revisions of the noise model and a replication, all on the same day and before the
seal:

| step | record | noise model / change | verdict |
|---|---|---|---|
| v1 | `data/2026-09-02_nk_stage0_gate.json` | first model | PROCEED; an internal audit then ruled it a model artifact (relaxation omitted, a finite-twirl residue and a transpilation bug) |
| v2 | `data/2026-09-02_nk_stage0_gate_v2.json` | relaxation and finite (32-instance) gate twirling added | STOP on 3 of its 7 pessimistic rows |
| v3 | `data/2026-09-02_nk_stage0_gate_v3.json` | after an adjudication: heavy-hex instead of line layout; relaxation scheduled at the floors T1 >= 150 us, T2 >= 100 us (which first appear in this version and become the manifest's day-of floors), in place of the exposure multipliers behind two v2 STOP rows; a pulse-identical twin circuit for the post-mitigation residual check | STOP on 1 of 9 rows (eps = 3.5e-3, twirled, second seed, at the floors): residual check 0.026674 against the threshold 0.026667 |
| v3R2 | `data/2026-09-02_nk_stage0_v3R2.json` | replication of that row over ten seeds (the two v3 seeds included) under a rule fixed before any replica ran: pass if mean + 1.645 standard errors < threshold | PASS: mean 0.0239, upper confidence bound 0.0260 < 0.0267; two of the ten replicas individually exceeded the threshold (0.026674 and 0.0293) |

The manifest was sealed on the v3R2 pass. At the manifest floors the residual check ran at 74-100% of its threshold.
Two candidate calibration estimators were tried and rejected before the seal (both would fail the binding simulated
configuration). Here `eps` is a compound harness knob that moves the depolarizing and coherent-ZZ components together;
it is not a physical two-qubit error rate (the paper's forecasts use 2.5e-3).

## 2. The `n_k` amendment and its gate (2026-09-04 / 2026-09-05)

The single permitted amendment would have relaxed the floors to T1 >= 100 us and T2 >= 70 us, added a per-qubit
readout cap and widened the backends to `ibm_fez`, `ibm_marrakesh` and `ibm_kingston`, provided that three
pre-committed simulated verdict rows all passed (`data/2026-09-04_nk_v4_row_R1ro0.018.json`, `_R2`, `_R3`). The
worst-admissible row (T1 = 100 us, T2 = 70 us, eps = 3.5e-3) failed on 2026-09-05 at the tightest per-qubit readout
cap of the pre-committed ladder (1.8%; the looser rungs, 3.0% and 2.5%, are partial rows taken from run logs, with every
recorded replica above the threshold): upper confidence bound of the post-mitigation residual check 0.0287 against the
sealed threshold 0.0267. A row identical except eps = 2.85e-3 passed (0.0258); since eps is a compound knob, the pair
of rows does not locate a physical two-qubit error rate. The amendment was denied
(`data/2026-09-05_nk_v4_verdict.json`). Why an amendment was prepared rather than the sealed manifest run is now
documented by the day-of chain checks (paper, Appendix "Cost of a device-native estimator";
`data/2026-09-05_nk_device_attempt_record.json`): on 2026-09-03 no 12-qubit chain of `ibm_fez` met the sealed floors.

## 3. Re-seal of the amendment verdict file (2026-09-05)

The verdict file `data/2026-09-05_nk_v4_verdict.json` was corrected and re-sealed on the day of the verdict (current
SHA-256 `7849...325a`, in `data/2026-09-05_nk_v4_verdict.json.sha256`). Its first digest (`f1b7...be3d`) failed
verification because it had been computed over LF line endings while the file was written with CRLF. In the same
correction its `finding` text, which had asserted that tightening the readout cap did not move the residual check
(the ladder in the same file shows a 14.2% reduction against the 18.4% needed), was rewritten by hand, and caveats on
attributing the failure to a single noise axis were added. The verdict, amendment denied, is unchanged; the file's
`corrections` record lists only the text and the digest.

## 4. Superseded device noise model and the bootstrap

A simpler model of the L=8 coverage residual drew shots from |phi|^2 at the pooled retained count
(generator `src/analytic_coverage_leak.py`, output `data/_superseded/analytic_coverage_leak.json`). It underestimates
the spread of Delta_0 by factors of 34 to 195 relative to the raw-matched noiseless replicas
(`data/2026-09-27_delta0_reference_mc.json`, `budgets.<N>.sd_over_analytic_sigma`), and under its own threshold a sampler drawing exactly from |phi|^2 would exceed
that threshold at every budget (modeled false-positive rate 100%). Neither that model nor any z-score derived from it
is used. A nonparametric bootstrap of the retained counts (`src/bootstrap_device_curve.py`,
`data/bootstrap_device_curve.json`) is biased upward (its 50k interval [0.015, 0.101] excludes the point estimate
0.0044) and is not used as a test.

## 5. The estimator form O_k: chronology

In the committed records the operator form O_k = J^dagger (H - e_0)^k J is first used off an eigenstate, with the
leak-identity justification, in the shared-state test (`src/run_teeth_shared.py`,
`data/2026-08-26_teeth_shared_state.json`), after the sealed blinded test (which used an oracle-centered estimate and
does not depend on the form) and before the off-eigenstate example of the d=100 cut
(`src/estimator_form_offeigenstate.py`, key `R6_estimator_forms_off_eigenstate` of
`data/2026-09-27_theory_numerics.json`) was computed. That example is evidence for the choice, not its origin.

## 6. Tie-order history of the truncation scans

The L=6 ground-state probabilities come in exactly degenerate shells; which members of a shell a truncation keeps is a
convention. Earlier scans broke ties by floating-point noise or by NumPy `argsort` (which does not fix the order of
ties); since 2026-09-27 every scan uses the deterministic `lexsort` order (descending |psi_0|^2, ties by ascending
sector index).

- d=97 catch: `argsort` gave |dm1| = 0.155, |dm2| = 2.560; `lexsort` (used by `src/verify.py` and
  `src/interval_battery.py`) gives 0.186 / 2.729. The paper prints the `lexsort` values.
- d=98: an earlier run with ties broken by floating-point noise found |mhat_1 - mbar_1| = 0.225 (a miss); the
  deterministic order gives 0.284 (a rejection); 300 random shell choices give 0.172-0.295.
- Monte Carlo power: an earlier grid-only analysis (`data/2026-08-24_interval_moment_mc.json`) showed a single
  low-power point near 53% coverage; the scan at every d (`data/2026-09-27_interval_moment_mc.json`) supersedes it.

## 7. The Gauss-law figure: first admixture rule

The Gauss-law figure was first made with the unphysical eigenstate of H of largest overlap |<k|M|psi_0>|, with
M = sum_n (-1)^n n_n the staggered occupation. That rule is ill-posed: M is gauge-invariant, so every such overlap
vanishes (the largest is 1.8e-15), and the pick, set by rounding noise, happened to be a state with g_u = 2. The
explicit rule now used (`src/export_gauss_state.py`, `data/2026-09-27_gauss_state.json`) gives the same plotted curve.

## 8. Items the records do not document

- The choice of the U/t=8 current measure for the Christoffel figures (the blinded catalog is U/t in {3,4,5,6}); the
  paper reports the in-catalog U/t=4 check alongside.
- The joint-sampling study (`data/2026-09-01_lever1_oracle.json`, seal `data/2026-09-01_lever1_prereg_seal.json`)
  records the oracle stage only, as point estimates without confidence intervals, and records no verdict on its six
  sealed negative rules.

## 9. Paper restructure map (2026-10-02)

The 2026-10-02 restructure moved, and did not delete, results. Former location -> new location in `paper/main.tex`:

| content | before | after |
|---|---|---|
| leak identity, Eq. (app-leak) | App. B | Proposition 1, Sec. 2.3 (proof in App. B) |
| O_k vs commutator flip (d=100) | App. A | Sec. 2.3 (details App. A) |
| Fig. teeth (idealized + circular control) | Sec. 4 | merged into Fig. teeth_shared (Sec. 3.1) |
| Figs. akw, sqw, collective_screen | Sec. 3 | App. D (Illustrations) |
| Fig. momentcone, SDP layer (Sec. 9.4) | Sec. 2.2 / 9.4 | App. E |
| Theorem 1, Prop. (Markov-Krein), Fig. bracketing/christoffel | Sec. 7 | App. F (summary in Sec. 5.1) |
| blinded deviations (i)-(iv), id 162/184, legend counts | Sec. 5 | App. C |
| bond-basis estimator, QWC counts, Fig. circuit, noise forecast, n_k channel, sealed n_k paragraph | Sec. 6 / App. B / App. F | App. J |
| companion L=6 run (former Fig. heron) and the 2026-08-21 jobs | Sec. 6 | one paragraph in Sec. 6; details App. I (Fig. heron no longer input; `paper/figs/fig_heron.tex`, `data/heron_spectral.json` kept) |
| L=24 scalability check | Sec. 9.3 | App. K |
| "deployed" (rule, order, battery, thresholds) | throughout | "default" |

## 10. Earlier versions of the manuscript

Each version is in git history (commit in brackets); `_superseded/paper_2026-09-02/` keeps an upload
bundle of the 2026-09-02 draft, which was not submitted.

| date | version | commit |
|---|---|---|
| 2026-09-05 | text of 31 pp | `c4f0d77` |
| 2026-09-27 | recomputations as committed scripts with dated outputs (phase A), every section rewritten against them (phase B), paper-wide consistency pass over text, figure sources, code wording and bibliography (phase C) | `118e495`, `9a296bc`, `422e116` |
| 2026-09-28 | port to the SciPost Physics Core class (63 pp); length pass to 56 pp with no result removed; referee and visual review applied (65 pp) | `d1dab52`, `b2cf04d`, `70a3a6a` |
| 2026-10-02 | restructure (section 9), retitle to "... for sample-based spectral functions", the 2026-10-02 recomputations (`docs/REPRODUCE.md`, Dated recomputations) | working tree |

Statements of the 2026-09-05 text that the 2026-09-27 revision corrected (each checked against the data
files named):

- **Blinded test** (`data/2026-08-24_blind_harness_score.json`, `data/2026-09-27_blind_addenda.json`). The
  sealed primary endpoint (TPR 101/188 = 53.7% [46.6, 60.7] at FPR 1/112) was not reported; only one of the
  two failed sealed secondary predictions was reported as failed (the m₀-preserving subclass, sealed as
  designed-blind and caught 20/34, was presented as a success); the low FPR was not attributed to its cause
  (unbiased simulated estimator, threshold budgeting a 2% bias); the Hankel branch, evaluated on the
  reconstruction's own moments, cannot fire (0/300); one-node Krylov passes only because τ₂ includes the
  2% bias term, and Krylov with n_l ≥ 2 matches m₀–m₂ exactly; the seal was said to exclude any tuning after
  an outcome, but it is self-attested; the calibrated scope was not stated.
- **"Same samples" and hardware** (`data/2026-09-27_delta0_reference_mc.json`,
  `data/heron_counts_matched_L8_*.json`). The text said each moment was estimated from the same
  computational-basis samples and one moment on hardware; same-shot estimation was done only in simulation
  and the device supplies only the sampled support. The device circuits were described as a Trotterized
  evolution of c†|ψ₀⟩ with configuration recovery and readout-error mitigation; they prepare a product
  determinant and apply R_xx·R_yy hopping without Jordan–Wigner strings, the support is defined by
  post-selection alone, and no readout-error mitigation was applied. "Shot noise is negligible" (from the
  analytic model, now in `data/_superseded/`) and the single-seed comparison were replaced by 2000-replica
  reference distributions; Δ₀ rewards support spreading (a uniform in-sector sampler beats the device); the
  "clean pass" at 50k has relative L1 error 0.35; the retention-matched control quoted one draw at the 95.4th
  percentile of its distribution; the circuit-figure caption called the simulated bond-basis estimator "the
  deployed run"; the companion's `L=6` run was post-selected in reversed bit order.
- **`n_k` experiment.** The gate-negative outcome (amendment denied on 2026-09-05; not run) was not reported,
  and `n_k` was called the next experiment.
- **Theory and limits** (`data/2026-09-27_theory_numerics.json`). "3.8× power margin" and "conditioning is
  benign" overstated the power (zero guaranteed power at m₀–m₂; κ = 50 on rescaled nodes); "necessary ∘
  sufficient" was wrong (Hankel and Hausdorff are both necessary; `m₂ ≤ 573.7` used [0.148, 52.21]t); the
  Christoffel caption quoted grid maxima; `m₋₁` was called the sharpest moment; the estimator forms differ
  off an eigenstate; O = c(H−E₀)c† was called positive semidefinite (λ_min = −0.818); ⟨ρ_π⟩ = 0 was
  attributed to half filling instead of reflection symmetry.
- **Shot-budget interval test** (`data/2026-09-27_interval_battery.json`). The "single exception at d=98"
  is a tie-break artefact; the tie-robust lone-m₁ miss is d=97, and at 2% bias the whole battery passes 12
  truncated reconstructions.
- **Figure captions and framing.** `fig_separating` counted 24 instances in the shaded region, not 26;
  `fig_gausslaw` treated ε as a weight; `fig_sqw` rasters were not at their true extents; the overlap with
  the companion lacked citations and a statement of what this work adds.

The 2026-10-02 recomputations corrected, in turn, statements of the 2026-09-28 text: the noise forecast's
reference depth (292 CZ, one deployed circuit count, was called the companion's sampler depth; the median
over the deployed `L=8` circuits is 295, and the companion's 12-qubit circuits need 213;
`data/sampler_cz_count.json`); the date of the `L=8` jobs (created 2026-08-29 UTC, not 2026-08-28); the
`L=24` starting value Δ₀ = 44.2 (a degenerate subspace ground level, range 40.1–44.5;
`data/2026-10-02_beyond_ed_dmrg_catch.json`); the uniform-noise mixture called monotone at 50k (it has a
minimum near p ≈ 0.93), lowest near p = 0.75–0.81 at 30k and 16k (fine-grid minima p ≈ 0.88 and 0.77),
and increasing with p at 4k (0.270 → 0.326; it is flat to 0.3% up to p = 0.25, with a minimum 0.269 near
p ≈ 0.09) (`exact_E_fine_grid` of `data/2026-09-27_delta0_reference_mc.json`); and the reason the sealed `n_k` job was
not run, now recorded from the logs and re-queried calibrations
(`data/2026-09-05_nk_device_attempt_record.json`). The false-positive rate at the noise scale of the sealed
text (40/112) and the full-support Rayleigh–Ritz pass (relative L1 0.93) were computed on 2026-10-02 and
are reported in the paper.
