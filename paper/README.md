# The manuscript, revised on 2026-09-27

This folder holds the manuscript, in the SciPost submission class (`SciPost.cls`, 2021-08 version, the one with
the PhysCore option; `SciPost_bibstyle.bst`). The revtex4-2 build of the same text was 38 pp; the SciPost build
is 63 pp before the length pass. It is not on arXiv or in a journal; the intended venue
is SciPost Physics Core. The 2026-09-05 text (31 pp) was revised on 2026-09-27 in three steps: the
recomputations (phase A; scripts and dated outputs in [`../docs/REPRODUCE.md`](../docs/REPRODUCE.md)), every
section rewritten against them (phase B), and a paper-wide consistency pass over the text, the figure
sources, the code and the bibliography, each change checked by a second reader (phase C). Captions and
figures now agree, and `main.pdf` is built from this text.

**Not done yet (before submission):**
- a length pass (the introduction is about 2.4k words; several captions exceed 120 words) and the port to
  the SciPost template, after which every figure needs a new visual check at single-column width;
- the companion (arXiv:2608.16436) is cited with its v3 title and v3-only content, and v3 is announced on
  2026-09-28; if it is not public at submission, the title, the version pin and the two sentences that
  rely on it (Introduction; Sec. role) must be softened;
- the Zenodo DOI in the Data availability statement.

## Defects of the 2026-09-05 text, corrected in the revision

Each item describes the 2026-09-05 text. The revised `main.tex` was checked for every item on 2026-09-27.

**Blinded test** (`data/2026-08-24_blind_harness_score.json`, `data/2026-09-27_blind_addenda.json`)
- The sealed primary endpoint is not reported: TPR 101/188 = 53.7% [46.6, 60.7] at FPR 1/112 = 0.9%
  [0.2, 4.9].
- Only one of the two failed sealed secondary predictions is reported as failed; the m₀-preserving
  spurious-feature subclass (sealed as designed-blind, caught 20/34) is presented as a success.
- The 0.9% FPR is set by construction: the simulated estimator is unbiased while the threshold budgets a
  2% bias (modelled FPR 1.3% unbiased, 8.0% with +2% bias).
- The Hankel branch is evaluated on the reconstruction's own moments and cannot fire (0/300 by
  construction).
- One-node Krylov passes only inside the 2% bias budget; Krylov with n_l ≥ 2 matches m₀–m₂ exactly.
- The text says the seal means no threshold could be tuned after seeing an outcome. The seal is
  self-attested (no third-party timestamp) and cannot by itself exclude pre-seal iteration. The calibrated
  scope is not stated: L=6 doped ring, U/t ∈ {3, 4, 5, 6}; truncations keep 75.1–98.7% of the
  ground-state weight; spurious atoms are point atoms with 2.1–24.9% of m₀.

**"Same samples" and hardware** (`data/2026-09-27_delta0_reference_mc.json`, `data/heron_counts_matched_L8_*.json`)
- The abstract and introduction say each moment is estimated from the same computational-basis samples
  and that one moment is evaluated on hardware. Neither holds: same-shot estimation was done only in
  simulation (ρ ≈ 0.43, no power gain), and the device supplies only the sampled support
  (`m0_hat`/`m1_hat` are exact classical values).
- The device circuits are described as a Trotterized evolution of c†|ψ₀⟩ with configuration recovery. They
  prepare a product determinant (X gates, 5 up, 4 down) and apply R_xx·R_yy hopping without Jordan–Wigner
  strings; the support is defined by post-selection alone.
- The text names a readout-error-mitigation method that was not used: SamplerV2 ran with measurement
  twirling, Pauli gate twirling and dynamical decoupling; no readout-error mitigation was applied.
- "Shot noise is negligible" and the single-seed comparison are wrong. Against 2000 noiseless replicas at
  the device's raw shots, the device Δ₀ lies below all at 50k, at the 1.7th percentile at 30k, at the
  57.9th at 16k and above all at 4k; the analytic σ understated the spread 34–195-fold.
- Δ₀ rewards support spreading: a uniform in-sector sampler at the device's retained counts has lower Δ₀
  than the device at every budget. The text calls the 50k run a "clean pass"; its A(ω) has relative L1
  error 0.35 (and any pass threshold at L=8 would apply the L=6 sealed rule counterfactually, with an exact
  m̂₁).
- The retention-matched comparison quotes one noiseless draw (0.2032 at 50k). That draw sits at the 95.4th
  percentile of its distribution, and the comparator depends on how the kept shots are allocated over
  circuits (mean 0.150 with an equal split, 0.224 per circuit, at 50k).
- The circuit-figure caption says "the deployed run is L=6, 12 qubits"; the bond-basis estimator it draws
  was simulated only and never run on hardware.
- The companion's L=6 run was post-selected in reversed bit order, so all its retained determinants came
  from device errors.
- The `fig_device` caption (readout mitigation, one Aer seed, analytic σ) no longer matches the two-panel
  figure.

**`n_k` experiment** — the gate-negative outcome (amendment denied on 2026-09-05; not run) is not
reported, and `n_k` is still called the next experiment (`data/2026-09-05_nk_v4_verdict.json`).

**Theory and limits** (`data/2026-09-27_theory_numerics.json`)
- "3.8× power margin" and "conditioning is benign" overstate the power. Guaranteed power at the deployed
  order is zero because of the null-space dimension: a redistribution that keeps m₀–m₂ exactly moves 72% of
  the weight, and one that keeps m₀, m₁ with |Δ₂| = τ₂ moves 86%; both pass. κ of the node Vandermonde is
  50 on rescaled nodes; the raw-node 1.2×10⁶ is unit-dependent.
- "Necessary ∘ sufficient" is wrong: the Hankel and Hausdorff bounds are both necessary conditions, and the
  Wang/Mortimer composition is not implemented. `m₂ ≤ 573.7` used [0.148, 52.21]t; the weighted support
  gives 147.6 (oracle), the a-priori supports 579.8 or 243.3.
- The Christoffel caption quotes grid maxima 0.9998 and 0.875; the exact values are 1 (at the mean) and
  0.882. The "~20–30 moments" statement has no source. A passing reconstruction is bracketed only by the
  tolerance-inflated atom: 3.2–7.7× W₁ at t = 14–20 and 10–11× at the figure's t* on the figure's U/t=8
  measure (3.6–5.7× at t = 12–20 on the U/t=4 measure).
- `m₋₁` is not the sharpest moment under the deployed error model: m₁ and m₂ fire at smaller injected
  fractions on both the ring and the open chain.
- The estimator forms differ off an eigenstate (4.558 vs 3.717 at 44.4% coverage), enough to change a
  verdict; the deployed form must be fixed.
- The operator O = c(H−E₀)c† is called positive semidefinite; it is not (lowest eigenvalue −0.818). The
  bond-estimator probe shift ⟨ρ_q⟩ (checked at q = π) vanishes by reflection symmetry, not by half filling
  (that ground state has N = 4 on 6 sites).

**Shot-budget interval test** (`data/2026-09-27_interval_battery.json`, `data/2026-09-27_interval_moment_mc.json`)
- The "single honest exception at d=98" is a tie-break artefact (m₁ alone rejects on 57% of 300 tie
  choices). The tie-robust lone-m₁ miss is d=97. At 2% bias the whole battery passes 12 truncated
  reconstructions, three at 46–51% coverage, so the battery does not close the blind spot. Panel (b) of
  `fig_momentmc` now plots every depth.

**Figure captions**
- `fig_separating`: 26 instances lie in the shaded region, not 24 (m₁ on 16, m₂ alone on 10), including
  truncation id 162 and the clean false positive id 184; the exemplar id 222 is weight-preserving; the
  thresholds include the 2% bias term (`data/2026-09-27_separating_counts.json`).
- `fig_gausslaw`: ε is an amplitude, not a weight fraction (the unphysical weight is 0.155 at ε = 0.3);
  the admixture is a minimal ΣG² = 2 violation (`data/2026-09-27_gauss_state.json`).
- `fig_sqw`: the rasters are now placed at their true extents; the markers are filled; the charge window
  is clipped at 10t.

**Framing and citations** — the wording about corroboration and knowledge; the overlap with the
companion's v3 (captured weight, moment exactness through order 2K+1, the leakage identity, "weight alone
cannot certify"), which lacked citations and a statement of what this work adds; several bibliography
entries.

## Build

`make paper` compiles `main.tex` with the committed `main.bbl` (no BibTeX); after editing `refs.bib`, run
`pdflatex main && bibtex main && pdflatex main && pdflatex main` to refresh it. The source of
`figs/fig_circuit_qtk.pdf` is `figs/src/fig_circuit_qtk.tex` (`pdflatex fig_circuit_qtk.tex`). Files no
figure reads, and the figure files replaced on 2026-09-27, are in `figs/_superseded/`.
