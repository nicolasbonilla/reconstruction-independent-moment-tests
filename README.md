<div align="center">

# Reconstruction-independent moment tests
### *for quantum-computed dynamical spectra — a necessary-condition screen*

**Nicolás Bonilla Vargas** &nbsp;[![ORCID](https://img.shields.io/badge/ORCID-0009--0006--6155--4391-A6CE39?logo=orcid&logoColor=white)](https://orcid.org/0009-0006-6155-4391)

[![Manuscript](https://img.shields.io/badge/manuscript-2026--09--05%20text%2C%20revision%20in%20progress-lightgrey.svg)](paper/README.md)
[![License: MIT](https://img.shields.io/badge/code-MIT-green.svg)](LICENSE)
[![License: CC BY 4.0](https://img.shields.io/badge/paper-CC--BY--4.0-lightgrey.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](requirements.txt)

</div>

> **Status (2026-09-27).** The manuscript is under revision. It is **not on arXiv** and not in a journal.
> On 2026-09-27 the numbers the revision needs were recomputed as committed scripts with dated outputs
> (`data/2026-09-27_*.json`; "phase A", listed in [`docs/REPRODUCE.md`](docs/REPRODUCE.md)), and the text
> in [`paper/`](paper/) (now 38 pp) was revised against them: several statements of the 2026-09-05 text are
> corrected or withdrawn (listed in [`paper/README.md`](paper/README.md)). The length pass and the port to
> the SciPost template are still to come. Where text and data disagree, the code and data in this
> repository are the reference.

> **One sentence.** A reconstruction-independent, **necessary-condition** screen for sample-based quantum
> spectral functions: low-order spectral moments m₀, m₁, m₂ (and m₋₁ on open chains) are evaluated as
> ground-state operator expectations, never read back off the reconstructed A(ω), and compared with the
> moments of A(ω). A pass is consistent with the reconstruction; it is not proof of it.

This repository holds the code and data behind the paper. Numbers and figures trace to a script and a
data file listed in [`docs/REPRODUCE.md`](docs/REPRODUCE.md) and [`docs/FILE_INDEX.md`](docs/FILE_INDEX.md);
what is still missing is listed under [Known gaps](#known-gaps). Most non-hardware results regenerate from
source with `numpy`/`scipy`/`mpmath`; the device re-analysis and the noise-model forecasts also need
`qiskit`/`qiskit-aer`, and the `L=24` DMRG check needs `physics-tenpy`.

> **Companion works.** Data and pipeline — [arXiv:2608.16436](https://arxiv.org/abs/2608.16436); its v3 is
> titled *Captured weight and boundary leakage bound the error of sample-based spectral functions*
> (v1–v2: *Dynamical spectral functions from bitstring-sampled quantum subspaces*). Review — *Machine
> learning for sample-based quantum diagonalization: a review of generative configuration recovery and the
> classical-simulability frontier* ([arXiv:2608.05314](https://arxiv.org/abs/2608.05314)).

**Relation to the companion.** Several ingredients used here are results of the companion's v3: the
captured weight w (the coverage residual of this work is Δ₀ = m₀(1 − w)); exactness of the moments
through order 2K+1 for subspaces that contain the Krylov space (which is why Krylov reconstructions with
n_l ≥ 2 Lanczos vectors match m₀–m₂ exactly and pass by construction); the leakage identity
μ₂ − μ₂^S = ‖Q_S H P_S φ‖² (the operator leak behind the shared-state estimator); and the statement that
captured weight alone cannot certify a reconstruction. What this work adds: the moment battery used as a
screen whose estimates do not come from A(ω); a hash-sealed, blinded calibration of it, with the addenda
that qualify that calibration; the exact limits of what a pass means at the deployed order; a reference
Monte Carlo for the coverage residual on its own `L=8` `ibm_fez` counts; and a pre-registered `n_k` device
experiment that was gate-negative.

Figures: see `paper/`; README thumbnails will be regenerated from the native PGFPlots sources after the
revision.

---

## What the work shows

1. **A reconstruction-independent screen.** A reconstructed spectrum `A(ω)` is checked against
   **necessary conditions**: static ground-state spectral moments `mₖ = ∫ωᵏA(ω)dω = ⟨Oₖ⟩`, with the
   deployed first-moment form `⟨J(H−e₀)J⟩`, evaluated **independently** of the reconstruction and never read
   back off `A(ω)`. In this work that independent evaluation is done in simulation, or classically on the
   sampled subspace state; it is never done on hardware. The estimator forms `⟨J(H−e₀)J⟩` and
   `½⟨[J,[H,J]]⟩` agree only on an exact eigenstate: on a truncated state (L=6, 44.4% coverage) they give
   4.558 and 3.717, enough to change a verdict.

2. **Independence matters.** On subspace truncation the independent estimator fires while a back-computed
   (circular) control stays blind — an exact classical simulation of the doped `L=12` Hubbard ring, not a
   hardware result.

3. **Several error classes and domains.** The same construction catches a spurious low-frequency peak that
   preserves the total weight `m₀` but shifts `m₁`; a separate state-level check catches a
   gauge-symmetry-breaking admixture through the U(1) Gauss-law identity `⟨∑G²⟩ = 0`. The negative moment
   `m₋₁` gives no detection advantage under the deployed error model once noise is accounted for (see
   Outcomes).

4. **What a pass means.** In exact arithmetic the classical Christoffel/Markov–Krein bounds limit the weight
   a moment-matching measure can misplace, and the blind spot closes at order N−1. At the deployed order
   (`m₀`–`m₂`) the guaranteed power against fixed-pole redistributions is **zero**, because the moment rows
   leave a null space (dimension 2 for the 5-pole L=6 measure): a redistribution that keeps `m₀`, `m₁`,
   `m₂` exactly moves 72% of the weight, and one that keeps `m₀`, `m₁` with `|Δ₂| = τ₂` moves 86%; both pass.

5. **A shot-budget form, with its misses.** The exact-moment idealization is replaced by a shot-budget
   interval test at the `ibm_fez` 50k-shot budget. The joint `(m₀,m₁,m₂)`+Hankel battery catches some
   lone-`m₁` misses (the tie-robust case is `d=97`), but a scan of every truncation depth shows that at the
   deployed 2% bias it still passes 12 truncated reconstructions (see Outcomes).

**Scope.** The screen **rejects; it does not certify**. Its conditions are **necessary, not sufficient**.
It targets one error class (subspace truncation and weight misassignment), not device noise in general.
**No quantum advantage is claimed.** On IBM Heron the device supplies only the sampled support set; no
spectral moment is estimated on hardware, and the discriminating off-diagonal moments are simulation-only.

---

## Outcomes, including the negative ones

**Blinded test (hash-sealed simulation, 300 instances).**
- **Sealed primary endpoint:** TPR = 101/188 = **53.7%** (Wilson 95% [46.6, 60.7]%) at FPR = 1/112 = 0.9%
  ([0.2, 4.9]%) (`data/2026-08-24_blind_harness_score.json`; seal `data/prereg.sha256`). Per class:
  truncation 56/56, spurious features 45/61, under-converged Krylov **0/71**.
- **Two sealed secondary predictions failed:** one-node Krylov, expected to be caught, was caught 0/15; the
  m₀-preserving spurious-feature subclass, sealed as designed-blind, was caught 20/34 (59%, Wilson
  [42, 74]%). The register was not revised.
- **The low FPR is set by construction.** The simulated estimator was the exact moment plus unbiased
  Gaussian noise, while the frozen threshold budgets a 2% bias. Redrawing the 112 clean instances under the
  frozen rule gives a modelled FPR of 1.3% unbiased, 8.0% with a realized +2% bias and 11.5% with −2%
  (`data/2026-09-27_blind_addenda.json`). The 0.9% checks the implementation, not robustness to mitigation
  bias. The nominal joint target was 1−0.95³ = 14.3%.
- **Scope of the calibration:** L=6 doped ring, current probe, U/t ∈ {3, 4, 5, 6}; truncations kept 15.6–74.2%
  of the ground-state support (75.1–98.7% of its weight); spurious atoms carried 2.1–24.9% of m₀ and are
  point atoms in the generator (the prereg text says "Gaussian atom").
- **Post hoc addenda** (`data/2026-09-27_blind_addenda.json`): the 56 sealed truncations were rebuilt exactly
  (from a stored tie order that every run re-validates against the sealed record; the script aborts
  otherwise) and rescored post hoc with the shared-state estimator; all 56 are rejected. One-node Krylov misses m₂ by 2.35–4.19%
  and passes only inside the 2% bias budget (largest residual 0.986 τ₂). The Hankel branch of the rule is
  evaluated on the reconstruction's own moments, so it cannot fire (0/300 by construction).
- **The seal is self-attested** (no third-party timestamp); it fixes the protocol but cannot by itself
  exclude pre-seal iteration.

**"Same samples" was realized only in simulation.** `joint_covariance_composition.py` and
`lever1_flip_experiment.py` draw one set of shots that gives both the independent estimate and the
reconstruction: correlation ρ ≈ 0.43 at 84% coverage and no power gain at matched FPR
(`data/2026-09-01_lever1_oracle.json`). Their local estimator uses the exact ground-state amplitudes, so it
is not device-realizable as written. The `L=24` DMRG check shares the sampled subspace. On hardware it was
never done: the `m0_hat`/`m1_hat` fields of the device records are exact classical values.

**Hardware (`ibm_fez`, L=8, jobs of 2026-08-29): the device supplies only the sampled support.**
- **What ran.** Seven circuits (`src/hardware_matched_job_L8.py`) prepare a product determinant with X gates
  (5 up, 4 down electrons) and apply R_xx·R_yy hopping with RZZ on-site layers, without Jordan–Wigner
  strings; they sample the addition sector and are not the Hubbard propagator of c†|ψ₀⟩. Bitstrings are
  read in the correct order and post-selected to (N↑, N↓) = (5, 4), with no configuration recovery. SamplerV2
  ran with measurement twirling, Pauli gate twirling and dynamical decoupling; no readout-error mitigation
  was applied.
- **Coverage residual.** Δ₀ = 0.0044 / 0.032 / 0.123 / 0.373 at 50k/30k/16k/4k shots per circuit, as coverage
  of the 3920-configuration sector falls from 94% to 35%.
- **Reference Monte Carlo** (`src/delta0_reference_mc.py`, `data/2026-09-27_delta0_reference_mc.json`;
  2000 multinomial replicas of the exact noiseless output distributions at the device's raw shots). The
  device Δ₀ lies below all 2000 replicas at 50k, at the 1.7th percentile at 30k, at the 57.9th at 16k and
  above all 2000 at 4k. Device noise widens the support at high budgets and narrows it at low ones: the
  device |S| is 3682 against 3379 ± 14 noiseless at 50k and 1379 against 2027 ± 21 at 4k, outside all 2000
  replicas at every budget. The superseded analytic shot-noise model understated the spread 34–195-fold
  against this comparator (`data/_superseded/analytic_coverage_leak.json`); no z-scores are quoted.
- **Δ₀ rewards support spreading.** A uniform in-sector sampler at the device's retained counts gives
  E[Δ₀] = 3.4×10⁻¹⁰ / 1.5×10⁻⁶ / 5.2×10⁻⁴ / 0.093, and its reconstruction has relative L1 error below
  4×10⁻⁹ (the Lanczos floor) at 50k and 30k, 0.014 at 16k and 0.90 at 4k, against the device's
  0.35 / 0.56 / 0.90 / 1.12. It beats the device on both measures at every budget (at 4k, 0.90 ± 0.08
  over 10 replicas against the device's 1.12). Δ₀ shows that the diagonal channel runs on device data; it
  cannot tell a good device from noise, and it does not certify A(ω).
- **Necessary is not sufficient on device data.** At 50k both residuals are small (Δ₀ = 0.0044,
  Δ₁ = 0.019) while the reconstructed A(ω) has relative L1 error 0.35 (`falsifier.rel_L1`). The device's
  0.35 is below all 80 noiseless replicas at equal raw shots (0.566 ± 0.032); at 16k and 4k the device line
  shape is worse than noiseless (96th and 100th percentile).
- **Retention-matched control.** The device keeps 23.5–24.0% of its shots after post-selection; at 50k it
  keeps 92.4% of circuit k=0 (a single determinant when noiseless) against 11.3–13.1% of the others, so the
  comparator depends on how the kept shots are allocated. At the device's kept count per circuit the
  noiseless Δ₀ is 0.224 / 0.279 / 0.343 / 0.443 (device below all replicas at 50k–16k, 0.6th percentile at
  4k); with an equal split over circuits it is 0.150 / 0.205 / 0.273 / 0.404. The committed single-draw
  control (0.203 / 0.248 / 0.329 / 0.418, `data/2026-09-05_retention_matched_control.json`) sits at the
  95.4 / 93.6 / 99.65 / 78.9th percentiles of the equal-split distribution, so its 46× ratio at 50k is a
  high draw; the MC ratios are 34× (equal split) and 51× (per circuit) at 50k, and 1.08–1.19× at 4k.
- **The companion's L=6 run** (`fig_heron`) was post-selected in reversed bit order, so all its retained
  determinants came from device errors.

**`n_k` device experiment: gate-negative, not run.** The Stage-0 gate passed and manifest v1 was sealed on
2026-09-02 (`data/manifest_nk_device_v1.sha256`); its day-of rule requires a 12-qubit chain with
T1 ≥ 150 µs and T2 ≥ 100 µs on every qubit. The single permitted amendment (floors lowered to 100/70 µs, with
readout and CZ caps) failed its sealed noisy-simulation gate on 2026-09-05 (AMENDMENT-DENIED; row
R1ro0.018, ε = 3.5×10⁻³: canary UCB 0.02867 > sealed threshold 0.02667), so manifest v1 stands and has not
been run. A row identical except ε = 2.85×10⁻³ passed (R2); ε is a compound noise knob, not a single
physical parameter (`data/2026-09-05_nk_v4_verdict.json`).

**Shot-budget interval test** (L=6, U/t=4, 50k shots; deterministic truncation order since 2026-09-27;
`data/2026-09-27_interval_{moment_closure,battery,moment_mc}.json`).
- The "single exception at d=98" of the 2026-09-05 text does not survive. The cut at d=98 falls inside a
  24-fold degenerate |ψ₀|² shell; over 300 random choices of the kept shell members, `m₁` alone rejects on 57%.
- `d=97` is the tie-robust lone-`m₁` miss that `m₂` catches (`m₁` rejects on 0%, the battery on 100% of 300
  tie choices).
- Scanning every depth d = 1…185 at 2% bias: `m₁` alone misses 17 truncations; `m₂` catches 3 of them; the
  whole `(m₀,m₁,m₂)`+Hankel battery **passes 12 truncated reconstructions**, three of them at 46–51%
  coverage (d = 86, 87, 95). Over 300 tie choices the battery still passes d = 86 on 98.7% of them and
  d = 161–165 and 183–185 on at least 99.7%; the passes at d = 87, 95 and 160 depend on the tie order.
- Monte Carlo power of the lone `m₁` test (shot-only threshold) falls to the false-positive level in three
  narrow windows (lowest: 5.0% at d=162, 5.5% at d=86, 6.1% at d=95) that the 22-point grid stepped over.

**Other checks** (`data/2026-09-27_theory_numerics.json`, one key per item).
- The within-sector LP range of Δ₂ over m₀,m₁-preserving redistributions is [−0.92, 7.09]: the largest shift
  is 3.84 τ₂, and the range contains zero. The earlier committed case
  (`data/2026-08-24_within_sector_control.json`) has ‖δ‖₁ = 0.494 m₀, i.e. 24.7% of the weight transferred, and Δ₂ = 0.044 τ₂.
- The node Vandermonde has κ = 50 on nodes rescaled to [−1, 1] (1.2×10⁶ in raw units, a unit-dependent
  number). Zero guaranteed power comes from the null-space dimension, not from conditioning.
- A reconstruction that passes with tolerances is bracketed only by the tolerance-inflated extremal atom:
  3.6–5.7× W₁ at t = 12–20 on the U/t=4 measure; on the U/t=8 measure of the Christoffel figure,
  3.2–7.7× at t = 14–20 and 10–11× at the figure's t*. On that U/t=8 measure the exact max_t W_n = 1 (at
  the mean), 0.969, 0.886, 0.882 for n = 1…4.
- The Hausdorff bound m₂ ≤ 147.6 needs the weighted support [2.49, 13.92]t, which is oracle knowledge; the
  a-priori supports give 579.8 (full Fock range) or 243.3 (symmetry sector); the 573.7 of the 2026-09-05
  text used [0.148, 52.21]t. Both this bound and the Hankel bound are necessary conditions; the composition
  with the Wang/Mortimer certificates is not implemented.
- `m₋₁`: on the ring, `m₋₁ = ½⟨−T̂⟩ − D` is a 94% cancellation against the Kohn stiffness. At equal 2% bias,
  `m₁` and `m₂` fire at smaller injected fractions than `m₋₁` (ring 0.051 and 0.042 against 0.197; open
  chain 0.055 and 0.042 against 0.168). `m₋₁` fires first only for spurious poles at ≤ 0.5t.

---

## Known gaps

Not yet reproducible from this repository, or known defects:

- **Manuscript text.** `paper/main.tex` was revised on 2026-09-27 against the phase-A figures and numbers.
  It is not yet shortened or ported to the SciPost template, and it relies on the companion's arXiv v3; see
  [`paper/README.md`](paper/README.md).
- **`blind_addenda.py`, part (c):** the tie order of the sealed truncations was recovered from ARPACK's
  default start (first `eigsh` calls of a fresh process; scipy 1.13.1, numpy 1.26.4) and is stored in
  `data/2026-09-27_blind_addenda.json`. Reruns use the stored order after checking that it reproduces all 56
  sealed instances, and abort without writing if neither it nor the ARPACK rebuild does. The sealed files
  remain the canonical record.
- **`run_teeth_shared.py` (L=12, `argsort` order):** a rerun reproduces the two 5% crossings (d ≈ 15977 and
  43554) but moves the three smallest-d points of the sweep (m̄₁ by up to 6% at d = 180) against the
  committed 2026-08-26 record, which is kept (`data/2026-08-26_teeth_shared_state.json`,
  `paper/figs/teeth_shared.dat`).
- **`fig_sqw`:** the rasters are placed at their true extents, and the display windows clip 6.3% of the
  broadened charge-grid weight (above 10t; 8.5% at the worst q) and 0.25% of the spin weight (above 2.4t)
  (`src/check_sqw_extent.py` → `data/2026-09-27_sqw_extent_check.json`).
- **`paper/figs/bracket_maxW.dat`** (written by `export_bracketing_dat.py`, read by no figure) holds grid
  maxima that differ from the exact values in `christoffel_maxW.dat`; do not quote it.
- **Provenance hash:** `delta0_reference_mc.py` records the SHA-256 of
  `data/_superseded/analytic_coverage_leak.json`, which changes with line endings (a CRLF checkout gives a
  different hash for identical content).
- **Notebook:** the device reference Monte Carlo (`delta0_reference_mc.py`, qiskit, 7–20 min) is behind a
  switch that is off by default, so `make reproduce` does not rerun it (`make device` does).
- **Archive:** no Zenodo DOI yet; README thumbnails not regenerated.

---

## Repository structure

```
.
├── src/                                # exact-diagonalization engine + screen + figure exporters
│   ├── verify.py                       #   fast smoke test (moments + the joint battery catching a lone-m1 miss), seconds
│   ├── hubbard_ed.py                   #   vendored ED utilities (Jordan–Wigner Hubbard chain; numpy/scipy only)
│   ├── spectral_lanczos.py             #   sector Lanczos + Haydock: A(k,ω), S(q,ω), S^zz, current probe, teeth
│   ├── interval_moment*.py, interval_battery.py   # shot-budget interval test, MC, joint battery (deterministic order)
│   ├── blind_*.py, separating_counts.py#   the sealed blinded protocol (seal → generate → classify → score), addenda, figure recount
│   ├── hardware_matched_job_L8.py …    #   ibm_fez L=8 job, noiseless baseline, retention-matched control
│   ├── delta0_reference_mc.py, export_device_dat.py  # reference Monte Carlo for Δ₀ and the fig_device data
│   ├── export_heron_dat.py             #   fig_heron data from data/heron_spectral.json (formatting only)
│   ├── small_checks.py, within_sector_lp.py, …       # theory numerics of the revision (one JSON, one key per item)
│   ├── nk_*.py                         #   the pre-registered n_k device experiment (gate-negative; never run)
│   ├── run_*.py, export_*.py           #   demonstrations and one data-driven generator per figure fragment
│   ├── cache/*_L12.npz                 #   committed L=12 Lanczos caches read by export_akw/sqw/teeth
│   └── _superseded/                    #   quarantined prototypes and replaced versions (see its README)
├── data/                               # committed results (*.json), sealed records, real ibm_fez retained counts
│   ├── README.md                       #   which dated files supersede which (nothing is deleted)
│   └── _superseded/                    #   superseded model, precursor-program data, .dat copies (see its README)
├── notebooks/
│   └── 00_Reproduce_Everything.ipynb   #   narrated run of the committed generators, phase-A scripts included
├── paper/                              # the manuscript, revised 2026-09-27 (LaTeX source + PDF, 38 pp) — see paper/README.md
├── docs/
│   ├── REPRODUCE.md                    #   figure/number → script → exact command, phase-A outputs, known gaps
│   ├── FILE_INDEX.md                   #   every file, described
│   └── img/_superseded/                #   the retired 2026-08-22 README thumbnails
├── _superseded/                        # the 2026-09-02 arXiv bundle and upload guide (not submitted)
├── requirements.txt · Makefile · CITATION.cff · LICENSE
```

---

## Quick start

```bash
# 1. environment
python -m venv .venv && source .venv/bin/activate      # (Windows: .venv\Scripts\activate)
pip install -r requirements.txt

# 2. smoke test in seconds (numpy/scipy only): the current-probe moments, the shot-budget thresholds,
#    and the joint battery rejecting, via m2, a truncation that the lone first moment passes (d=97)
python src/verify.py        # or: make verify

# 3. the 2026-09-27 recomputations (numpy/scipy/mpmath; about 20 min in total)
make revision               # interval scripts, blinded-test addenda, theory numerics, Gauss-law data
make device                 # reference Monte Carlo for Δ₀ + fig_device data (needs qiskit; 7–20 min)

# 4. regenerate the native figure data (uses the committed src/cache/*_L12.npz)
make figures

# 5. everything: the narrated notebook (includes the 2026-09-27 scripts; about 45 min)
make reproduce

# 6. build the manuscript (needs a TeX distribution)
make paper                  # -> paper/main.pdf   (uses the committed main.bbl; no bibtex)
```

Scripts write their JSON into `data/` and their plotted data into `paper/figs/`, so `git status` after a
run shows whether a committed number moved. Sealed records (`data/prereg.*`, `data/blind_*.json`,
`data/manifest_nk_device_v1.*`, `data/2026-09-05_nk_v4_verdict.*`) are read-only: their writers refuse to
overwrite them.

---

## Reproducibility at a glance

| Layer | Reproducible here? | How |
|---|---|---|
| **The committed generators, one pass** | ✅ narrated, phase-A scripts included | `make reproduce` (the notebook) |
| **Exact diagonalization** (`A(k,ω)`, `S(q,ω)`, `S^zz`, the current-response tests, the teeth, the Christoffel bounds) | ✅ locally | `python src/<script>.py`; see `docs/REPRODUCE.md` |
| **Interval-moment test + joint battery** | ✅ deterministic order; every-depth scan and tie-break sensitivity | `make intervals` |
| **Blinded, pre-registered test** | ✅ re-scores the sealed record; addenda re-derive from it | `make blind` (seal `data/prereg.sha256`) |
| **Theory numerics of the revision** | ✅ one JSON, one key per item | `make theory` (≈15 min) |
| **Figures** | ✅ all data-driven figures | `make figures`; `fig_device` data from `make device` |
| **IBM Heron** | ✅ re-analysis of the retained L=8 counts (needs `qiskit`) | counts in `data/heron_counts_matched_L8_*.json` (acquired by `src/hardware_matched_job_L8.py`); reference Monte Carlo `src/delta0_reference_mc.py`; the companion's L=6 raw counts are not deposited |

Every result but the hardware acquisition is an **exact classical simulation**; re-acquiring device data
needs an IBM Quantum account.

---

## 🔒 Security note (please read before pushing)

This repository contains **no secrets**. Any IBM Quantum hardware re-run uses **your own** credentials
(`QiskitRuntimeService.save_account` locally, or an environment variable) — never commit a token.

---

## Citation

If you use this work, please cite it (see [`CITATION.cff`](CITATION.cff)):

```bibtex
@unpublished{BonillaVargas2026MomentTests,
  title  = {Reconstruction-independent moment tests for quantum-computed
            dynamical spectra},
  author = {Bonilla Vargas, Nicol\'as},
  note   = {Manuscript in revision; code and data at
            https://github.com/nicolasbonilla/reconstruction-independent-moment-tests},
  year   = {2026}
}
```

Once a DOI or journal reference exists, this entry and `CITATION.cff` will be updated.

## License

- **Code** (`*.py`, `notebooks/`, `Makefile`): [MIT](LICENSE).
- **Paper text and figures** (`paper/`): [CC-BY-4.0](LICENSE).
