<div align="center">

# Reconstruction-independent moment tests
### *for quantum-computed dynamical spectra — a transportable necessary-condition screen*

**Nicolás Bonilla Vargas** &nbsp;[![ORCID](https://img.shields.io/badge/ORCID-0009--0006--6155--4391-A6CE39?logo=orcid&logoColor=white)](https://orcid.org/0009-0006-6155-4391)

[![Draft](https://img.shields.io/badge/draft-2026--09--02%2C%20under%20revision-lightgrey.svg)](paper/README.md)
[![License: MIT](https://img.shields.io/badge/code-MIT-green.svg)](LICENSE)
[![License: CC BY 4.0](https://img.shields.io/badge/paper-CC--BY--4.0-lightgrey.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](requirements.txt)

</div>

> **Status (2026-09-26).** The manuscript is in revision and is **not yet on arXiv** or in a journal.
> [`paper/`](paper/) holds the **2026-09-02 draft** (30 pp), which contains statements now being corrected
> (see [`paper/README.md`](paper/README.md)); it is not the final text.

> **One sentence.** A reconstruction-independent, **necessary-condition** screen for sample-based quantum
> spectral functions: low-order spectral moments m₀, m₁, m₂ (and m₋₁ on open chains) are evaluated as
> ground-state operator expectations, never read back off the reconstructed A(ω), and compared with the
> moments of A(ω). A pass is consistent with the reconstruction; it is not proof of it.

This repository holds the code and data behind the paper. Numbers and figures trace to a script and a
data file listed in [`docs/REPRODUCE.md`](docs/REPRODUCE.md) and [`docs/FILE_INDEX.md`](docs/FILE_INDEX.md);
the few that do not yet have a committed generator are listed under [Known gaps](#known-gaps-being-fixed-with-the-revision).
Most non-hardware results regenerate from source with `numpy`/`scipy`/`mpmath`; the noise-model forecasts
also need `qiskit-aer`, and the `L=24` DMRG check needs `physics-tenpy`.

> **Companion works.** Data & pipeline — *Dynamical spectral functions from bitstring-sampled quantum
> subspaces* ([arXiv:2608.16436](https://arxiv.org/abs/2608.16436)). Review — *Machine learning for
> sample-based quantum diagonalization: a review of generative configuration recovery and the
> classical-simulability frontier* ([arXiv:2608.05314](https://arxiv.org/abs/2608.05314)).

Figures: see the draft PDF in `paper/`; thumbnails will be regenerated from the native PGFPlots sources
after the figure revision.

---

## What the paper shows

1. **A reconstruction-independent screen, made operational.** A reconstructed spectrum `A(ω)` is checked
   against a battery of **necessary-condition** tests: static ground-state spectral moments
   `mₖ = ∫ωᵏA(ω)dω = ⟨Oₖ⟩` with `Oₖ = J·adₕᵏ(J)` (so `m₀ = ⟨J²⟩`, `m₁ = ½⟨[J,[H,J]]⟩`), evaluated
   **independently** of the reconstruction (in this work: in simulation, or classically on the sampled
   subspace state) and **never read back off `A(ω)`**. A pass is corroboration, never proof.

2. **Independence is what gives the screen teeth.** On sample-based-diagonalization subspace truncation the
   independent estimator fires while a back-computed (circular) control stays identically blind — an exact
   classical simulation of the doped `L=12` Hubbard ring, not a hardware result.

3. **Transportability across error classes and domains.** The same construction catches a spurious
   low-frequency (Drude-like) peak that preserves the total weight `m₀` yet fails the mean-frequency `m₁`;
   a complementary state-level check catches a gauge-symmetry-breaking error via the exact U(1) Gauss-law
   identity `⟨∑G²⟩ = 0`.

4. **What a pass means, exactly.** In exact arithmetic, the classical Christoffel/Markov–Krein bounds limit
   the cumulative weight a moment-matching measure can misplace. At the deployed order (`m₀`–`m₂`) the
   battery has no guaranteed power against within-sector redistribution: a `K=1` redistribution moving 49%
   of the weight passes (`data/2026-08-24_within_sector_control.json`).

5. **A shot-budget form.** The exact-moment idealization is replaced by a shot-budget interval test; a joint
   `(m₀,m₁,m₂)` + Hankel–Stieltjes **feasibility battery** catches the accidental single-moment miss found
   near `d≈97–98` (the committed `d=98` case depends on tie-breaking and does not reproduce from a fresh
   clone; `d=97` does, see Known gaps), forecast at the real `ibm_fez` 50k-shot budget.

**Honest scope (stated throughout).** The screen **falsifies, it does not certify**; its conditions are
**necessary, not sufficient**; it targets one error class (subspace truncation / weight misassignment), not
device noise in general; and **no quantum learning advantage is claimed** — on IBM Heron the device
supplies only the sampled support set; no spectral moment is estimated on hardware, and the discriminating
off-diagonal moments are simulation-only.

---

## Outcomes, including the negative ones

- **Sealed primary endpoint** (blinded, hash-sealed simulation, 300 instances): TPR = 101/188 = **53.7%**
  (Wilson 95% [46.6, 60.7]%) at FPR = 1/112 = 0.9% ([0.2, 4.9]%), with an unbiased simulated estimator
  (`data/2026-08-24_blind_harness_score.json`; seal `data/prereg.sha256`). Per class: truncation 56/56,
  spurious features 45/61, under-converged Krylov **0/71**.
- **Two sealed secondary predictions failed:** one-node Krylov, expected to be caught, was caught 0/15; the
  m₀-preserving spurious-feature subclass, sealed as designed-blind, was caught 20/34 (reported post hoc;
  the register was not revised).
- **Hardware supplies only the sampled support.** On retained `ibm_fez` counts (L=8 addition sector) Δ₀
  grows from 0.004 to 0.37 as coverage falls from 94% to 35%. The `m0_hat`/`m1_hat` fields in
  `data/heron_counts_matched_L8_*.json` are exact classical values. At 50k shots the reconstructed A(ω)
  still has relative L1 error 0.35 (`falsifier.rel_L1`) while Δ₀ = 0.0044. The L=8 jobs
  (`src/hardware_matched_job_L8.py`) ran SamplerV2 with measurement twirling, Pauli gate twirling and
  dynamical decoupling; no readout-error mitigation was applied (the draft's "TREX" is a misnomer).
  Bitstrings were read in the correct qubit order, and post-selection alone defines the support (no
  configuration recovery): 82,697 of 350,000 shots kept at 50k.
- **Retention-matched control:** the device keeps ≈23.6% of shots after post-selection, the noiseless
  simulation 100%. At equal kept samples the noiseless Δ₀ is 0.203 / 0.248 / 0.329 / 0.418 against the
  device's 0.0044 / 0.032 / 0.123 / 0.373 (50k/30k/16k/4k per circuit; `src/retention_matched_control.py`,
  `data/2026-09-05_retention_matched_control.json`). The device's small Δ₀ comes from its noise-widened
  support (at 50k, |S| = 3682 on the device against 2721 noiseless at equal kept samples and 3381 at equal
  raw shots). The wider support also gives the device the closer line shape: relative L1 0.35, against 0.60
  for the noiseless run at equal raw shots (`data/heron_counts_matched_L8_DRYRUN.json`). A small Δ₀ still
  does not certify A(ω). The manuscript's statistics for this comparison are being redone.
- **n_k device experiment: gate-negative, not run.** The Stage-0 gate passed and manifest v1 was sealed on
  2026-09-02 (`data/manifest_nk_device_v1.sha256`). Its day-of rule requires a 12-qubit chain with
  T1 ≥ 150 µs and T2 ≥ 100 µs on every qubit. The single permitted amendment (floors lowered to 100/70 µs,
  with readout and CZ caps) failed its sealed noisy-simulation gate on 2026-09-05 (AMENDMENT-DENIED; row
  R1ro0.018: canary UCB 0.02867 > sealed threshold 0.02667), so manifest v1 stands; it has not been run.
  A row identical except ε = 2.85×10⁻³ passed (R2); ε is a compound noise knob.
  `data/2026-09-05_nk_v4_verdict.json`.

---

## Known gaps (being fixed with the revision)

These are **not** reproducible from this repository yet, or are being recomputed:

- `fig_gausslaw` plots `paper/figs/gauss_state.dat`, which has no committed generator.
- `fig_heron` and `fig_device`: no exporter yet for `paper/figs/heron_{exact,hw}.dat` (from
  `data/heron_spectral.json`) or for the coordinates typed inline in `fig_device.tex` (from the retained counts).
- `fig_separating`: `src/separating_counts.py` recounts the figure from the sealed record, but the export of
  `paper/figs/sep_*.dat` from it is pending. The recount differs from `data/separating_demonstration.json`
  and from the draft's caption: the shaded region holds 26 instances, not 24 — the 24 spurious-feature
  instances plus one determinant truncation (id 162, so not every truncation sits at `m₀` residual > 1) and
  the clean false positive (id 184); `m₁` fires on 16 and `m₂` alone on 10. The exemplar id 222 is
  weight-preserving (`params.preserve_m0 = true`), not total-weight-changing as the caption says. A
  corrected `separating_demonstration.json` is pending.
- The hardware sampling-noise analysis (a reference Monte Carlo for Δ₀) is being redone; the draft's
  analytic σ understates the spread.
- The `m₂` bracket `123.3 ≤ m₂ ≤ 573.7` (`necessary_sufficient_composition.py`) is being revised.
- The `d=98` single-moment miss does not reproduce from a fresh clone. The determinant truncation
  depends on tie-breaking: fresh runs of `interval_moment.py` and `interval_battery.py` both give
  `|Δm₁| = 0.254` at d=98 (committed: 0.2625 and 0.2250), above the 2%-bias interval 0.229, so `m₁`
  alone rejects d=98 at 2% bias (it still misses at 4%). Until 2026-09-26 both scripts printed a
  hard-coded "d=98 missed by m₁" message; they now print the verdict computed in the run. Other
  truncation depths shift too (for example `|Δm₁|` at d=58 moves from 3.736 to 3.447), but no other
  verdict of either script changes. The miss itself exists: `verify.py`'s auto-scan finds one at d=97
  (`|Δm₁| = 0.155`, caught by `m₂`).
- Figure fixes (including the `fig_sqw` raster extent) and regenerated README thumbnails.

---

## Repository structure

```
.
├── src/                                # exact-diagonalization engine + screen + figure exporters
│   ├── verify.py                       #   ★ fast smoke test (moments + the joint battery catching a lone-m1 miss), seconds
│   ├── hubbard_ed.py                   #   vendored ED utilities (Jordan–Wigner Hubbard chain; numpy/scipy only)
│   ├── spectral_lanczos.py             #   sector Lanczos + Haydock: A(k,ω), S(q,ω), S^zz, current probe, teeth
│   ├── interval_moment.py              #   shot-budget interval-moment forecast (ibm_fez 50k-shot budget)
│   ├── interval_battery.py             #   joint (m0,m1,m2)+Hankel battery
│   ├── blind_*.py, separating_counts.py#   the sealed blinded protocol (seal → generate → classify → score) + figure recount
│   ├── hardware_matched_job_L8.py …    #   ibm_fez L=8 job, noiseless baseline, retention-matched control
│   ├── nk_*.py                         #   the pre-registered n_k device experiment (gate-negative; never run)
│   ├── run_*.py, export_*.py           #   demonstrations and one data-driven generator per figure fragment
│   ├── cache/*_L12.npz                 #   committed L=12 Lanczos caches read by export_akw/sqw/teeth
│   └── _superseded/                    #   quarantined prototypes (see its README)
├── data/                               # committed results (*.json), sealed records, real ibm_fez retained counts
│   └── _superseded/                    #   superseded model, precursor-program data, .dat duplicates (see its README)
├── notebooks/
│   └── 00_Reproduce_Everything.ipynb   #   narrated run of the committed generators (optional cells need qiskit-aer)
├── paper/                              # 2026-09-02 draft under revision (LaTeX source + PDF, 30 pp) — see paper/README.md
├── docs/
│   ├── REPRODUCE.md                    #   figure/number → script → exact command, and the known gaps
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

# 2. reproduce the headline in SECONDS (numpy/scipy only):
#    the current-probe moments, the shot-budget interval, and the joint battery REJECTING an
#    auto-calibrated truncation that the lone first moment misses.
python src/verify.py        # or: make verify

# 3. run the committed generators in one narrated pass (optional cells need qiskit-aer)
jupyter notebook notebooks/00_Reproduce_Everything.ipynb    # or headless: make reproduce

# 4. regenerate the native figure data from the exact engine (uses the committed src/cache/*_L12.npz)
make figures

# 5. build the 2026-09-02 draft (needs a TeX distribution, e.g. TeX Live / MiKTeX)
make paper                              # -> paper/main.pdf   (uses the committed main.bbl; no bibtex)
```

Scripts write their JSON into `data/` and their plotted data into `paper/figs/`, so `git status` after a
run shows whether a committed number moved. Sealed records (`data/prereg.*`, `data/blind_*.json`,
`data/manifest_nk_device_v1.*`, `data/2026-09-05_nk_v4_verdict.*`) are read-only: their writers refuse to
overwrite them.

---

## Reproducibility at a glance

| Layer | Reproducible here? | How |
|---|---|---|
| **The committed generators, one pass** | ✅ narrated (see [Known gaps](#known-gaps-being-fixed-with-the-revision)) | `notebooks/00_Reproduce_Everything.ipynb` (or `make reproduce`) |
| **Exact diagonalization** (`A(k,ω)`, `S(q,ω)`, `S^zz`, the current-response falsifier, the teeth, the Christoffel bounds) | ✅ locally | `python src/<script>.py`; see `docs/REPRODUCE.md` |
| **Interval-moment forecast + joint battery** | ✅ (the committed `d=98` case does not reproduce; Known gaps) | `src/interval_moment.py`, `src/interval_battery.py` |
| **Blinded, pre-registered test** | ✅ re-scores the sealed record | `src/blind_score.py` on `data/blind_*.json` (seal `data/prereg.sha256`) |
| **Figures** | ✅ most; four figures lack a committed generator or exporter (`fig_gausslaw`, `fig_heron`, `fig_device`, `fig_separating`; Known gaps) | `make figures` (native pgfplots from the exact engine) |
| **IBM Heron** | ✅ re-analysis of retained counts (needs `qiskit-aer`) | L=8 retained `ibm_fez` counts (`data/heron_counts_matched_L8_*.json`, acquired by `src/hardware_matched_job_L8.py`; re-analysed by `src/bootstrap_device_curve.py` and `src/aer_noiseless_baseline.py`, which import its falsifier) + retention-matched control (`src/retention_matched_control.py`); the companion's L=6 raw counts are not deposited. That L=6 run (`fig_heron`) was post-selected in reversed qubit order, without configuration recovery, so its 300/300 coverage came from device errors (disclosed by the companion repository on 2026-09-26); the L=8 job here reads the bits in the correct order |

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

Once an arXiv identifier, DOI or journal reference exists, this entry and `CITATION.cff` will be updated.

## License

- **Code** (`*.py`, `notebooks/`, `Makefile`): [MIT](LICENSE).
- **Paper text and figures** (`paper/`): [CC-BY-4.0](LICENSE).
