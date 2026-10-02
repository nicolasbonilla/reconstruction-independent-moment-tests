<div align="center">

# Reconstruction-independent moment tests
### *for sample-based spectral functions*

**Nicolás Bonilla Vargas** &nbsp;[![ORCID](https://img.shields.io/badge/ORCID-0009--0006--6155--4391-A6CE39?logo=orcid&logoColor=white)](https://orcid.org/0009-0006-6155-4391)<br>
Universidad Nacional de Colombia · SRH University, Munich · Daita AI, Munich

[![Manuscript](https://img.shields.io/badge/manuscript-SciPost%20Physics%20Core%20submission%20version-lightgrey.svg)](paper/README.md)
[![License: MIT](https://img.shields.io/badge/code-MIT-green.svg)](LICENSE)
[![License: CC BY 4.0](https://img.shields.io/badge/paper%2C%20data%2C%20docs-CC--BY--4.0-lightgrey.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](requirements.txt)

</div>

> **Status (2026-10-02).** [`paper/`](paper/) holds the manuscript prepared for submission to SciPost
> Physics Core; it has not yet been submitted. Every number and figure in it traces to a script in `src/`
> and a data file in `data/` ([`docs/REPRODUCE.md`](docs/REPRODUCE.md)); the few that do not regenerate
> exactly from source are listed under [Known gaps](#known-gaps).

Sample-based quantum diagonalization reconstructs a dynamical spectral function `A(ω)` in a sampled
determinant subspace, and checking it against moments read off `A(ω)` is circular. This work compares
the reconstruction's moments `m̄ₖ` with full-operator moments `m̂ₖ = ⟨0_S|J†(H−e₀)ᵏJ|0_S⟩` evaluated on the
same subspace ground state. Their difference is exactly the probe's leak out of the subspace at order `k`.
The test is a **necessary condition**: failing it rejects a reconstruction, passing it does not certify one.
Every moment in this work is evaluated classically; on IBM Heron the device supplies only sampled supports.

**Companion study.** N. Bonilla Vargas, *Captured weight and boundary leakage bound the error of
sample-based spectral functions*, [arXiv:2608.16436](https://arxiv.org/abs/2608.16436). It bounds the
reconstruction error at exact-diagonalization scale given the exact ground state; its captured weight `w`
gives the coverage residual here, `Δ₀ = m₀(1 − w)`, the `k = 0` leak. This work adds tests computable from
the subspace state alone. Review of the field: *Machine learning for sample-based quantum diagonalization:
a review of generative configuration recovery and the classical-simulability frontier*,
[arXiv:2608.05314](https://arxiv.org/abs/2608.05314).

---

## Results

1. **The leak hierarchy** (paper, Sec. 2.3, Proposition 1). On the subspace ground state,
   `m̂ₖ − m̄ₖ` is the leak of `J|0_S⟩` out of the subspace at order `k`: the missed weight at `k = 0`, the
   squared boundary coupling at `k = 2` when no weight is missed, and zero through order `2K+1` when the
   subspace contains the order-`K` Krylov space of the probe state. It needs sparse matrix–vector products
   from `|0_S⟩` only. Off an eigenstate the operator form matters: on a truncated `L=6` ring
   `⟨J(H−e₀)J⟩ = 4.558` but `½⟨[J,[H,J]]⟩ = 3.717`; at `τ₁ = 0.229` the first rejects a reconstruction whose
   first moment is 11% low and the second passes it (`src/estimator_form_offeigenstate.py`, key
   `R6_estimator_forms_off_eigenstate` of `data/2026-09-27_theory_numerics.json`).
2. **What the moments detect** (Sec. 3). Under determinant truncation the zeroth-order leak (the
   companion's captured-weight deficit) is the operative channel. When weight moves in frequency at fixed
   total, `m₀` is blind and `m₁`, `m₂` reject it, in the doped `L=12` current response from a moved fraction
   of a few percent (key `R8_inverse_moment_sensitivity`). The optical f-sum moment `m₋₁` shifts most in relative terms on the ring, but a
   94% cancellation against the Kohn stiffness makes it reject last unless the spurious weight lies at
   `0.5t` or below (`run_inverse_moment_falsifier.py`, with `--r8` for the noise-normalized comparison).
3. **A pre-registered blinded calibration in simulation** (Sec. 4; oracle estimator, doped `L=6` ring,
   300 hash-sealed instances). Rejected: all 56 subspace truncations, 25 of 27 weight-changing and 20 of 34
   weight-preserving spurious features, and none of 71 under-converged Krylov reconstructions (the 56
   with two or more Lanczos nodes reproduce `m₀`–`m₂` exactly; the 15 one-node ones miss `m₂` by
   2.4–4.2%, within `τ₂`). False-positive rate: 1/112 as implemented, 40/112 = 35.7% (Wilson
   95% [27.4, 44.9]%) at the noise scale of the sealed text, under which noise also rejects 23 of the 71
   Krylov reconstructions (`src/blind_score.py`,
   `src/blind_addenda.py` → `data/2026-08-24_blind_harness_score.json`,
   `data/2026-09-27_blind_addenda.json`).
4. **Exact limits of a pass** (Sec. 5). With `m₀`–`m₂` the guaranteed power against fixed-pole weight
   redistributions is zero for any spectrum with four or more distinct poles; on the five-pole `L=6` current
   spectrum a redistribution that moves 72% of the weight passes (`within_sector_lp.py`). Rayleigh–Ritz on
   the full ground-state support of the calibration ring passes (`Δ₂ = 0.34 τ₂`) with a relative `L₁` error
   of 0.93; over `U/t = 3–8` it passes with 0.68–0.99 (`full_support_rr.py` →
   `data/2026-10-02_full_support_rr.json`). The shared-state leak scores correlate only moderately with the
   `L₁` error (Spearman 0.31–0.44).
5. **Coverage is not fidelity, on device data** (Sec. 6). On supports sampled by `ibm_fez` for an `L=8`
   addition sector, `Δ₀` grows from 0.0044 to 0.37 as coverage falls from 94% to 35%, but it rewards support
   spreading: at 50k shots per circuit the device lies below all 2000 noiseless replicas of its own circuits,
   a uniform in-sector sampler leaves less weight uncovered than the device at every budget, and
   `Δ₀ = 0.0044` coexists with a relative `L₁` error of 0.35 (`delta0_reference_mc.py` →
   `data/2026-09-27_delta0_reference_mc.json`).

**Scope.** The test rejects; it does not certify. Thresholds are calibrated in simulation for the current
probe of one doped `L=6` ring, on synthetic errors of classically known truth; all systems are at
exact-diagonalization scale (`L ≤ 12`, plus an `L=24` DMRG check of code scalability). The realizable
estimator is the deterministic shared-state one, whose residual is the leak, not the distance to the true
moment (at `L=24` its `m̂₀` lies 10–15.5% below the DMRG `m₀`); the shot-noise rejection rates come from the
oracle estimator. No moment has been estimated from device counts. For a surrogate bond-basis probe the
forecast depolarizing bias of `m₁` at 295 CZ is 37.6% (static local noise) to 116.4% (global) of the signal
without mitigation (`data/2026-08-30_offdiag_noise_forecast.json`). The sealed on-device `n_k` test was not
run: on the calibrations checked (2026-09-03 to 2026-09-05) no 12-qubit chain met its sealed coherence
floors, and no job was submitted
(`data/2026-09-05_nk_device_attempt_record.json`). Each limit is stated with its evidence in the paper,
Sec. 7.2.

---

## Device jobs

All on `ibm_fez`, SamplerV2 with measurement twirling, Pauli gate twirling and dynamical decoupling; no
readout-error mitigation. Only the sampled support enters any result; no moment is estimated from counts.

| Created (UTC) | Job IDs | What | Deposit | Role in the paper |
|---|---|---|---|---|
| 2026-08-29 | `da983qse74ec73ajfr2g` (50k), `da98k16rbfbs73chckvg` (30k), `da98icsjbipc73fej78g` (16k), `da98cksjbipc73fej1tg` (4k) | `L=8`, 16 qubits, 7 number-conserving sampler circuits per job | raw counts `data/heron_counts_matched_L8_<id>.json`; submitted circuits in `data/sampler_cz_count.json` | Sec. 6, Fig. device |
| 2026-08-21 | `da44t8u1vhnc73fkie30` (U = 4), `da459juaa69c739jcb80` (U = 3.5) | `L=6`, 12 qubits, 7 circuits × 35 000 shots | raw REST records `data/heron_gatetp_2026-08-21/`; reanalysis `src/gate_tp_analysis.py` → `data/2026-10-02_gatetp_reanalysis.json` | none: the retained support is the whole 300-configuration sector under either bit order, so a support-based test returns the same verdict without the device data |
| — | none | sealed `n_k` test (manifest v1) | `data/2026-09-05_nk_device_attempt_record.json`, `src/nk_record/` | not run; 0 jobs, 0 shots |

The companion's own `L=6` run was post-selected in reversed bit order, so every determinant it retained
came from a device error; its spectra are in `data/heron_spectral.json` (copied from the companion's
deposit), its raw counts are not deposited here.

---

## Known gaps

Full list with numbers in [`docs/REPRODUCE.md`](docs/REPRODUCE.md#known-gaps). In short:

- **Blinded truncations** regenerate from the seed only with the stored tie order (degenerate `|ψ₀|²`
  shells); `blind_addenda.py` re-validates it against the sealed record on every run. The sealed files are
  the canonical record.
- **`run_teeth_shared.py`** (`L=12`, `argsort` order): a rerun reproduces both 5% crossings but moves the
  three smallest-`d` points of the sweep; the committed 2026-08-26 record is the one plotted.
- **`fig_sqw`** display windows clip 6.3% of the broadened charge weight and 0.25% of the spin weight.
- **Not deposited:** the companion's raw `L=6` counts; the Colab notebook that submitted the 2026-08-21 jobs
  (its analysis is re-implemented and checked in `gate_tp_analysis.py`).
- **Archive:** the Zenodo DOI is added at submission.

---

## Repository structure

```
.
├── src/                                # engine, tests, generators (one per number family and figure)
│   ├── verify.py                       #   smoke test + record check of the numbers main.tex prints from five records
│   ├── hubbard_ed.py, spectral_lanczos.py  # ED utilities; sector Lanczos + Haydock (A(k,ω), S(q,ω), S^zz, current probe)
│   ├── estimator_form_offeigenstate.py #   O_k against the commutator form off an eigenstate
│   ├── interval_moment*.py, interval_battery.py   # shot-budget interval test, every-depth scan, Monte Carlo
│   ├── blind_*.py, separating_counts.py#   sealed blinded protocol (seal → generate → classify → score), addenda, recount
│   ├── within_sector_lp.py, christoffel_tolerance_lp.py, full_support_rr.py, small_checks.py  # limits of a pass
│   ├── hardware_matched_job_L8.py, delta0_reference_mc.py, export_device_dat.py  # L=8 device run, reference samplers
│   ├── sampler_cz_count.py, offdiag_*.py, bond_moment_estimator.py  # device-native estimator cost and noise forecast
│   ├── gate_tp_analysis.py             #   reanalysis of the 2026-08-21 L=6 jobs
│   ├── beyond_ed_dmrg_catch.py         #   L=24 DMRG code-scalability check (physics-tenpy)
│   ├── nk_*.py, nk_record/             #   sealed n_k test (gate, seal, job) and the record of its device attempt
│   ├── run_*.py, export_*.py           #   demonstrations and one data-driven generator per figure fragment
│   ├── cache/*_L12.npz                 #   committed L=12 Lanczos caches
│   └── _superseded/                    #   quarantined prototypes and replaced versions (see its README)
├── data/                               # dated results (*.json), sealed records, raw ibm_fez counts and job records
│   ├── README.md                       #   which file supersedes which; generator of each record
│   ├── heron_gatetp_2026-08-21/        #   raw records of the two 2026-08-21 jobs (byte for byte)
│   └── _superseded/                    #   superseded model, precursor data, earlier string versions
├── notebooks/00_Reproduce_Everything.ipynb   # narrated run of the generators
├── paper/                              # the manuscript (SciPost Physics Core class), see paper/README.md
├── docs/
│   ├── REPRODUCE.md                    #   figure/number → script → command; sealed-record checks; known gaps
│   ├── FILE_INDEX.md                   #   every file, described
│   └── PROCESS_RECORD.md               #   design history of the sealed tests and of the manuscript
├── _superseded/                        # an unsubmitted 2026-09-02 upload bundle (record only)
├── requirements.txt · Makefile · CITATION.cff · .zenodo.json · LICENSE
```

---

## Quick start

```bash
python -m venv .venv && source .venv/bin/activate      # (Windows: .venv\Scripts\activate)
pip install -r requirements.txt

python src/verify.py        # or: make verify   (seconds; ALL CHECKS PASS expected)
make revision               # interval tests, blinded test + addenda, theory numerics incl. full-support RR, Gauss law (~20 min)
make figures                # the .dat fragments of every data-driven figure
make device                 # reference Monte Carlo for Δ₀ + fig_device data (qiskit; 7–30 min)
make hardware               # CZ counts, 295-CZ noise forecast, 2026-08-21 reanalysis (qiskit, qiskit-aer, qiskit-ibm-runtime; no QPU)
make dmrg                   # L=24 DMRG check (physics-tenpy; ~50 min)
make reproduce              # the narrated notebook (~45–60 min)
make paper                  # pdflatex, bibtex, pdflatex ×2 -> paper/main.pdf
```

Scripts write JSON into `data/` and plotted data into `paper/figs/`. After a run,
`git diff --ignore-cr-at-eol --stat data paper/figs` lists the outputs that changed; run time, run date and
environment fields change on every run, numbers should not. Sealed records (`data/prereg.*`,
`data/blind_*.json`, `data/manifest_nk_device_v1.*`, `data/2026-09-05_nk_v4_verdict.*`) are read-only: their
writers refuse to overwrite them.

## Reproducibility at a glance

| Layer | How |
|---|---|
| Exact diagonalization (`A(k,ω)`, `S(q,ω)`, current response, truncation sweeps, Christoffel bounds) | `python src/<script>.py`; map in `docs/REPRODUCE.md` |
| Leak identity, estimator forms, limits of a pass | `make theory` |
| Interval test and joint battery | `make intervals` |
| Blinded, pre-registered test | `make blind` (seal `data/prereg.sha256`) |
| Figures | `make figures`; `fig_device` data from `make device` |
| Device records (no QPU) | `make device`, `make hardware`; raw counts and job records in `data/` |
| Sealed `n_k` test | `make nkrecord` (offline); the gate scripts need `make nk-sealed-dir` first |

Every result other than the device sampling is an exact classical computation. Re-acquiring device data
needs an IBM Quantum account of your own.

## Security note

This repository contains no credentials. A device re-run uses your own saved account
(`QiskitRuntimeService.save_account` locally) or an environment variable; never commit a token. The IBM job
records in `data/heron_gatetp_2026-08-21/` carry the submitting account's `user_id` (an IBMid account
identifier), as returned by the API.

## Citation

See [`CITATION.cff`](CITATION.cff).

```bibtex
@unpublished{BonillaVargas2026MomentTests,
  title  = {Reconstruction-independent moment tests for sample-based spectral functions},
  author = {Bonilla Vargas, Nicol\'as},
  note   = {Code and data: https://github.com/nicolasbonilla/reconstruction-independent-moment-tests},
  year   = {2026}
}
```

## License

- **Code** (`*.py`, `notebooks/`, `Makefile`): [MIT](LICENSE).
- **Paper, data and documentation** (`paper/`, `data/`, `docs/`, READMEs): [CC BY 4.0](LICENSE).
- `paper/SciPost.cls` and `paper/SciPost_bibstyle.bst` are SciPost's template files under their own terms.
