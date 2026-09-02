# Reproduce every number and figure

Every non-hardware result in the paper is an exact classical computation reproducible from this
repository with `numpy`/`scipy` alone. Run scripts from `src/` (`cd src && python <script>.py`), or run
the whole pipeline end to end in [`../notebooks/00_Reproduce_Everything.ipynb`](../notebooks/00_Reproduce_Everything.ipynb).

## Fastest check (seconds)

```bash
python src/verify.py        # or: make verify
```
Recomputes the current-probe moments `m₀, m₁, m₂` of the doped `L=6, U/t=4` ring, the shot-budget
interval at the real `ibm_fez` 50k-shot budget, and shows the joint `(m₀,m₁,m₂)`+Hankel battery
**rejecting** the truncation (`d=98`) that the lone first moment misses. Expected: **ALL CHECKS PASS**.

## Figure → engine → data

| Figure (paper) | What it shows | Engine / script | Data fragment(s) |
|---|---|---|---|
| Fig. 1 `fig_hero` | the reconstruction-independent moment loop (schematic) | native TikZ | — |
| Fig. 2 `fig_momentcone` | a moment sequence refutable by geometry | `export_momentcone_dat.py` | `momentcone_point.dat` |
| Fig. 3 `fig_akw` | `A(k,ω)` Mott map + per-`k` screen | `spectral_lanczos.run_akw` → `export_akw_dat.py` | `akw_*.dat`, `akw_true/wrong.png` |
| Fig. 4 `fig_sqw` | `S(q,ω)` gapped / `S^zz` gapless | `spectral_lanczos.struct_factors` → `export_sqw_dat.py` | `sqw_*.dat`, `dcp_*.dat` |
| Fig. 5 `fig_falsifier` | wrong distribution of correct total weight | `run_sumrule_falsifier.py` → `export_falsifier_dat.py` | `falsifier_*.dat` |
| Fig. 6 `fig_teeth` | independence vs circular control | `spectral_lanczos.run_teeth` → `export_teeth_dat.py` | `teeth.dat` |
| Fig. 7 `fig_gausslaw` | U(1) Gauss-law cross-domain falsifier | `run_gausslaw_falsifier.py` | `gauss.dat` |
| Fig. 8 `fig_bracketing` | the Gauss–Radau bracketing mechanism | `run_moment_bound_theorem.py` → `export_bracketing_dat.py` | `bracket_*.dat` |
| Fig. 9 `fig_christoffel` | the certified `Wₙ(t)` miss-distance | `run_moment_bound_theorem.py` → `export_christoffel_dat.py` | `christoffel_*.dat` |
| Fig. 10 `fig_heron` | the screen on real IBM Heron data | `run_heron_screen.py` | `heron_*.dat`, `heron_spectral.json` |

## Key numbers → script

| Number | Value | Where |
|---|---|---|
| current-probe moments `(m₀,m₁,m₂)` (L=6, U/t=4) | `0.5542, 4.2271, 33.3887` | `verify.py`, `interval_moment.py` |
| shot-budget 95% interval `z·δ₁` (ibm_fez, 2% bias) | `0.229` (5.4% of m₁) | `verify.py`, `interval_moment.py` |
| joint battery closes the `d=98` blind spot | `|Δm₂| = 3.28 > 1.85` | `verify.py`, `interval_battery.py` |
| single-particle Mott edges / gap (L=12) | `μ± = ±2.484t`, `Δ ≈ 4.97t` | `spectral_lanczos.run_akw` |
| charge onset `Δc` (L=12) | `≈ 5.71t` | `spectral_lanczos.struct_factors` |

## Hardware

The IBM Heron result runs the screen on the **cached** spectral data of the companion study
([arXiv:2608.16436](https://arxiv.org/abs/2608.16436)), shipped here as `data/heron_spectral.json`
(`ibm_fez`, 12 qubits, 50k shots, TREX + Pauli twirling + dynamical decoupling, full `|S|=300/300`
sector). Re-acquiring the device data needs an IBM Quantum account; the classical screen does not.
