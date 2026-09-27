# Superseded scripts — not used to reproduce the manuscript

These scripts are earlier prototypes. They do **not** generate the manuscript's figures and in places
contradict its numbers. They are kept only as a record (moved here 2026-09-26; nothing was deleted).
Their paths are not maintained: `FIGSRC` pointed at a machine-local copy of the companion's figures and
is now `None`, and outputs would land under `src/06_results/`, `src/07_figures/` (gitignored).

| Script | Why it is stale or circular | Live source for that figure / number |
|---|---|---|
| `run_moment_bound_theorem.py` | Prints `max_t W_4 = 0.445` on a coarse 60-point grid with no degeneracy aggregation. The exact value is `0.882243` (root-finding, since 2026-09-27: `paper/figs/christoffel_maxW.dat`, key `R9_christoffel_max_Wn` of `data/2026-09-27_theory_numerics.json`); the 400-point grid used before gave `0.875`. | `export_christoffel_dat.py` (exact `max_t Wₙ`), `export_bracketing_dat.py`, `markov_krein_window.py` (`fig_christoffel`, `fig_bracketing`) |
| `run_heron_screen.py` | Circular: a reconstruction-vs-reconstruction comparison (moments integrated from the reconstructed `A(ω)` against an injected distortion of the same reconstruction), not a device-error catch. It also reads a stale 37-point `heron_hw.dat` from a machine-local path. | `fig_heron` plots `paper/figs/heron_{exact,hw}.dat`, written by `src/export_heron_dat.py` from `data/heron_spectral.json` (byte-identical to the committed files) |
| `make_akw_figure.py` | Stale `L=6` raster (`η = 0.28`, Gaussian-blur corruption); a different figure from the manuscript's `L=12` panel. | `spectral_lanczos.py` → `export_akw_dat.py` (`fig_akw`) |
| `make_sqw_figure.py` | Stale `L=6` raster (`η_c = 0.30`, `η_s = 0.16`) on penalty-projected ED. | `spectral_lanczos.py` → `export_sqw_dat.py` (`fig_sqw`) |
| `make_paper_figs.py` | Produced the 2026-08-22 prototype rasters, now in `docs/img/_superseded/` (see the README there). | the native PGFPlots figures in `paper/figs/*.tex` |
| `export_falsifier_dat.py` | Writes `falsifier_{spectrum,sweep}.dat` for `fig_falsifier`, which is no longer in the manuscript; those files are in `paper/figs/_superseded/`. | — |

`paper_style.py` (in `src/`) is the plotting identity these raster scripts used; no live script needs it.

| Script (added 2026-09-27) | Why superseded | Live source |
|---|---|---|
| `necessary_sufficient_composition_2026-08.py` | The 2026-08 version of `src/necessary_sufficient_composition.py`; it generated `data/necessary_sufficient_composition.json` (kept unchanged as the record). It called the Hausdorff/support bound "sufficient (Wang/Mortimer SDP)" and took the support from all full-Fock eigenvalues including zero-weight states (m2 <= 573.68). | `src/necessary_sufficient_composition.py` (plan R5/B11): two necessary conditions, four support choices from a-priori to oracle (weighted support gives m2 <= 147.6); key `R5_hausdorff_weighted_support` of `data/2026-09-27_theory_numerics.json` |
