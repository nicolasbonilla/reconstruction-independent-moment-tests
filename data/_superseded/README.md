# Superseded data (kept as a record, not evidence)

Moved here 2026-09-26 with `git mv`; nothing was deleted. The 2026-09-05 text still quotes the σ of
`analytic_coverage_leak.json` (9.3×10⁻⁵ at 50k to 9.2×10⁻⁴ at 4k) as "shot noise is negligible"; the
revision replaces it with the Monte Carlo spread. No other manuscript number rests on these files. One
live script reads one of them: `src/delta0_reference_mc.py` (2026-09-27) reads `analytic_coverage_leak.json`
only to report how far its σ falls below the Monte Carlo spread. Records superseded later but kept in
`data/` itself are listed in [`../README.md`](../README.md).

- `analytic_coverage_leak.json` — superseded shot-noise model of the coverage leak (`fp_rate = 1.0`); retained as a record, not evidence. Its σ is 34–195× smaller than the spread of 2000 noiseless replicas at the device's raw shots (`data/2026-09-27_delta0_reference_mc.json`, `budgets.<N>.sd_over_analytic_sigma`). Generator: `src/analytic_coverage_leak.py` (re-runs write here).

## `precursor_program_2026-08-17/`

- `2026-08-17_faf_L8_U8.json` — precursor-program data; not used by the manuscript; generator not in this repository.
- `2026-08-17_gate0_scaling_L8.json` — precursor-program data; not used by the manuscript; generator not in this repository.
- `2026-08-17_gate1a_powercheck.json` — precursor-program data; not used by the manuscript; generator not in this repository.
- `2026-08-17_groundtruth_L8_U8.json` — precursor-program data; not used by the manuscript; generator not in this repository.
- `2026-08-17_window_scan_L8.json` — precursor-program data; not used by the manuscript; generator not in this repository.
- `2026-08-18_chigap_bakeoff.json` — precursor-program data; not used by the manuscript; generator not in this repository.
- `2026-08-18_chigap_partial.json` — precursor-program data; not used by the manuscript; generator not in this repository.
- `2026-08-18_contrib1_calibration.json` — precursor-program data; not used by the manuscript; generator not in this repository.
- `2026-08-18_gate_A.json` — precursor-program data; not used by the manuscript; generator not in this repository.
- `2026-08-18_gate_B.json` — precursor-program data; not used by the manuscript; generator not in this repository.

## `dat_duplicates_2026-08-22/`

Copies, as of 2026-08-22, of `.dat` files in `paper/figs/`. Checked 2026-09-27 (line endings ignored):
- 19 are identical to the live file in `paper/figs/`, which is the authoritative copy (read by a figure,
  written by a committed exporter, or both). One of them, `bracket_maxW.dat`, is written by
  `src/export_bracketing_dat.py` but read by no figure, and its grid maxima differ from the exact values in
  the live `christoffel_maxW.dat`; do not quote it.
- 2 (`christoffel_maxW.dat`, `christoffel_poles.dat`) are **no longer** identical to the live files, which
  `src/export_christoffel_dat.py` rewrote on 2026-09-27 (exact maxima; weighted atoms only). They equal the
  replaced versions kept as `paper/figs/_superseded/christoffel_maxW_grid400.dat` and
  `paper/figs/_superseded/christoffel_poles_with_zero_weight_rows.dat`.
- 4 (`falsifier_spectrum`, `falsifier_sweep`, `gauss`, `heron_sweep`) duplicate orphans that no figure
  reads, now in `paper/figs/_superseded/` (see its README).

- `akw_edges.dat` — duplicate of `paper/figs/akw_edges.dat` (authoritative).
- `akw_extent.dat` — duplicate of `paper/figs/akw_extent.dat` (authoritative).
- `akw_gapedge.dat` — duplicate of `paper/figs/akw_gapedge.dat` (authoritative).
- `akw_screen.dat` — duplicate of `paper/figs/akw_screen.dat` (authoritative).
- `bracket_F.dat` — duplicate of `paper/figs/bracket_F.dat` (authoritative).
- `bracket_env_n2.dat` — duplicate of `paper/figs/bracket_env_n2.dat` (authoritative).
- `bracket_env_n4.dat` — duplicate of `paper/figs/bracket_env_n4.dat` (authoritative).
- `bracket_maxW.dat` — duplicate of `paper/figs/bracket_maxW.dat` (read by no figure; see above).
- `christoffel_Wn.dat` — duplicate of `paper/figs/christoffel_Wn.dat` (authoritative).
- `christoffel_floor.dat` — duplicate of `paper/figs/christoffel_floor.dat` (authoritative).
- `christoffel_maxW.dat` — the 400-point-grid maxima; equals `paper/figs/_superseded/christoffel_maxW_grid400.dat` (the live `paper/figs/christoffel_maxW.dat` now holds the exact maxima).
- `christoffel_poles.dat` — with three zero-weight rows; equals `paper/figs/_superseded/christoffel_poles_with_zero_weight_rows.dat` (the live file lists weighted atoms only).
- `dcp_lower.dat` — duplicate of `paper/figs/dcp_lower.dat` (authoritative).
- `dcp_upper.dat` — duplicate of `paper/figs/dcp_upper.dat` (authoritative).
- `falsifier_spectrum.dat` — duplicate of `paper/figs/_superseded/falsifier_spectrum.dat` (an orphan: no figure reads it).
- `falsifier_sweep.dat` — duplicate of `paper/figs/_superseded/falsifier_sweep.dat` (an orphan: no figure reads it).
- `gauss.dat` — duplicate of `paper/figs/_superseded/gauss.dat` (an orphan: no figure reads it).
- `heron_exact.dat` — duplicate of `paper/figs/heron_exact.dat` (authoritative).
- `heron_hw.dat` — duplicate of `paper/figs/heron_hw.dat` (authoritative).
- `heron_sweep.dat` — duplicate of `paper/figs/_superseded/heron_sweep.dat` (an orphan: no figure reads it).
- `momentcone_point.dat` — duplicate of `paper/figs/momentcone_point.dat` (authoritative).
- `sqw_deltac.dat` — duplicate of `paper/figs/sqw_deltac.dat` (authoritative).
- `sqw_extent.dat` — duplicate of `paper/figs/sqw_extent.dat` (authoritative).
- `sqw_spinpeak.dat` — duplicate of `paper/figs/sqw_spinpeak.dat` (authoritative).
- `teeth.dat` — duplicate of `paper/figs/teeth.dat` (authoritative).
