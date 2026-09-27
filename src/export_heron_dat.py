# -*- coding: utf-8 -*-
"""Write the fig_heron coordinates from committed JSON only (plan item R11).

Reads data/heron_spectral.json (the companion's L=6, U/t=4 reconstructed spectra on a 600-point grid,
eta = 0.15; ibm_fez job d9s16avpemts73ct6g8g) and writes, in paper/figs/:

  heron_exact.dat   w A     the exact target A(omega)            (JSON keys grid, A_exact)
  heron_hw.dat      w Ahw   the reconstruction from that run     (JSON keys grid, A_hw)

with w to 4 decimals and A to 5 decimals, the format of the committed files (this script reproduces them
byte for byte). No computation beyond formatting; the numbers are the JSON's.

What the figure can and cannot show (see docs/REPRODUCE.md, Hardware): in that run the retained subspace
covered the whole 300-configuration sector, so A_hw equals A_exact to rounding (JSON hw_relL1 = 0), which is
exact by coverage, not a fidelity result. The run was post-selected in reversed bit order, so all of its
retained determinants came from device errors; it used 5x10^4 shots per circuit on 7 circuits (3.5x10^5 in
total), measurement twirling, Pauli gate twirling and dynamical decoupling, and no readout-error mitigation;
its raw counts are not deposited.

Run from src/:  python export_heron_dat.py     (under a second)
"""
import os, json

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_JSON = os.path.normpath(os.path.join(HERE, "..", "data", "heron_spectral.json"))
OUT = os.path.normpath(os.path.join(HERE, "..", "paper", "figs"))


def write_dat(path, header, xs, ys):
    with open(path, "w", encoding="utf-8") as f:
        f.write(header + "\n")
        for x, y in zip(xs, ys):
            f.write(f"{x:.4f} {y:.5f}\n")


def main():
    d = json.load(open(SRC_JSON, encoding="utf-8"))
    grid, a_ex, a_hw = d["grid"], d["A_exact"], d["A_hw"]
    if not (len(grid) == len(a_ex) == len(a_hw)):
        raise SystemExit("heron_spectral.json: grid and spectra have different lengths")
    write_dat(os.path.join(OUT, "heron_exact.dat"), "w A", grid, a_ex)
    write_dat(os.path.join(OUT, "heron_hw.dat"), "w Ahw", grid, a_hw)
    dmax = max(abs(x - y) for x, y in zip(a_ex, a_hw))
    print(f"backend {d['backend']}, job {d['job_id']}: L={d['L']}, U={d['U']}, eta={d['eta']}, "
          f"sector {d['nsector']}, |S| = {d['hw_S']}, stored relative L1 {d['hw_relL1']}")
    print(f"{len(grid)} grid points in [{grid[0]:.4f}, {grid[-1]:.4f}]; max |A_hw - A_exact| = {dmax:.1e} "
          f"(equal by full-sector coverage, not a fidelity result)")
    print("wrote heron_exact.dat, heron_hw.dat ->", OUT)


if __name__ == "__main__":
    main()
