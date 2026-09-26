# -*- coding: utf-8 -*-
"""Build fig_falsifier assets from the L=12 sector-Lanczos/Haydock current cache (spectral_lanczos.py).
The true current spectral function A_J(omega) of the doped 1D Hubbard ring is essentially all mid-IR /
Mott-scale weight (omega ~ 10t); the coherent Drude sits at omega=0 (charge stiffness) and is not part of
the regular spectrum shown. The WRONG reconstruction invents a SPURIOUS low-omega (Drude-like) peak that
steals a fraction of the mid-IR weight -- a documented reconstruction artifact -- preserving the total
weight m0=<J^2> EXACTLY while failing the first current moment m1 (wrong mean frequency)."""
import os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')
OUT = os.path.normpath(os.path.join(HERE, '..', 'paper', 'figs'))
L = int(os.environ.get('FALS_L', '12'))
OM_SPLIT = 6.0          # divide the (empty) low-omega region from the mid-IR/Mott band
OM_SPUR = 2.5           # where the spurious low-omega peak is invented
ETA = 0.35


def main():
    d = np.load(os.path.join(CACHE, f'current_L{L}.npz'))
    om, w = d['om'], d['w']; m0 = float(d['m0']); m1 = float(d['m1'])
    hi = om >= OM_SPLIT
    W_hi = w[hi].sum()

    def wrong_poles(f):
        w_w = w.copy(); w_w[hi] *= (1.0 - f)
        return np.append(om, OM_SPUR), np.append(w_w, f * W_hi)   # spurious peak carries the stolen weight

    wgrid = np.linspace(0, 15, 560)      # match the axis xmax exactly -> no curve drawn past the frame

    def broaden(oms, wts):
        A = np.zeros_like(wgrid)
        for e, wt in zip(oms, wts):
            if wt > 1e-10:
                A += wt * (ETA / np.pi) / ((wgrid - e) ** 2 + ETA ** 2)
        return A

    A_true = broaden(om, w)
    ow, ww = wrong_poles(0.40)
    A_wrong = broaden(ow, ww)
    with open(os.path.join(OUT, 'falsifier_spectrum.dat'), 'w') as f:
        f.write("omega Atrue Awrong\n")
        for x, a, b in zip(wgrid, A_true, A_wrong):
            f.write(f"{x:.4f} {a:.6f} {b:.6f}\n")

    # residual sweep: fraction of mid-IR weight misplaced into the spurious low-omega peak
    with open(os.path.join(OUT, 'falsifier_sweep.dat'), 'w') as f:
        f.write("f r0 r1\n")
        for fr in np.linspace(0.0, 0.5, 11):
            ow, ww = wrong_poles(fr)
            m0w = float(ww.sum()); m1w = float((ow * ww).sum())
            r0 = abs(m0w - m0) / m0 * 100.0
            r1 = abs(m1w - m1) / m1 * 100.0
            f.write(f"{fr*100:.3f} {r0:.6f} {r1:.4f}\n")

    print(f"L={L}: m0=<J^2>={m0:.4f}  m1={m1:.4f}  poles={len(om)}")
    print(f"mid-IR weight (om>={OM_SPLIT})={W_hi:.4f}  low-om weight={m0-W_hi:.4f}  A_true max={A_true.max():.3f} A_wrong max={A_wrong.max():.3f}")
    print("wrote falsifier_spectrum.dat, falsifier_sweep.dat ->", OUT)


if __name__ == '__main__':
    main()
