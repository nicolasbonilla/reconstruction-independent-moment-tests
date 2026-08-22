# -*- coding: utf-8 -*-
"""Build fig_sqw assets from the sector-Lanczos/Haydock cache (spectral_lanczos.py): raster-field PNGs
(paperhot, display gamma) for the charge S(q,omega) and spin S^zz(q,omega) structure factors, the des
Cloizeaux-Pearson two-spinon edges, the charge-gap onset Delta_c, and the EXTRACTED spinon-peak points
(dominant peak per q) that track the dCP dispersion. Each channel on its own frequency window (charge ~U,
spin ~J=4t^2/U) so neither the charge band nor the gapless spinon continuum is crushed."""
import os, sys
import numpy as np
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')
OUT = os.path.normpath(os.path.join(HERE, '..', 'paper', 'figs'))
L = int(os.environ.get('SQW_L', '12'))
GAMMA = 0.62
PAPERHOT = LinearSegmentedColormap.from_list('paperhot', [
    (22/255, 8/255, 5/255), (90/255, 22/255, 11/255), (176/255, 54/255, 12/255),
    (217/255, 142/255, 50/255), (246/255, 196/255, 83/255), (251/255, 239/255, 199/255)])


def main():
    d = np.load(os.path.join(CACHE, f'sqw_L{L}.npz'))
    qq, wc, ws, Sc, Ss = d['qq'], d['wc'], d['ws'], d['Sc'], d['Ss']
    U = float(d['U']); WC_MAX = float(d['wc_max']); WS_MAX = float(d['ws_max'])
    qfine = np.linspace(qq.min(), qq.max(), 320)
    for Sd, name in [(Sc, 'sqw_charge'), (Ss, 'sqw_spin')]:
        Sfine = interp1d(qq, Sd, axis=0, kind='linear')(qfine)
        vmax = float(np.percentile(Sd, 99.5))
        norm = np.clip(Sfine / vmax, 0, 1) ** GAMMA
        plt.imsave(os.path.join(OUT, name + '.png'), norm.T, cmap=PAPERHOT, vmin=0, vmax=1, origin='lower')
    # charge-gap onset: lowest omega carrying >5% of the charge peak
    prof_c = Sc.max(axis=0); onset = float(wc[np.argmax(prof_c > 0.05 * prof_c.max())])
    with open(os.path.join(OUT, 'sqw_deltac.dat'), 'w') as f:
        f.write("onset\n%.4f\n" % onset)
    # des Cloizeaux-Pearson two-spinon edges (effective Heisenberg J=4t^2/U)
    J = 4.0 / U
    qpi = np.linspace(qq.min(), qq.max(), 240); qrad = qpi * np.pi
    lower = (np.pi * J / 2) * np.abs(np.sin(qrad)); upper = np.pi * J * np.abs(np.sin(qrad / 2))
    with open(os.path.join(OUT, 'dcp_lower.dat'), 'w') as f:
        f.write("q om\n"); [f.write(f"{a:.4f} {b:.4f}\n") for a, b in zip(qpi, lower)]
    with open(os.path.join(OUT, 'dcp_upper.dat'), 'w') as f:
        f.write("q om\n"); [f.write(f"{a:.4f} {b:.4f}\n") for a, b in zip(qpi, upper)]
    # extracted dominant spinon peak per q (tracks the dCP dispersion) -> data-on-guide dots
    with open(os.path.join(OUT, 'sqw_spinpeak.dat'), 'w') as f:
        f.write("q om\n")
        for iq, q in enumerate(qq):
            om_pk = float(ws[np.argmax(Ss[iq])])
            f.write(f"{q:.4f} {om_pk:.4f}\n")
    with open(os.path.join(OUT, 'sqw_extent.dat'), 'w') as f:
        f.write("qmin qmax wc_max ws_max\n%.4f %.4f %.4f %.4f\n" % (qq.min(), qq.max(), WC_MAX, WS_MAX))
    print("L=%d | charge onset Delta_c=%.3f | J=%.3f | windows charge[0,%.1f] spin[0,%.1f] | nq=%d"
          % (L, onset, J, WC_MAX, WS_MAX, len(qq)))
    print("wrote sqw_charge/spin.png, dcp_*.dat, sqw_deltac/spinpeak/extent.dat ->", OUT)


if __name__ == '__main__':
    main()
