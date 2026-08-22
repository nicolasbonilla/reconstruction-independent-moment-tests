# -*- coding: utf-8 -*-
"""Build fig_akw assets from the sector-Lanczos/Haydock spectral cache (spectral_lanczos.py):
raster-field PNGs (paperhot) for the TRUE A(k,omega) and a WRONG reconstruction (illustrative injected
corruption modeled on the determinant-truncation weight-misassignment family of fig_teeth), the per-k
independent-vs-reconstruction first-moment residual (the screen), the Mott-gap edges mu+/mu-, and the
band-edge anchor points for the gap arrow. Native pgfplots overlays axes/labels/colorbar/guides.

The corruption rebalances upper/lower Hubbard-band INTENSITY per momentum so the reconstructed band
centroid flattens to (1-alpha)*eps_k: preserves per-k m0 EXACTLY and the empty gap A(k,E_F)=0, fails m1."""
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
L = int(os.environ.get('AKW_L', '12'))
GAMMA = 0.80                       # mild display gamma: reveal the fainter high-|omega| Hubbard structure
ALPHA = 0.25                       # corruption strength: reconstructed centroid -> (1-alpha)*eps_k
PAPERHOT = LinearSegmentedColormap.from_list('paperhot', [
    (22/255, 8/255, 5/255), (90/255, 22/255, 11/255), (176/255, 54/255, 12/255),
    (217/255, 142/255, 50/255), (246/255, 196/255, 83/255), (251/255, 239/255, 199/255)])


def main():
    d = np.load(os.path.join(CACHE, f'akw_L{L}.npz'))
    kk, ww, A = d['kk'], d['wg'], d['A']
    gap = float(d['gap']); mu_plus = gap / 2.0; mu_minus = -gap / 2.0     # edges measured from E_F (mid-gap)
    U = float(d['U'])
    # exact, reconstruction-INDEPENDENT first moment (anticommutator sum rule, half filling): m1hat(k)=eps_k
    krad = kk * np.pi
    m1hat = -2.0 * np.cos(krad)

    # ---- illustrative injected corruption: per-k inter-Hubbard-band centroid flattening ----
    Aw = A.copy()
    lo = ww < 0.0; hi = ww > 0.0
    for ik in range(Aw.shape[0]):
        col = Aw[ik]
        Wlo = col[lo].sum(); Whi = col[hi].sum()
        if Wlo < 1e-12 or Whi < 1e-12:
            continue
        S = Wlo + Whi
        clo = (col[lo] * ww[lo]).sum() / Wlo
        chi = (col[hi] * ww[hi]).sum() / Whi
        target = (1.0 - ALPHA) * m1hat[ik]
        dd = chi - clo
        a = S * (chi - target) / (Wlo * dd)
        b = S * (target - clo) / (Whi * dd)
        col[lo] *= a; col[hi] *= b

    # ---- screen: independent m1hat vs the WRONG reconstruction's first moment ----
    def m1field(M):
        w = M / (M.sum(axis=1, keepdims=True) + 1e-12)
        return (w * ww[None, :]).sum(axis=1)
    r = np.abs(m1hat - m1field(Aw))
    m0err = float(np.abs(A.sum(axis=1) - Aw.sum(axis=1)).max())

    # wrap the periodic point k=2pi == k=0 so heatmaps + screen span a clean [0,2]
    kk_d = np.append(kk, 2.0)
    A_d = np.vstack([A, A[:1]]); Aw_d = np.vstack([Aw, Aw[:1]]); r_d = np.append(r, r[0])

    with open(os.path.join(OUT, 'akw_screen.dat'), 'w') as f:
        f.write("k resid\n"); [f.write(f"{a:.4f} {b:.4f}\n") for a, b in zip(kk_d, r_d)]
    with open(os.path.join(OUT, 'akw_edges.dat'), 'w') as f:
        f.write("mu_plus mu_minus gap\n%.4f %.4f %.4f\n" % (mu_plus, mu_minus, gap))
    # anchored band-edge points for the gap arrow: the direct Mott gap sits at k_F=pi/2
    with open(os.path.join(OUT, 'akw_gapedge.dat'), 'w') as f:
        f.write("k om\n0.5 %.4f\n0.5 %.4f\n" % (mu_plus, mu_minus))

    # ---- raster fields (paperhot, shared vmax, mild gamma), interpolated in k for continuity ----
    kfine = np.linspace(kk_d.min(), kk_d.max(), 360)
    vmax = float(np.percentile(A, 99.6))
    for M, name in [(A_d, 'akw_true'), (Aw_d, 'akw_wrong')]:
        Mf = interp1d(kk_d, M, axis=0, kind='linear')(kfine)
        norm = np.clip(Mf / vmax, 0, 1) ** GAMMA
        plt.imsave(os.path.join(OUT, name + '.png'), norm.T, cmap=PAPERHOT, vmin=0, vmax=1, origin='lower')
    with open(os.path.join(OUT, 'akw_extent.dat'), 'w') as f:
        f.write("kmin kmax wmin wmax\n%.4f %.4f %.4f %.4f\n" % (kk_d.min(), kk_d.max(), ww.min(), ww.max()))

    # PH-symmetry sanity: A(k,w)=A(k+pi,-w) at half filling
    Aroll = np.roll(A[:, ::-1], L // 2, axis=0)
    ph = float(np.abs(A - Aroll).sum() / (A.sum() + 1e-12))
    print("L=%d | mu+=%.3f mu-=%.3f gap=%.3f | PH-asym=%.3f | m0err=%.2e | vmax=%.4f | screen max=%.3ft"
          % (L, mu_plus, mu_minus, gap, ph, m0err, vmax, r.max()))
    print("Delta_1(k):", " ".join("%.2f" % x for x in r))
    print("wrote akw_true/wrong.png, akw_screen/edges/gapedge/extent.dat ->", OUT)


if __name__ == '__main__':
    main()
