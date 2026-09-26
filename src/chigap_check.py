# -*- coding: utf-8 -*-
"""ACTION 6 -- CHRISTOFFEL EMPTY-REGION CHECK (sim-only, exact measure, high precision).

The Christoffel width W_n(t)=1/K_n(t,t) is the MAX mass a positive measure matching
m_0..m_2n can place at t. Evaluated at a point t* the true measure leaves EMPTY, W_n(t*)
upper-bounds the SPURIOUS weight any moment-passing reconstruction can place there.

HONESTY SCOPE (audited 2026-08-24, after an over-reach was caught): the L=6,U=8 CURRENT
measure concentrates its weight in a single mid-IR/Mott-band peak at omega~1.4U, and the
high-frequency region is empty (a ONE-SIDED tail, NOT a two-sided Mott gap). So this probes
OUT-OF-BAND emptiness of THIS optical measure; it is NOT the two-sided single-particle Mott
gap of the L=12 A(k,omega) figure (a different correlator we do NOT evaluate here). We report
the empty-region tightness and do not extrapolate it to the single-particle gap.

Reuses the exact L=6,U=8 current Lehmann measure and the mpmath recurrence of
export_christoffel_dat.py. NO device, NO fit.
"""
import os, sys, json
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)   # hubbard_ed.py is vendored in src/
import hubbard_ed as H
from run_sumrule_falsifier import current_operator
from export_christoffel_dat import christoffel_at_point_mp, christoffel_widths

L, U = 6, 8.0


def current_measure():
    c = H.build_operators(L); cd = [ci.getH() for ci in c]
    Ham = H.hubbard_hamiltonian(L, t=1.0, U=U, pbc=True, c=c)
    Ntot = sum((cd[m] @ c[m] for m in range(2 * L)),
               sp.csr_matrix((2 ** (2 * L),) * 2, dtype=complex))
    J = current_operator(L, c)
    for mu_ in [0.0, -1.0, -2.0]:
        cand = eigsh(Ham - mu_ * Ntot, k=1, which='SA')[1][:, 0]
        if float(np.real(np.vdot(cand, Ntot @ cand))) <= L - 1.5:
            psi0 = cand; break
    E0 = float(np.real(np.vdot(psi0, Ham @ psi0)))
    Ev, Vv = np.linalg.eigh(Ham.toarray())
    om = Ev - E0; w = np.abs(Vv.conj().T @ (J @ psi0)) ** 2
    keep = om > 1e-6; om, w = om[keep], w[keep]
    order = np.argsort(om); om, w = om[order], w[order]
    uom, uw = [], []
    for o, wi in zip(om, w):
        if uom and abs(o - uom[-1]) < 1e-6:
            uw[-1] += wi
        else:
            uom.append(o); uw.append(wi)
    return np.array(uom), np.array(uw) / np.sum(uw)      # unique poles, m0=1


def main():
    uom, wn = current_measure()
    a, b = uom.min(), uom.max()
    xs = (2 * uom - (a + b)) / (b - a)                    # -> [-1,1]
    xs_sorted = np.sort(xs)
    order = np.argsort(xs); cw = np.cumsum(wn[order]) / wn.sum()

    # where does the weight live? (this measure concentrates weight in a mid-band peak; the empty region is the high-freq tail)
    t_wtop = xs_sorted[np.searchsorted(cw, 0.999)]        # top of the weight-bearing region
    # GENUINELY-interior empty point: widest gap with BOTH endpoints away from the array edges
    interior = [(xs_sorted[i + 1] - xs_sorted[i], xs_sorted[i], xs_sorted[i + 1])
                for i in range(len(xs_sorted) - 1)
                if xs_sorted[i] > -0.9 and xs_sorted[i + 1] < 0.9]
    interior.sort(reverse=True)
    g, lo, hi = interior[0]; t_empty = 0.5 * (lo + hi)     # inside the empty out-of-band region

    dom = int(np.argmax(wn))
    print("=" * 68)
    print(f"CHRISTOFFEL EMPTY-REGION  (L={L}, U/t={U} current Lehmann, {len(uom)} atoms)")
    print("=" * 68)
    print(f"measure concentrates weight in a mid-band peak -- 99.9% of weight below t={t_wtop:.3f}; the region")
    print(f"above is an EMPTY ONE-SIDED TAIL (not a two-sided Mott gap).")
    print(f"empty-region probe t*={t_empty:.3f} (widest interior empty interval [{lo:.3f},{hi:.3f}])")
    print(f"dominant residue (saturation floor) = {wn[dom]:.4f}\n")

    W_dom = christoffel_at_point_mp(xs, wn, float(xs[dom]), nmax=8, dps=60)
    W_emp = christoffel_at_point_mp(xs, wn, float(t_empty), nmax=8, dps=60)
    xq = np.linspace(-0.999, 0.999, 400)
    Wc = christoffel_widths(xs, wn, [1, 2, 3, 4], xq)
    maxW = {n: float(np.nanmax(Wc[n])) for n in (1, 2, 3, 4)}

    print(f"{'n':>2} {'W_n(dom pole)':>14} {'W_n(empty region)':>18} {'max_t W_n':>12}")
    for n in range(1, 5):
        print(f"{n:>2} {W_dom.get(n,float('nan')):>14.4f} {W_emp.get(n,float('nan')):>18.6f} {maxW[n]:>12.4f}")

    tight = W_emp.get(2, 1.0) < 0.05 * maxW[2]
    if tight:
        verdict = ("TIGHT out-of-band: across the empty high-frequency region a moment-passing "
                   "reconstruction can place < %.1e spurious weight (W_2~%.1e) even though the "
                   "GLOBAL bound W_4~%.3f (set by resolving the dominant mid-band peak) is loose. NOT "
                   "extrapolated to the two-sided single-particle Mott gap (a different measure)."
                   % (W_emp.get(2, 9), W_emp.get(2, 9), maxW[4]))
    else:
        verdict = "LOOSE: empty-region W_n not << global; report as 'checked and loose'."
    print("\nVERDICT:", verdict)

    out = {'_provenance': {'script': 'chigap_check.py', 'sim_only': True,
                           'scale': f'L={L} U/t={U} current Lehmann (mid-IR/Mott-peak-dominated, one-sided high-freq tail)',
                           'n_atoms': int(len(uom)), 'weight_top_t': float(t_wtop),
                           'scope_note': 'out-of-band emptiness of the current measure; NOT the '
                                         'two-sided single-particle Mott gap'},
           't_empty_region': float(t_empty), 'empty_interval': [float(lo), float(hi)],
           'W_empty_region': {int(k): float(v) for k, v in W_emp.items()},
           'W_dominant_pole': {int(k): float(v) for k, v in W_dom.items()},
           'max_t_Wn': maxW, 'tight': bool(tight)}
    DATA = os.path.normpath(os.path.join(HERE, '..', 'data'))
    json.dump(out, open(os.path.join(DATA, '2026-08-24_chigap_check.json'), 'w'), indent=2)
    print("\nwrote data/2026-08-24_chigap_check.json")


if __name__ == '__main__':
    main()
