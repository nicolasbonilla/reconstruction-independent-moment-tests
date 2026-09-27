# -*- coding: utf-8 -*-
"""Export native-figure .dat for fig_christoffel with the board's physics corrections:
(1) aggregate degenerate eigenvalues by unique frequency (the '4 poles' were one 6-fold atom);
(2) dense-grid (400 pts) W_n(t) curves (christoffel_Wn.dat);
(3) W_n at the dominant pole via mpmath high precision (double underflows) -> saturation floor = residue.
(4) [2026-09-27, R9] max_t W_n EXACTLY, by root-finding on the degree-2n polynomial K_n(t) = 1/W_n(t)
    (mpmath, 60 digits): christoffel_maxW.dat now holds these values (the 400-point-grid maxima it held
    before, 0.999843/0.964330/0.885610/0.875125, are kept in paper/figs/_superseded/christoffel_maxW_grid400.dat
    and in the JSON below). max_t W_n is invariant under the affine rescaling t <-> omega.
(5) [2026-09-27, m1, m2] Gauss-Radau rules with a node at t OUTSIDE the weighted hull (all weights > 0,
    atom = W_n(t)); the rescaled frame [-1,1] <-> [0.148, 52.21] t (full-Fock range incl. zero-weight states)
    vs the weighted support; W_1 = Cantelli.
Writes paper/figs/christoffel_{Wn,maxW,floor,poles}.dat and the keys R9_christoffel_max_Wn,
m1_radau_outside_hull, m2_rescaled_frame of data/2026-09-27_theory_numerics.json. Run: cd src && python export_christoffel_dat.py"""
import os, sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)   # hubbard_ed.py is vendored in src/
import hubbard_ed as H
from run_sumrule_falsifier import current_operator
OUT = os.path.normpath(os.path.join(HERE, '..', 'paper', 'figs'))
L, U = 6, 8.0


def christoffel_widths(x, wts, orders, tq):
    nmax = max(orders)
    pkm1_n = np.zeros_like(x); pk_n = np.ones_like(x)
    pkm1_t = np.zeros_like(tq); pk_t = np.ones_like(tq)
    norm2 = np.empty(nmax + 2); norm2[0] = float(np.sum(wts))
    K_t = pk_t ** 2 / norm2[0]
    out = {0: 1.0 / K_t}
    for k in range(nmax):
        alpha = float(np.sum(wts * x * pk_n ** 2) / norm2[k])
        beta = norm2[k] / norm2[k - 1] if k >= 1 else 0.0
        pkp1_n = (x - alpha) * pk_n - beta * pkm1_n
        pkp1_t = (tq - alpha) * pk_t - beta * pkm1_t
        norm2[k + 1] = float(np.sum(wts * pkp1_n ** 2))
        if norm2[k + 1] <= 1e-300:
            break
        K_t = K_t + pkp1_t ** 2 / norm2[k + 1]
        pkm1_n, pk_n = pk_n, pkp1_n
        pkm1_t, pk_t = pk_t, pkp1_t
        if (k + 1) in orders:
            out[k + 1] = 1.0 / K_t
    return out


def christoffel_at_point_mp(x, wts, tstar, nmax, dps=60):
    """High-precision Christoffel W_n(tstar)=1/K_n via monic Stieltjes recurrence in mpmath."""
    mp.mp.dps = dps
    xs = [mp.mpf(float(v)) for v in x]; ws = [mp.mpf(float(v)) for v in wts]
    t = mp.mpf(float(tstar))
    pkm1 = [mp.mpf(0)] * len(xs); pk = [mp.mpf(1)] * len(xs)
    pkm1_t = mp.mpf(0); pk_t = mp.mpf(1)
    norm2 = mp.fsum(ws)
    K = pk_t ** 2 / norm2
    Wn = {0: 1.0 / K}
    prev_norm2 = norm2
    for k in range(nmax):
        alpha = mp.fsum([ws[i] * xs[i] * pk[i] ** 2 for i in range(len(xs))]) / prev_norm2
        beta = prev_norm2 / norm2_km1 if k >= 1 else mp.mpf(0)
        pkp1 = [(xs[i] - alpha) * pk[i] - beta * pkm1[i] for i in range(len(xs))]
        pkp1_t = (t - alpha) * pk_t - beta * pkm1_t
        nn = mp.fsum([ws[i] * pkp1[i] ** 2 for i in range(len(xs))])
        if nn <= mp.mpf(10) ** (-dps + 5):
            break
        K = K + pkp1_t ** 2 / nn
        pkm1, pk = pk, pkp1
        pkm1_t, pk_t = pk_t, pkp1_t
        norm2_km1 = prev_norm2; prev_norm2 = nn
        Wn[k + 1] = float(1.0 / K)
    return Wn


def _pmul(p, q):
    r = [mp.mpf(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            r[i + j] += x * y
    return r


def _padd(p, q):
    n = max(len(p), len(q))
    return [(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)]


def _peval(p, t):
    return mp.fsum([c * mp.mpf(t) ** k for k, c in enumerate(p)])


def christoffel_K_poly(a, b, n):
    """K_n(t) = sum_{j<=n} pi_j(t)^2/||pi_j||^2 as low-first mpmath coefficients (monic recurrence a, b of
    markov_krein_window.recurrence; ||pi_j||^2 = b_0 b_1 ... b_j with b_0 = m0)."""
    pis = [[mp.mpf(1)], [-a[0], mp.mpf(1)]]
    for k in range(1, n):
        nxt = _padd(_pmul([-a[k], mp.mpf(1)], pis[k]), [-b[k] * x for x in pis[k - 1]])
        pis.append(nxt)
    K = [mp.mpf(0)]; nrm = mp.mpf(1)
    for j in range(n + 1):
        nrm = nrm * b[j]
        K = _padd(K, [x / nrm for x in _pmul(pis[j], pis[j])])
    return K


def exact_max_W(a, b, n):
    """global max over t in R of W_n(t) = 1/K_n(t): real critical points of K_n (roots of K_n')."""
    K = christoffel_K_poly(a, b, n)
    dK = [k * K[k] for k in range(1, len(K))]
    roots = mp.polyroots(dK[::-1], maxsteps=2000, extraprec=400) if len(dK) > 2 else [-dK[0] / dK[1]]
    best = None
    for r in roots:
        if abs(mp.im(r)) < mp.mpf(10) ** (-30):
            rr = mp.re(r); Wv = 1 / _peval(K, rr)
            if best is None or Wv > best[1]:
                best = (rr, Wv)
    return float(best[0]), float(best[1]), len(roots)


def theory_numerics(uom, uw, wn, xs, a_e, b_e, low_levels, runtime0):
    """R9 + m1 + m2 records (called from main)."""
    from small_checks import record_theory_numerics
    import markov_krein_window as mkw
    mp.mp.dps = 60
    idx = np.argsort(wn)[::-1]
    weighted = wn > 1e-12
    n_w = int(weighted.sum())
    tmean = float(np.sum(wn * xs)); mu = float(np.sum(wn * uom)); var = float(np.sum(wn * uom ** 2)) - mu ** 2
    exact = {}; grid = {}
    for n in (1, 2, 3, 4):
        a, b = mkw.recurrence(xs, wn, n)
        targ, Wmax, nroots = exact_max_W(a, b, n)
        exact[str(n)] = dict(max_W=Wmax, argmax_t_rescaled=targ, argmax_omega=a_e + (targ + 1) * (b_e - a_e) / 2,
                             argmax_inside_minus1_1=bool(-1 <= targ <= 1), n_critical_points=nroots,
                             check_W_at_argmax_recurrence=float(mkw.christoffel_W(a, b, targ, n)))
    for npts in (400, 4000, 40000):
        tq = np.linspace(-0.999, 0.999, npts)
        Wc = christoffel_widths(xs, wn, [1, 2, 3, 4], tq)
        grid[str(npts)] = {str(n): float(np.nanmax(Wc[n])) for n in (1, 2, 3, 4)}
    a1, b1 = mkw.recurrence(xs, wn, 1)
    W1_mean = float(mkw.christoffel_W(a1, b1, tmean, 1))
    Wdom = {}
    for n in (1, 2, 3, 4, 5, 6):
        a, b = mkw.recurrence(xs, wn, n)
        Wdom[str(n)] = float(mkw.christoffel_W(a, b, float(xs[idx[0]]), n))
    R9 = dict(measure='L=6, U/t=8 doped current (m0=1), rescaled to [-1,1] with the full-Fock range',
              n_atoms_all=int(len(uom)), n_weighted_atoms=n_w,
              method='roots of dK_n/dt, K_n = 1/W_n the degree-2n Christoffel polynomial (mpmath 60 digits)',
              exact_max=exact, grid_max=grid,
              W1_at_mean=W1_mean, W1_at_mean_equals_m0=bool(abs(W1_mean - 1.0) < 1e-12),
              t_mean_rescaled=tmean, omega_mean=mu,
              W_at_dominant_atom=Wdom, dominant_residue=float(wn[idx[0]]),
              dominant_atom_omega=float(uom[idx[0]]), dominant_atom_t=float(xs[idx[0]]),
              christoffel_maxW_dat='rewritten with exact_max (was grid 400: see grid_max["400"])',
              plan_expectation='max W4 = 0.875 on a 400-point grid, 0.882 on 4000; max W1 = m0 = 1 at the mean')
    # ---- m1: Gauss-Radau with a node outside the weighted hull ----
    hull = [float(xs[weighted].min()), float(xs[weighted].max())]
    rad = []
    for n in (1, 2, 3, 4):
        a, b = mkw.recurrence(xs, wn, n)
        for tout in (-1.5, -0.95, -0.6, 0.0, 0.6, 1.5):
            nodes, wts = mkw.gauss_radau(a, b, 1.0, tout, n)
            at = [float(wt) for nd, wt in zip(nodes, wts) if abs(float(nd) - tout) < 1e-9]
            lam = float(mkw.christoffel_W(a, b, tout, n))
            merr = max(abs(float(np.sum(wn * xs ** p)) - float(mp.fsum([wt * nd ** p for nd, wt in zip(nodes, wts)])))
                       for p in range(2 * n + 1))
            rad.append(dict(n=n, t=tout, outside_weighted_hull=bool(tout < hull[0] or tout > hull[1]),
                            all_weights_positive=bool(all(float(x) > 0 for x in wts)),
                            atom_at_t=at[0] if at else None, W_n_t=lam,
                            atom_equals_W_rel_err=(abs(at[0] - lam) / lam) if at else None,
                            max_moment_err_deg_le_2n=merr, nodes=[float(x) for x in nodes]))
    m1 = dict(claim='the Christoffel-Markov-Stieltjes equality needs no convex-hull hypothesis: a Gauss-Radau rule '
                    'with a node at any real t exists with positive weights and its atom at t equals W_n(t)',
              weighted_hull_rescaled=hull, rows=rad,
              all_positive=bool(all(r['all_weights_positive'] for r in rad)),
              max_atom_rel_err=float(max(r['atom_equals_W_rel_err'] for r in rad)),
              max_moment_err=float(max(r['max_moment_err_deg_le_2n'] for r in rad)))
    # ---- m2: rescaled frame ----
    om_w = uom[weighted]
    tmap = lambda t: a_e + (t + 1) * (b_e - a_e) / 2
    W1n = lambda om: 1.0 / (1.0 + (om - mu) ** 2 / var)
    cant = []
    for om_t in (12.0, 14.0, 20.0, tmap(0.3), tmap(0.6), tmap(0.9)):
        trs = (2 * om_t - (a_e + b_e)) / (b_e - a_e)
        cant.append(dict(omega=om_t, t_rescaled=trs, W1_recurrence=float(mkw.christoffel_W(a1, b1, trs, 1)),
                         cantelli_sigma2_over_sigma2_plus_dev2=var / (var + (om_t - mu) ** 2),
                         abs_diff=abs(float(mkw.christoffel_W(a1, b1, trs, 1)) - var / (var + (om_t - mu) ** 2))))
    m2 = dict(rescaling_range_omega=[a_e, b_e],
              rescaling_note='[a,b] = min/max positive-frequency eigenvalue of the FULL Fock space (zero-weight '
                             'states included), not the weighted support',
              weighted_support_omega=[float(om_w.min()), float(om_w.max())],
              weighted_support_rescaled=hull,
              fraction_of_minus1_1_beyond_last_weighted_pole=(1 - hull[1]) / 2,
              t_star_to_omega={str(t): tmap(t) for t in (0.3, 0.6, 0.9)},
              mean_omega=mu, variance=var, W1_at_omega_14=W1n(14.0),
              W1_equals_cantelli=cant,
              lowest_levels_H_minus_0N=low_levels,
              ground_triplet_note='the lowest level of H (mu=0) is 3-fold (spin triplet at N=4); the aggregated '
                                  'current measure does not depend on M (J is a spin scalar; checked in R4: '
                                  'sector N_up=2 vs mixed-M ground vector)')
    rt = time.time() - runtime0
    record_theory_numerics('R9_christoffel_max_Wn', R9, 'export_christoffel_dat.py', rt, seed=None)
    record_theory_numerics('m1_radau_outside_hull', m1, 'export_christoffel_dat.py', rt, seed=None)
    record_theory_numerics('m2_rescaled_frame', m2, 'export_christoffel_dat.py', rt, seed=None)
    return exact


def main():
    t_start = time.time()
    c = H.build_operators(L); cd = [ci.getH() for ci in c]
    Ham = H.hubbard_hamiltonian(L, t=1.0, U=U, pbc=True, c=c)
    Ntot = sum((cd[m] @ c[m] for m in range(2 * L)), sp.csr_matrix((2 ** (2 * L),) * 2, dtype=complex))
    J = current_operator(L, c)
    for mu_ in [0.0, -1.0, -2.0]:
        cand = eigsh(Ham - mu_ * Ntot, k=1, which='SA',
                     v0=np.random.default_rng(0).standard_normal(Ham.shape[0]) + 0j)[1][:, 0]   # fixed start
        if float(np.real(np.vdot(cand, Ntot @ cand))) <= L - 1.5:
            psi0 = cand; break
    E0 = float(np.real(np.vdot(psi0, Ham @ psi0)))
    Jpsi = J @ psi0
    Ev, Vv = np.linalg.eigh(Ham.toarray())
    om = Ev - E0; w = np.abs(Vv.conj().T @ Jpsi) ** 2
    keep = om > 1e-6; om, w = om[keep], w[keep]

    # AGGREGATE degenerate eigenvalues by unique frequency (board fix)
    order = np.argsort(om); om, w = om[order], w[order]
    uom, uw = [], []
    for o, wi in zip(om, w):
        if uom and abs(o - uom[-1]) < 1e-6:
            uw[-1] += wi
        else:
            uom.append(o); uw.append(wi)
    uom = np.array(uom); uw = np.array(uw)
    wn = uw / uw.sum()                                   # normalized residues (m0=1)
    a, b = uom.min(), uom.max()
    xs = (2 * uom - (a + b)) / (b - a)                   # rescale to [-1,1]

    # dominant + next atoms
    idx = np.argsort(wn)[::-1]
    print(f"distinct atoms: {len(uom)}; dominant residue={wn[idx[0]]:.4f} at t={xs[idx[0]]:.4f} (omega={uom[idx[0]]:.4f})")
    print(f"next residue={wn[idx[1]]:.4f} at t={xs[idx[1]]:.4f} (omega={uom[idx[1]]:.4f})")
    with open(os.path.join(OUT, 'christoffel_poles.dat'), 'w') as f:
        f.write("t residue\n")
        for i in idx[:8]:
            if wn[i] > 1e-12:          # 2026-09-27: list only weighted atoms (zero-weight rows were noise-ordered)
                f.write(f"{xs[i]:.6f} {wn[i]:.6f}\n")

    # (a) dense-grid W_n(t) for n=1..4
    xq = np.linspace(-0.999, 0.999, 400)
    orders = [1, 2, 3, 4]
    Wc = christoffel_widths(xs, wn, orders, xq)
    with open(os.path.join(OUT, 'christoffel_Wn.dat'), 'w') as f:
        f.write("t W1 W2 W3 W4\n")
        for j in range(len(xq)):
            f.write(f"{xq[j]:.5f} " + " ".join(f"{np.clip(Wc[n][j],0,1):.6f}" for n in orders) + "\n")
    print("400-grid max_t W_n:", {n: round(float(np.nanmax(Wc[n])), 6) for n in orders})
    # (b) EXACT max_t W_n (R9, root-finding) -> christoffel_maxW.dat ; also R9/m1/m2 JSON keys
    low = np.sort(eigsh(Ham, k=4, which='SA', v0=np.random.default_rng(0).standard_normal(Ham.shape[0]) + 0j,
                        tol=1e-12)[0]).tolist()
    exact = theory_numerics(uom, uw, wn, xs, float(a), float(b), low, t_start)
    with open(os.path.join(OUT, 'christoffel_maxW.dat'), 'w') as f:
        f.write("n maxW\n")
        for n in orders:
            f.write(f"{n} {exact[str(n)]['max_W']:.6f}\n")
    print("exact max_t W_n:", {n: round(exact[str(n)]['max_W'], 6) for n in orders})

    # (c) W_n at the dominant pole via mpmath -> saturates at residue
    tstar = float(xs[idx[0]])
    Wn_pole = christoffel_at_point_mp(xs, wn, tstar, nmax=10, dps=60)
    with open(os.path.join(OUT, 'christoffel_floor.dat'), 'w') as f:
        f.write("n Wn_at_pole\n")
        for n in sorted(Wn_pole):
            if n >= 1:
                f.write(f"{n} {Wn_pole[n]:.6f}\n")
    print("W_n at dominant pole (mpmath):", {n: round(Wn_pole[n], 4) for n in sorted(Wn_pole) if n >= 1})
    print("residue floor:", round(float(wn[idx[0]]), 4))
    print("wrote christoffel_*.dat ->", OUT)


if __name__ == '__main__':
    main()
