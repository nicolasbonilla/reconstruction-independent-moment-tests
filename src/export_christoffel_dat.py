# -*- coding: utf-8 -*-
"""Export native-figure .dat for fig_christoffel with the board's physics corrections:
(1) aggregate degenerate eigenvalues by unique frequency (the '4 poles' were one 6-fold atom);
(2) dense-grid (400 pts) max_t W_n(t);
(3) W_n at the dominant pole via mpmath high precision (double underflows) -> saturation floor = residue."""
import os, sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.normpath(os.path.join(HERE, '..', '..', '00_shared', 'lib'))
sys.path.insert(0, LIB); sys.path.insert(0, HERE)
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


def main():
    c = H.build_operators(L); cd = [ci.getH() for ci in c]
    Ham = H.hubbard_hamiltonian(L, t=1.0, U=U, pbc=True, c=c)
    Ntot = sum((cd[m] @ c[m] for m in range(2 * L)), sp.csr_matrix((2 ** (2 * L),) * 2, dtype=complex))
    J = current_operator(L, c)
    for mu_ in [0.0, -1.0, -2.0]:
        cand = eigsh(Ham - mu_ * Ntot, k=1, which='SA')[1][:, 0]
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
            f.write(f"{xs[i]:.6f} {wn[i]:.6f}\n")

    # (a) dense-grid W_n(t) for n=1..4  + (b) dense max_t W_n
    xq = np.linspace(-0.999, 0.999, 400)
    orders = [1, 2, 3, 4]
    Wc = christoffel_widths(xs, wn, orders, xq)
    with open(os.path.join(OUT, 'christoffel_Wn.dat'), 'w') as f:
        f.write("t W1 W2 W3 W4\n")
        for j in range(len(xq)):
            f.write(f"{xq[j]:.5f} " + " ".join(f"{np.clip(Wc[n][j],0,1):.6f}" for n in orders) + "\n")
    with open(os.path.join(OUT, 'christoffel_maxW.dat'), 'w') as f:
        f.write("n maxW\n")
        for n in orders:
            f.write(f"{n} {float(np.nanmax(Wc[n])):.6f}\n")
    print("dense max_t W_n:", {n: round(float(np.nanmax(Wc[n])), 4) for n in orders})

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
