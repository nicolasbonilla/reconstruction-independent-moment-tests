# -*- coding: utf-8 -*-
"""Export native-figure .dat for NEW fig_bracketing: the mechanism of Theorem 1. Finitely many moments
do not pin a spectrum -- they BRACKET its cumulative weight F(t) between the extremal (Gauss-Radau
principal-representation) measures, and the bracket width IS the Christoffel function W_n(t), which
contracts as n grows. Measure = doped Hubbard L=6 U=8 current-operator Lehmann density (same as fig_christoffel)."""
import os, sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.normpath(os.path.join(HERE, '..', '..', '00_shared', 'lib'))
sys.path.insert(0, LIB); sys.path.insert(0, HERE)
import hubbard_ed as H
from run_sumrule_falsifier import current_operator
OUT = os.path.normpath(os.path.join(HERE, '..', 'paper', 'figs'))
L, U = 6, 8.0


def recurrence(x, w, K):
    """Monic three-term recurrence coeffs from the discrete measure {(x_i,w_i)}.
    Returns alpha[0..K], beta[0..K] (beta[0]=m0, beta[k]=norm2[k]/norm2[k-1]), and norm2[0..K]."""
    n2 = [float(np.sum(w))]
    alpha = []; beta = [n2[0]]
    pkm1 = np.zeros_like(x); pk = np.ones_like(x)
    for k in range(K + 1):
        ak = float(np.sum(w * x * pk ** 2) / n2[k]); alpha.append(ak)
        bk = (n2[k] / n2[k - 1]) if k >= 1 else 0.0
        if k >= 1:
            beta.append(n2[k] / n2[k - 1])
        pkp1 = (x - ak) * pk - (bk if k >= 1 else 0.0) * pkm1
        if k < K:
            n2.append(float(np.sum(w * pkp1 ** 2)))
        pkm1, pk = pk, pkp1
    return np.array(alpha), np.array(beta), np.array(n2)


def monic_eval(alpha, beta, t, n):
    """Evaluate monic pi_{n-1}(t), pi_n(t)."""
    pm1 = 0.0; p0 = 1.0
    for k in range(n):
        p1 = (t - alpha[k]) * p0 - (beta[k] if k >= 1 else 0.0) * pm1
        pm1, p0 = p0, p1
    return pm1, p0     # pi_{n-1}(t), pi_n(t)


def radau_bracket(alpha, beta, n, tgrid):
    """Gauss-Radau (n+1 nodes, one fixed at t) principal representations of F(t).
    Returns Flo(t)=sum weights strictly below t, Fhi(t)=Flo+weight_at_t (=W_n(t))."""
    Flo = np.zeros_like(tgrid); Fhi = np.zeros_like(tgrid)
    m0 = beta[0]
    for it, t in enumerate(tgrid):
        pnm1, pn = monic_eval(alpha, beta, t, n)
        if abs(pn) < 1e-300:
            Flo[it] = np.nan; Fhi[it] = np.nan; continue
        aR = t - beta[n] * pnm1 / pn                        # modified last diagonal
        d = np.array([alpha[k] for k in range(n)] + [aR])   # size n+1
        e = np.sqrt(np.array([beta[k] for k in range(1, n + 1)]))
        Tm = np.diag(d) + np.diag(e, 1) + np.diag(e, -1)
        ev, V = np.linalg.eigh(Tm)
        wts = m0 * (V[0, :] ** 2)
        j = int(np.argmin(np.abs(ev - t)))                  # the node fixed at t
        w_at_t = wts[j]
        Flo[it] = float(np.sum(wts[ev < t - 1e-9]))
        Fhi[it] = Flo[it] + w_at_t
    return Flo, Fhi


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
    # aggregate degenerate frequencies
    order = np.argsort(om); om, w = om[order], w[order]
    uom, uw = [], []
    for o, wi in zip(om, w):
        if uom and abs(o - uom[-1]) < 1e-6:
            uw[-1] += wi
        else:
            uom.append(o); uw.append(wi)
    uom = np.array(uom); uw = np.array(uw); wn = uw / uw.sum()
    a, b = uom.min(), uom.max()
    xs = (2 * uom - (a + b)) / (b - a)

    # true cumulative F(t) staircase
    o2 = np.argsort(xs); xss = xs[o2]; wss = wn[o2]
    Fcum = np.cumsum(wss)
    tgridF = np.linspace(-1, 1, 600)
    Ftrue = np.array([wss[xss <= t].sum() for t in tgridF])
    with open(os.path.join(OUT, 'bracket_F.dat'), 'w') as f:
        f.write("t F\n"); [f.write(f"{t:.4f} {v:.6f}\n") for t, v in zip(tgridF, Ftrue)]

    alpha, beta, n2 = recurrence(xs, wn, 6)
    tq = np.linspace(-0.985, 0.985, 240)
    for n, name in [(2, 'bracket_env_n2'), (4, 'bracket_env_n4')]:
        Flo, Fhi = radau_bracket(alpha, beta, n, tq)
        with open(os.path.join(OUT, name + '.dat'), 'w') as f:
            f.write("t Flo Fhi\n")
            for t, lo, hi in zip(tq, Flo, Fhi):
                if np.isfinite(lo):
                    f.write(f"{t:.4f} {max(lo,0):.6f} {min(hi,1):.6f}\n")
        width = np.nanmax(Fhi - Flo)
        print(f"n={n}: max bracket width (=max W_n) = {width:.4f}")
    # max W_n vs n for the inset
    with open(os.path.join(OUT, 'bracket_maxW.dat'), 'w') as f:
        f.write("n maxW\n")
        for n in [1, 2, 3, 4]:
            Flo, Fhi = radau_bracket(alpha, beta, n, tq)
            f.write(f"{n} {float(np.nanmax(Fhi - Flo)):.6f}\n")
    print("wrote bracket_*.dat ->", OUT)


if __name__ == '__main__':
    main()
