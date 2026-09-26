# -*- coding: utf-8 -*-
"""TWO-SIDED MARKOV-KREIN INTERVAL-MASS BOUND (sim-only, exact measure, high precision).

The paper deploys only the ONE-SIDED Christoffel/Gauss-Radau point bound
    |F_rec(t) - F_true(t)| <= W_n(t),      W_n(t) = 1/sum_{j=0}^n p_j(t)^2,
on the cumulative spectral weight F(t)=mu((-inf,t]).  Classical moment theory
(Chebyshev-Markov-Stieltjes inequalities; Krein-Nudelman, Karlin-Studden,
Dette-Studden) gives a strictly stronger TWO-SIDED bound on the mass
    mu([c,d]) = F(d) - F(c)
that a positive measure sharing m_0..m_{2n} can place on ANY frequency window [c,d]:

    F_-(t) <= F_true(t) <= F_+(t)   for every moment-matching measure,   with
    F_+(t) - F_-(t) = W_n(t)        (the paper's extremal-spread theorem), so
    max(0, F_-(d)-F_+(c)) <= mu([c,d]) <= min(m0, F_+(d)-F_-(c)),
    interval width = W_n(c) + W_n(d).

F_-(t), F_+(t) are the lower/upper principal (Gauss-Radau) representations with a
node fixed at t: F_-(t) = mass of that rule strictly left of t, F_+(t) = F_-(t) +
[its atom at t], and the atom at t equals the Christoffel weight lambda_n(t)=W_n(t)
(Karlin-Studden Ch.II-IV; Dette-Studden Ch.1,3).  This is a rigorous outer bracket;
the SHARP Krein-Nudelman window extremum (attained by a single principal
representation with nodes at both c and d, equivalently the exact-moment optimum of
the moment SDP of Wang et al. / Mortimer) is no wider -- so bracketing is guaranteed.

HONESTY SCOPE (same measure and caveats as chigap_check.py): the L=6, U/t=8 CURRENT
Lehmann measure concentrates ~97% of its weight in a single mid-IR/Mott-band peak at
omega~1.4U; the high-frequency region is an EMPTY ONE-SIDED tail (NOT a two-sided Mott
gap).  Consequently the two-sided interval is TIGHT for any window whose BOTH endpoints
fall in weight-free frequency (out-of-band, or the gaps flanking the peak) and LOOSE only
when a window edge cuts through the dominant peak, where its width equals the same
one-sided max_t W_4 ~ 0.875 the main text already reports.  Not extrapolated to the
two-sided single-particle Mott gap (a different correlator).  NO device, NO fit.

Self-verifying: (V1) a Gauss-Radau node lands exactly on the prescribed t and the rule
reproduces m_0..m_{2n}; (V2) its atom at t equals W_n(t) from the Christoffel recurrence;
(V3) every reported interval brackets the exact window mass.
"""
import os, sys, json
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)   # hubbard_ed.py is vendored in src/
import hubbard_ed as H
from run_sumrule_falsifier import current_operator

mp.mp.dps = 60
L, U = 6, 8.0


def current_measure():
    """Exact L=6, U/t=8 doped-current Lehmann measure: unique poles (rescaled) and residues (m0=1)."""
    c = H.build_operators(L); cd = [ci.getH() for ci in c]
    Ham = H.hubbard_hamiltonian(L, t=1.0, U=U, pbc=True, c=c)
    Ntot = sum((cd[m] @ c[m] for m in range(2 * L)), sp.csr_matrix((2 ** (2 * L),) * 2, dtype=complex))
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
    return np.array(uom), np.array(uw) / np.sum(uw)


# ---------- classical moment machinery (mpmath, high precision) ----------
def recurrence(xs, ws, n):
    """Monic three-term recurrence coeffs a_0..a_n, b_0..b_n (b_0=m0, b_k=||pi_k||^2/||pi_{k-1}||^2)
    from the discrete measure -- these are functions of m_0..m_{2n} only (Hankel PD)."""
    xs = [mp.mpf(float(v)) for v in xs]; ws = [mp.mpf(float(v)) for v in ws]
    N = len(xs)
    pkm1 = [mp.mpf(0)] * N; pk = [mp.mpf(1)] * N
    a = []; b = []; nrm = [mp.fsum(ws)]
    b.append(nrm[0])
    for k in range(n + 1):
        ak = mp.fsum([ws[i] * xs[i] * pk[i] ** 2 for i in range(N)]) / nrm[k]
        a.append(ak)
        if k == n:
            break
        pkp1 = [(xs[i] - ak) * pk[i] - (b[k] if k >= 1 else mp.mpf(0)) * pkm1[i] for i in range(N)]
        nk1 = mp.fsum([ws[i] * pkp1[i] ** 2 for i in range(N)])
        nrm.append(nk1); b.append(nk1 / nrm[k])
        pkm1, pk = pk, pkp1
    return a, b


def _monic_vals(a, b, t, n):
    t = mp.mpf(t); pim1 = mp.mpf(0); pi = mp.mpf(1); vals = [pi]
    for k in range(n):
        pip1 = (t - a[k]) * pi - (b[k] if k >= 1 else mp.mpf(0)) * pim1
        vals.append(pip1); pim1, pi = pi, pip1
    return vals


def christoffel_W(a, b, t, n):
    """W_n(t) = 1 / sum_{j=0}^n p_j(t)^2 (orthonormal p_j)."""
    pis = _monic_vals(a, b, t, n)
    nrm = mp.mpf(1); K = mp.mpf(0)
    for j in range(n + 1):
        nrm = nrm * b[j]            # ||pi_j||^2 = m0 * prod_{i=1}^j b_i
        K += pis[j] ** 2 / nrm
    return 1 / K


def gauss_radau(a, b, m0, t_fix, n):
    """(n+1)-point Gauss-Radau rule with one node fixed at t_fix (Golub/Gautschi):
    modify the last monic diagonal a_n -> a_n* = t_fix - b_n * pi_{n-1}(t)/pi_n(t)."""
    pis = _monic_vals(a, b, t_fix, n)
    an_star = mp.mpf(t_fix) - b[n] * pis[n - 1] / pis[n]
    diag = [a[k] for k in range(n)] + [an_star]
    off = [mp.sqrt(b[k]) for k in range(1, n + 1)]
    m = n + 1
    M = mp.zeros(m)
    for i in range(m):
        M[i, i] = diag[i]
        if i < m - 1:
            M[i, i + 1] = off[i]; M[i + 1, i] = off[i]
    E, V = mp.eigsy(M)
    nodes = [E[i] for i in range(m)]
    wts = [m0 * V[0, i] ** 2 for i in range(m)]
    return nodes, wts


def envelope_F(a, b, m0, t, n):
    """Lower/upper cumulative-weight envelopes at interior node t: F_-(t), F_+(t), and the atom W_n(t)."""
    nodes, wts = gauss_radau(a, b, m0, t, n)
    tf = mp.mpf(t); tol = mp.mpf('1e-9')
    Fminus = mp.mpf(0); atom = mp.mpf(0)
    for nd, wt in zip(nodes, wts):
        if nd < tf - tol:
            Fminus += wt
        elif abs(nd - tf) <= tol:
            atom += wt
    return float(Fminus), float(Fminus + atom), float(atom), (nodes, wts)


def window_bound(a, b, m0, c, d, n):
    """Two-sided Markov-Krein bracket on mu([c,d]) from m_0..m_{2n}."""
    Fmc, Fpc, _, _ = envelope_F(a, b, m0, c, n)
    Fmd, Fpd, _, _ = envelope_F(a, b, m0, d, n)
    lo = max(0.0, Fmd - Fpc)
    hi = min(m0, Fpd - Fmc)
    return lo, hi


def main():
    uom, wn = current_measure()
    lo_e, hi_e = uom.min(), uom.max()
    xs = (2 * uom - (lo_e + hi_e)) / (hi_e - lo_e)      # rescale support to [-1,1]
    m0 = 1.0
    idx = np.argsort(wn)[::-1]
    tpk = float(xs[idx[0]])
    print("=" * 74)
    print(f"TWO-SIDED MARKOV-KREIN INTERVAL-MASS BOUND  (L={L}, U/t={U} current Lehmann)")
    print("=" * 74)
    print(f"{len(uom)} atoms; dominant residue {wn[idx[0]]:.4f} at t={tpk:.4f} (omega={uom[idx[0]]:.3f}), "
          f"next {wn[idx[1]]:.4f} at t={float(xs[idx[1]]):.4f}")
    cum = np.cumsum(wn[np.argsort(xs)])
    print(f"~{100*cum[np.searchsorted(np.sort(xs), -0.35)]:.1f}% of the weight lies below t=-0.35 "
          f"(mid-IR/Mott peak); the region t>0 is an empty one-sided tail\n")

    results = {'_provenance': {'script': 'markov_krein_window.py', 'sim_only': True,
                               'measure': f'L={L} U/t={U} doped-current Lehmann (mid-IR-peak-dominated, one-sided tail)',
                               'n_atoms': int(len(uom)),
                               'bound': 'two-sided Chebyshev-Markov-Stieltjes / Krein-Nudelman envelope bracket',
                               'refs': ['KreinNudelman1977', 'KarlinStudden1966', 'DetteStudden1997']},
               'orders': {}}

    windows = {                              # (c, d) rescaled to [-1,1]
        'across_peak':  (-0.75, -0.35),      # both edges in the gaps flanking the mid-IR peak (contains it)
        'out_of_band':  (0.30, 0.90),        # empty high-frequency tail
    }

    for n in (1, 4):                         # n=1: deployed battery m0,m1,m2 ; n=4: transportable-max m0..m8
        a, b = recurrence(xs, wn, n)
        # ---------- self-verification ----------
        tt = -0.3
        Fm, Fp, atom, (nodes, wts) = envelope_F(a, b, m0, tt, n)
        v1_hit = float(min(abs(float(nd) - tt) for nd in nodes))
        v1_mom = max(abs(float(np.sum(wn * xs ** p)) - float(sum(wt * (nd ** p) for nd, wt in zip(nodes, wts))))
                     for p in range(0, 2 * n + 1))
        v2 = abs(atom - float(christoffel_W(a, b, tt, n)))
        assert v1_hit < 1e-8 and v1_mom < 1e-6 and v2 < 1e-8, "self-verification failed"

        # ---------- worst-case (peak-splitting) width = max_t W_n ----------
        grid = np.linspace(-0.999, 0.999, 400)
        Wg = np.array([float(christoffel_W(a, b, float(t), n)) for t in grid])
        targ = float(grid[int(np.argmax(Wg))]); Wmax = float(Wg.max())
        lo_pk, hi_pk = window_bound(a, b, m0, targ, 0.95, n)

        print(f"--- n={n}  (moments m_0..m_{2 * n}"
              + ("; deployed battery" if n == 1 else "; transportable maximum") + ") ---")
        print(f"  [V1] radau node-hit err={v1_hit:.1e}, moment-repro err(deg<=2n)={v1_mom:.1e}  "
              f"[V2] atom==W_n err={v2:.1e}")
        print(f"  max_t W_n = {Wmax:.4f} at t={targ:.3f};  window edge on the peak -> "
              f"mu in [{lo_pk:.4f}, {hi_pk:.4f}] (width {hi_pk - lo_pk:.4f} = max_t W_n)")
        order_rec = {'moments_used': f'm0..m{2 * n}', 'max_t_Wn': Wmax,
                     'peak_edge_window': {'edges': [targ, 0.95], 'interval': [lo_pk, hi_pk],
                                          'width': hi_pk - lo_pk}, 'windows': {}}
        for name, (c, d) in windows.items():
            lo, hi = window_bound(a, b, m0, c, d, n)
            true_mass = float(np.sum(wn[(xs >= c) & (xs <= d)]))
            Wc = float(christoffel_W(a, b, c, n)); Wd = float(christoffel_W(a, b, d, n))
            brackets = (lo - 1e-9) <= true_mass <= (hi + 1e-9)
            assert brackets, f"bracket failed on {name}"
            tag = 'TIGHT' if (hi - lo) < 0.05 else 'loose'
            print(f"  window {name:12s} [{c:+.2f},{d:+.2f}]: mu in [{lo:.4f}, {hi:.4f}]  "
                  f"width {hi - lo:.4f}  true {true_mass:.4f}  ({tag}; W_n(c)+W_n(d)={Wc + Wd:.4f})")
            order_rec['windows'][name] = {'edges': [c, d], 'interval': [lo, hi], 'width': hi - lo,
                                          'true_mass': true_mass, 'brackets': bool(brackets)}
        results['orders'][f'n={n}'] = order_rec
        print()

    outdir = os.path.normpath(os.path.join(HERE, '..', 'data'))
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, '2026-08-28_markov_krein_window.json')
    json.dump(results, open(outpath, 'w'), indent=2)
    print("all self-verifications passed; wrote", os.path.relpath(outpath, os.path.join(HERE, '..')))


if __name__ == '__main__':
    main()
