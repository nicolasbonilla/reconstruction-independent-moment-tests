# -*- coding: utf-8 -*-
r"""Sector-basis Lanczos + Haydock continued-fraction spectral engine for the 1D Hubbard ring (PBC),
adapted from the validated Paper-2 akw_lanczos.py and EXTENDED to the dynamical structure factors.
Scales far past dense ED: ground state by sparse Lanczos in the (up-string x dn-string) sector, each
Green's function by ~nl Haydock steps of sparse mat-vecs. Matrix-free ground state for large L.

  A(k,w)      = -1/pi Im[ <0|c_k (z-(H-E0))^{-1} c_k^dag|0> + <0|c_k^dag (z-(E0-H))^{-1} c_k|0> ]
  S(q,w)      = -1/pi Im  <0| rho_q^dag (z-(H-E0))^{-1} rho_q |0>     (charge, N-conserving)
  S^zz(q,w)   = -1/pi Im  <0| Sz_q^dag  (z-(H-E0))^{-1} Sz_q  |0>     (spin,   N-conserving)
Validated against exact sector diagonalization at L=6 (Paper 2)."""
import os, sys, time, numpy as np, scipy.sparse as sp
from itertools import combinations
from scipy.sparse.linalg import eigsh, LinearOperator
t0 = time.time(); log = lambda *a: print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')


def strings(L, n):
    S = []
    for c in combinations(range(L), n):
        m = 0
        for b in c:
            m |= (1 << b)
        S.append(m)
    S.sort(); idx = {m: i for i, m in enumerate(S)}
    return S, idx


def hop(L, n, t=1.0):
    """single-spin hopping -t sum_<ij>(c^dag_i c_j + h.c.), PBC, in the n-string basis."""
    S, idx = strings(L, n); D = len(S); rows = []; cols = []; val = []
    for a, m in enumerate(S):
        for i in range(L):
            j = (i + 1) % L
            for (p, q) in ((i, j), (j, i)):
                if (m >> q) & 1 and not (m >> p) & 1:
                    m2 = (m & ~(1 << q)) | (1 << p)
                    lo, hi = min(p, q), max(p, q)
                    mask = m & (((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1))
                    sign = -1.0 if (bin(mask).count('1') & 1) else 1.0
                    rows.append(idx[m2]); cols.append(a); val.append(-t * sign)
    return sp.csr_matrix((val, (rows, cols)), shape=(D, D)), S, idx


def cdag_map(L, n):
    """c^dag_j : n-string -> (n+1)-string, per site j (list of sparse matrices)."""
    Sa, ia = strings(L, n); Sb, ib = strings(L, n + 1)
    ops = []
    for j in range(L):
        rows = []; cols = []; val = []
        for a, m in enumerate(Sa):
            if not (m >> j) & 1:
                m2 = m | (1 << j)
                sign = -1.0 if (bin(m & ((1 << j) - 1)).count('1') & 1) else 1.0
                rows.append(ib[m2]); cols.append(a); val.append(sign)
        ops.append(sp.csr_matrix((val, (rows, cols)), shape=(len(Sb), len(Sa))))
    return ops


def occ_matrix(S, L):
    return np.array([[(m >> i) & 1 for i in range(L)] for m in S], dtype=float)


def sector_H(L, U, nup, ndn, t=1.0):
    """matrix-free H for the (nup,ndn) sector: Tu(x)I + I(x)Td + U*doublons."""
    Tu, Su, iu = hop(L, nup, t); Td, Sd, idd = hop(L, ndn, t)
    Du = len(Su); Dd = len(Sd)
    upocc = occ_matrix(Su, L); dnocc = occ_matrix(Sd, L)
    diagU = U * (upocc @ dnocc.T).ravel()

    def matvec(x):
        X = x.reshape(Du, Dd)
        Y = Tu @ X + (Td @ X.T).T
        return Y.ravel() + diagU * x
    H = LinearOperator((Du * Dd, Du * Dd), matvec=matvec, dtype=float)
    return H, Su, iu, Sd, idd, Du, Dd


def ground_state(L, U, nup, nd, t=1.0, explicit_max=12):
    """GS energy+vector of the (nup,nd) sector. Explicit sparse H for small L (fast ARPACK);
    matrix-free LinearOperator for large L (keeps RAM bounded, no multi-GB kron matrix)."""
    Hlin, Su, iu, Sd, idd, Du, Dd = sector_H(L, U, nup, nd, t)
    if L <= explicit_max:
        Tu, _, _ = hop(L, nup, t); Td, _, _ = hop(L, nd, t)
        upocc = occ_matrix(Su, L); dnocc = occ_matrix(Sd, L)
        diagU = U * (upocc @ dnocc.T).ravel()
        H = (sp.kron(Tu, sp.identity(Dd), format='csr')
             + sp.kron(sp.identity(Du), Td, format='csr') + sp.diags(diagU)).tocsr()
        E, V = eigsh(H, k=1, which='SA', v0=np.random.default_rng(0).standard_normal(H.shape[0]))
    else:
        E, V = eigsh(Hlin, k=1, which='SA', ncv=10, maxiter=2000,
                     v0=np.random.default_rng(0).standard_normal(Hlin.shape[0]))
    return float(E[0]), V[:, 0], Du, Dd


def current_string(L, n, t=1.0):
    """single-spin current J_s = -i t sum_i (c^dag_{i+1} c_i - c^dag_i c_{i+1}), PBC, in the n-string basis."""
    S, idx = strings(L, n); D = len(S); rows = []; cols = []; val = []
    for a, m in enumerate(S):
        for i in range(L):
            j = (i + 1) % L
            for (p, q, coeff) in ((j, i, -1j * t), (i, j, +1j * t)):   # c^dag_p c_q
                if (m >> q) & 1 and not (m >> p) & 1:
                    m2 = (m & ~(1 << q)) | (1 << p)
                    lo, hi = min(p, q), max(p, q)
                    mask = m & (((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1))
                    sign = -1.0 if (bin(mask).count('1') & 1) else 1.0
                    rows.append(idx[m2]); cols.append(a); val.append(coeff * sign)
    return sp.csr_matrix((val, (rows, cols)), shape=(D, D))


def haydock_poles(Hfun, v0, nl):
    """Ritz (Lanczos) pole representation of the spectral density of v0 under H:
    returns (theta_j, w_j) with sum w_j = <v0|v0> and sum theta_j w_j = <v0|H|v0> (exact to order 2*nl)."""
    nrm2 = np.vdot(v0, v0).real
    if nrm2 < 1e-14:
        return np.array([0.0]), np.array([0.0])
    v = v0 / np.sqrt(nrm2); vp = np.zeros_like(v); b = 0.0; a = []; bs = []
    for _ in range(nl):
        w = Hfun(v); an = np.vdot(v, w).real; a.append(an)
        w = w - an * v - b * vp
        bn = np.sqrt(np.vdot(w, w).real)
        if bn < 1e-10:
            break
        bs.append(bn); vp = v; v = w / bn; b = bn
    a = np.array(a); bs = np.array(bs[:len(a) - 1])
    T = np.diag(a) + np.diag(bs, 1) + np.diag(bs, -1)
    theta, Sv = np.linalg.eigh(T)
    return theta, nrm2 * (Sv[0, :] ** 2)


def run_current(L=12, U=8.0, t=1.0, nl=220, nup=None, nd=None, explicit_max=12):
    """current spectral function A_J(omega) poles for a DOPED Hubbard ring (fixed sector = fixed filling)."""
    if nup is None:
        nup = nd = (L // 3)                       # ~2/3 filling (matches the L=6 doped illustration)
    E0, psi0, Du, Dd = ground_state(L, U, nup, nd, t, explicit_max)
    log(f"[cur] L={L} U={U} filling={2*nup}/{L} per spin nup=nd={nup}: GS dim={Du*Dd}  E0={E0:.6f}")
    Ju = current_string(L, nup, t); Jd = current_string(L, nd, t)
    Psi = psi0.reshape(Du, Dd)
    JPsi = (Ju @ Psi) + (Psi @ Jd.T)              # J|0> = (Ju (x) I + I (x) Jd)|0>, stays in the sector
    v0 = JPsi.ravel()
    Hfun = sector_H(L, U, nup, nd, t)[0]
    om, w = haydock_poles(lambda x: Hfun.matvec(x) - E0 * x, v0, nl)
    m0 = float(w.sum()); m1 = float((om * w).sum())
    keep = w > 1e-10; om, w = om[keep], w[keep]
    os.makedirs(CACHE, exist_ok=True)
    np.savez(os.path.join(CACHE, f'current_L{L}.npz'), L=L, U=U, nup=nup, nd=nd, om=om, w=w, m0=m0, m1=m1, E0=E0)
    log(f"[cur] m0=<J^2>={m0:.4f}  m1={m1:.4f}  poles={len(om)}  om-range=[{om.min():.2f},{om.max():.2f}]")
    return om, w, m0, m1


def _sub_gs(Hsub):
    """ground state (E0, vec) of a subspace matrix -- dense for tiny d, sparse Lanczos otherwise."""
    d = Hsub.shape[0]
    if d < 24:
        e, V = np.linalg.eigh(Hsub.toarray())
        return float(e[0]), V[:, 0]
    e, V = eigsh(Hsub.tocsr(), k=1, which='SA',
                 v0=np.random.default_rng(0).standard_normal(d))   # fixed start -> deterministic
    return float(e[0]), V[:, 0]


def _logspace_d(nsup, npts, dmin):
    lo, hi = np.log10(max(dmin, 4)), np.log10(nsup)
    ds = np.unique(np.round(10 ** np.linspace(hi, lo, npts)).astype(int))[::-1]
    return np.clip(ds, 4, nsup)


def run_teeth(L=12, U=8.0, t=1.0, npts=18, dmin=None, cov_min=0.15):
    """THE TEETH EXPERIMENT in the sector basis: subspace-truncation error vs the INDEPENDENT ground-state
    first moment m1_op. For each coverage f, keep the top-d determinants by |psi0|^2 (SQD ordering), and get
    the truncated first moment WITHOUT full diagonalization: m1_trunc=<J0_d|(H_d-E0_d)|J0_d> via one Lanczos
    ground state of the d-dim subspace. Scales to L=12 (245k determinants). Circular control is 0 by construction."""
    nup = nd = L // 3                                  # ~2/3 filling (matches fig_falsifier)
    Tu, Su, iu = hop(L, nup, t); Td, Sd, idd = hop(L, nd, t)
    Du = len(Su); Dd = len(Sd)
    upocc = occ_matrix(Su, L); dnocc = occ_matrix(Sd, L)
    diagU = U * (upocc @ dnocc.T).ravel()
    Hs = (sp.kron(Tu, sp.identity(Dd), format='csr') + sp.kron(sp.identity(Du), Td, format='csr')
          + sp.diags(diagU)).tocsr()
    Ju = current_string(L, nup, t); Jd = current_string(L, nd, t)
    Js = (sp.kron(Ju, sp.identity(Dd), format='csr') + sp.kron(sp.identity(Du), Jd, format='csr')).tocsr()
    log(f"[teeth] L={L} U={U} filling={2*nup}/{L}: sector dim={Du*Dd}, H nnz={Hs.nnz}")
    E0, psi0 = _sub_gs(Hs)
    Jp = Js @ psi0
    m1_op = float(np.real(np.vdot(Jp, Hs @ Jp)) - E0 * np.real(np.vdot(Jp, Jp)))   # independent estimator
    prob = np.abs(psi0) ** 2
    order = np.argsort(prob)[::-1]
    nsup = int((prob > 1e-14).sum()); order = order[:nsup]
    log(f"[teeth] E0={E0:.6f}  m1_op={m1_op:.6f}  n_support={nsup}")
    if dmin is None:
        dmin = max(30, int(round(cov_min * nsup)))
    ds = _logspace_d(nsup, npts, dmin)                 # log-spaced determinant counts (fair across L)
    rows = []
    for d in ds:
        d = int(d); idx = np.sort(order[:d])
        Hsub = Hs[idx][:, idx]; Jsub = Js[idx][:, idx]
        e0, g = _sub_gs(Hsub)
        Jg = Jsub @ g
        m1_trunc = float(np.real(np.vdot(Jg, Hsub @ Jg)) - e0 * np.real(np.vdot(Jg, Jg)))
        r_ind = abs(m1_trunc - m1_op) / abs(m1_op)
        rows.append((d / nsup, d, r_ind, 0.0))
        log(f"[teeth]   d={d:6d} (cov={100*d/nsup:5.1f}%)  r_indep={100*r_ind:7.2f}%")
    os.makedirs(CACHE, exist_ok=True)
    np.savez(os.path.join(CACHE, f'teeth_L{L}.npz'), L=L, U=U, m1_op=m1_op, nsup=nsup,
             cov=np.array([r[0] for r in rows]), d=np.array([r[1] for r in rows]),
             rind=np.array([r[2] for r in rows]), rcirc=np.array([r[3] for r in rows]))
    log(f"[teeth] WROTE cache/teeth_L{L}.npz")
    return rows, m1_op, nsup


def haydock(Hfun, v0, nl, z):
    """continued fraction <v0|(z-H)^{-1}|v0> via nl Lanczos steps (v0 need not be normalized)."""
    nrm2 = np.vdot(v0, v0).real
    if nrm2 < 1e-14:
        return np.zeros_like(z)
    v = v0 / np.sqrt(nrm2); vp = np.zeros_like(v); b = 0.0; a = []; bs = []
    for _ in range(nl):
        w = Hfun(v); an = np.vdot(v, w).real; a.append(an)
        w = w - an * v - b * vp
        bn = np.sqrt(np.vdot(w, w).real); bs.append(bn)
        if bn < 1e-10:
            break
        vp = v; v = w / bn; b = bn
    a = np.array(a); bs = np.array(bs)
    G = z - a[-1]
    for m in range(len(a) - 2, -1, -1):
        G = (z - a[m]) - bs[m] ** 2 / G
    return nrm2 / G


def run_akw(L=12, U=8.0, t=1.0, eta=0.15, nl=180, nw=700, wmax=9.0, explicit_max=12):
    nup = nd = L // 2
    E0, psi0, Du, Dd = ground_state(L, U, nup, nd, t, explicit_max)
    Psi = psi0.reshape(Du, Dd)
    log(f"[akw] L={L} U={U}: GS dim={Du*Dd}  E0={E0:.6f}")
    # single-particle Mott gap from the N+/-1 ground sectors -> anchored band edges mu+/-
    Ep, _, _, _ = ground_state(L, U, nup + 1, nd, t, explicit_max)
    Em, _, _, _ = ground_state(L, U, nup - 1, nd, t, explicit_max)
    mu_plus = Ep - E0; mu_minus = E0 - Em; gap = mu_plus - mu_minus
    mu = 0.5 * (mu_plus + mu_minus)                       # E_F at mid-gap (= U/2 by PH symmetry)
    log(f"[akw] mu+={mu_plus:.4f} mu-={mu_minus:.4f} gap={gap:.4f}  E_F=mid-gap")
    Ha = sector_H(L, U, nup + 1, nd, t)[0]
    Hr = sector_H(L, U, nup - 1, nd, t)[0]
    cdU = cdag_map(L, nup)                                 # c^dag_j : nup -> nup+1
    cU = cdag_map(L, nup - 1)                              # c_j = cU[j].T : nup -> nup-1
    wg = np.linspace(-wmax, wmax, nw)                      # omega measured from E_F
    zt = (wg + mu) + 1j * eta
    ks = list(range(L))
    A = np.zeros((len(ks), nw))
    for n in ks:
        k = 2 * np.pi * n / L
        ph = np.exp(1j * k * np.arange(L)) / np.sqrt(L)
        Dua = cdU[0].shape[0]
        add = np.zeros((Dua, Dd), dtype=complex)
        for j in range(L):
            add += ph[j] * (cdU[j] @ Psi)
        Gp = haydock(lambda x: Ha.matvec(x) - E0 * x, add.ravel(), nl, zt)
        Dur = cU[0].shape[1]
        rem = np.zeros((Dur, Dd), dtype=complex)
        for j in range(L):
            rem += np.conj(ph[j]) * (cU[j].T @ Psi)
        Gm = haydock(lambda x: E0 * x - Hr.matvec(x), rem.ravel(), nl, zt)
        A[n] = -(1.0 / np.pi) * np.imag(Gp + Gm)
        log(f"[akw]   k/pi={2*n/L:.3f}  int A dw={np.trapz(A[n], wg):.3f}")
    kk = np.array(ks) / L * 2.0
    os.makedirs(CACHE, exist_ok=True)
    np.savez(os.path.join(CACHE, f'akw_L{L}.npz'), L=L, U=U, eta=eta, kk=kk, wg=wg, A=A,
             E0=E0, mu_plus=mu_plus, mu_minus=mu_minus, gap=gap)
    log(f"[akw] WROTE cache/akw_L{L}.npz  ({len(ks)} k-points, nw={nw})")
    return kk, wg, A, gap


def struct_factors(L=12, U=8.0, t=1.0, eta_c=0.20, eta_s=0.12, nl=180,
                   nwc=700, nws=700, wc_max=11.0, ws_max=2.6, explicit_max=12):
    """charge S(q,w) and spin S^zz(q,w): N-conserving Haydock in the ground sector."""
    nup = nd = L // 2
    E0, psi0, Du, Dd = ground_state(L, U, nup, nd, t, explicit_max)
    Psi = psi0.reshape(Du, Dd)
    log(f"[sqw] L={L} U={U}: GS dim={Du*Dd}  E0={E0:.6f}")
    Hgs = sector_H(L, U, nup, nd, t)[0]
    Hfun = lambda x: Hgs.matvec(x) - E0 * x
    # site-basis number and Sz occupation matrices (Du,L) and (Dd,L)
    Su, _ = strings(L, nup); Sd, _ = strings(L, nd)
    upocc = occ_matrix(Su, L); dnocc = occ_matrix(Sd, L)     # (Du,L),(Dd,L)
    wc = np.linspace(0, wc_max, nwc); zc = wc + 1j * eta_c
    ws = np.linspace(0, ws_max, nws); zs = ws + 1j * eta_s
    qs = list(range(1, L))                                   # q=2pi m/L, m=1..L-1
    Sc = np.zeros((len(qs), nwc)); Ss = np.zeros((len(qs), nws))
    for iq, m in enumerate(qs):
        q = 2 * np.pi * m / L
        phj = np.exp(-1j * q * np.arange(L))                 # (L,)
        # rho_q|0> and Sz_q|0> acting site-wise: n_j = nup_j (x) I + I (x) ndn_j ; Sz_j = 1/2(nup_j - ndn_j)
        # density: (rho_q Psi)_{ab} = sum_j phj[j]*(upocc[a,j]+dnocc[b,j]) Psi_{ab}
        rq = (upocc @ phj)[:, None] * Psi + Psi * (dnocc @ phj)[None, :]
        sq = 0.5 * ((upocc @ phj)[:, None] * Psi - Psi * (dnocc @ phj)[None, :])
        # subtract ground-state expectation (q!=0 so <rho_q>=<Sz_q>=0, but keep for safety)
        Sc[iq] = -(1.0 / np.pi) * np.imag(haydock(Hfun, rq.ravel(), nl, zc))
        Ss[iq] = -(1.0 / np.pi) * np.imag(haydock(Hfun, sq.ravel(), nl, zs))
        log(f"[sqw]   q/pi={2*m/L:.3f}  int S_c={np.trapz(Sc[iq], wc):.3f}  int S_zz={np.trapz(Ss[iq], ws):.3f}")
    qq = np.array(qs) / L * 2.0
    os.makedirs(CACHE, exist_ok=True)
    np.savez(os.path.join(CACHE, f'sqw_L{L}.npz'), L=L, U=U, qq=qq, wc=wc, ws=ws, Sc=Sc, Ss=Ss,
             eta_c=eta_c, eta_s=eta_s, wc_max=wc_max, ws_max=ws_max)
    log(f"[sqw] WROTE cache/sqw_L{L}.npz")
    return qq, wc, ws, Sc, Ss


if __name__ == "__main__":
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    what = sys.argv[2] if len(sys.argv) > 2 else 'both'
    if what in ('akw', 'both'):
        run_akw(L=L, U=8.0, eta=0.16, nl=180, nw=700, wmax=9.0)
    if what in ('sqw', 'both'):
        struct_factors(L=L, U=8.0, eta_c=0.20, eta_s=0.12, nl=180)
    log("DONE")
