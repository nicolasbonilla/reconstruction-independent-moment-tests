# -*- coding: utf-8 -*-
r"""
BEYOND-ED DMRG CHECK -- the shared-state moment test on a finite-sampling SQD subspace at sector
dimensions where full exact diagonalization is infeasible, with DMRG/MPS as the independent reference.
================================================================================================

Model: 1D Fermi-Hubbard chain, OPEN boundary, U/t given, half filling (nup = ndn = L/2, S_z=0).
Probe O = total current  J = -i t sum_i (c^dag_{i+1} c_i - h.c.).

Roles (kept separate):
  * DMRG (TeNPy) role (a): the ground state |0> from which the subspace S is SAMPLED (Born-rule
    perfect sampling of occupation configurations ~ |<config|0>|^2 at a finite shot budget).
  * DMRG (TeNPy) role (b): the independent reference.  Moments of the same probe,
    M_k = <0| J (H-E0)^k J |0>, from the current MPO applied to the MPS.  DMRG is never used in the test.

The test (classical post-processing of the SAMPLED subspace with FULL sparse operators, as in
run_teeth_shared.py / Eq.(app-leak)):
    |0_S>  = ground state of  H_S = Pi_S H Pi_S     (subspace Hamiltonian on the sampled configs)
    e0     = <0_S| H |0_S>                          (shared energy reference)
    m_bar_k = <0_S| J_S (H_S - e0)^k J_S |0_S>      reconstruction moment (within-S Lehmann)
    m_hat_k = <0_S| J   (H   - e0)^k J   |0_S>      full-operator moment (J, H reach outside S)
    Delta_k = |m_hat_k - m_bar_k|                   (Delta_0 is the weight of J|0_S> outside S)

Degeneracy: if the lowest level of H_S is degenerate (gap < DEGEN_TOL), |0_S> is not unique and
every moment above depends on the chosen vector.  The script then reports [min, max] of each quantity
over the ground manifold and leaves such rows out of the monotonicity and correlation summaries.

m_hat is computed with on-the-fly configuration-basis operators acting on the sampled-subspace ground
state (no diagonalization of the full sector, no DMRG).  The configuration-basis operator conventions
are checked in-run against full ED at L=8.

Scope: 1D Hubbard is DMRG-solvable, so this checks the test code at sector dimensions beyond full ED
(L=20: 3.4e10, L=24: 7.3e12) against DMRG references; it is not a beyond-classical claim.

Run:  python src/beyond_ed_dmrg_catch.py   (needs physics-tenpy; set BEYOND_ED_CACHE to relocate the
DMRG/sample cache, default src/cache) -> data/<DATE>_beyond_ed_dmrg_catch.json
"""
import os, sys, json, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from itertools import combinations
from collections import Counter

t0 = time.time()
log = lambda *a: print(f"[{time.time()-t0:8.1f}s]", *a, flush=True)
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, '..', 'data'))
CACHE = os.environ.get('BEYOND_ED_CACHE', os.path.join(HERE, 'cache'))   # DMRG/sample cache
DATE = '2026-10-02'


# ============================ vectorized config-basis operators ============================
def popcount64(x):
    """population count of a non-negative int64 numpy array."""
    x = x.astype(np.int64)
    x = x - ((x >> 1) & 0x5555555555555555)
    x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F
    return (x * 0x0101010101010101) >> 56


def between_mask(p, q):
    lo, hi = (p, q) if p < q else (q, p)
    return ((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1)


def _pack(up, dn, L):
    return (up.astype(np.int64) << L) | dn.astype(np.int64)


def _sig(s, betw):
    """(-1)^(popcount(s & betw)) as +-1 float, vectorized."""
    return 1.0 - 2.0 * (popcount64(s & betw) & 1)


def offdiag_terms(up, dn, L, t, kind):
    """All off-diagonal moves of H-hopping (kind='hop') or the current (kind='cur') acting on the
    config set (up[],dn[]).  Returns (src_idx, new_up, new_dn, amp) concatenated over spins/bonds.
    Up-string and dn-string carry independent JW signs (kron(Tu,I)+kron(I,Td) convention, verified
    against TeNPy to machine precision)."""
    n = len(up); ar = np.arange(n)
    SI = []; NU = []; ND = []; VAL = []
    for i in range(L - 1):
        j = i + 1
        if kind == 'hop':
            dirs = ((i, j, -t + 0j), (j, i, -t + 0j))       # c^dag_p c_q, p=first
        else:
            dirs = ((j, i, -1j * t), (i, j, +1j * t))       # J = -it(c^dag_{i+1}c_i) + h.c.
        for (p, q, coeff) in dirs:
            betw = between_mask(p, q); bp = 1 << p; bq = 1 << q
            # up spin
            m = (((up >> q) & 1) == 1) & (((up >> p) & 1) == 0)
            idx = ar[m]
            if idx.size:
                nu = (up[idx] & ~bq) | bp
                amp = coeff * _sig(up[idx], betw)
                SI.append(idx); NU.append(nu); ND.append(dn[idx]); VAL.append(amp)
            # dn spin
            m = (((dn >> q) & 1) == 1) & (((dn >> p) & 1) == 0)
            idx = ar[m]
            if idx.size:
                nd = (dn[idx] & ~bq) | bp
                amp = coeff * _sig(dn[idx], betw)
                SI.append(idx); NU.append(up[idx]); ND.append(nd); VAL.append(amp)
    if not SI:
        z = np.array([], dtype=np.int64)
        return z, z, z, np.array([], dtype=complex)
    return (np.concatenate(SI), np.concatenate(NU), np.concatenate(ND),
            np.concatenate(VAL).astype(complex))


def build_op_on_set(up, dn, keys_sorted, order, L, t, U, kind):
    """Sparse operator (H if kind='hop'+U diagonal, J if kind='cur') PROJECTED onto the config set
    whose packed keys are keys_sorted (sorted) with `order` mapping sorted->original index."""
    n = len(up)
    si, nu, nd, val = offdiag_terms(up, dn, L, t, kind)
    rows = []; cols = []; vals = []
    if si.size:
        nk = _pack(nu, nd, L)
        pos = np.searchsorted(keys_sorted, nk)
        pos = np.clip(pos, 0, len(keys_sorted) - 1)
        hit = keys_sorted[pos] == nk
        dst = order[pos[hit]]
        rows.append(dst); cols.append(si[hit]); vals.append(val[hit])
    if kind == 'hop' and U != 0.0:
        diag = U * popcount64(up & dn).astype(float)
        rows.append(np.arange(n)); cols.append(np.arange(n)); vals.append(diag.astype(complex))
    if rows:
        r = np.concatenate(rows); c = np.concatenate(cols); v = np.concatenate(vals)
    else:
        r = c = np.array([], dtype=int); v = np.array([], dtype=complex)
    return sp.csr_matrix((v, (r, c)), shape=(n, n))


def apply_current_vec(up, dn, amp, L, t):
    """J acting on a sparse vector supported on (up,dn) with amplitudes amp -> aggregated
    (out_up, out_dn, out_amp) over the WHOLE Hilbert space (reaches configs outside the set)."""
    si, nu, nd, val = offdiag_terms(up, dn, L, t, 'cur')
    contrib = val * amp[si]
    keys = _pack(nu, nd, L)
    uk, inv = np.unique(keys, return_inverse=True)
    out = np.zeros(len(uk), dtype=complex)
    np.add.at(out, inv, contrib)
    out_up = (uk >> L).astype(np.int64)
    out_dn = (uk & ((1 << L) - 1)).astype(np.int64)
    return out_up, out_dn, out, uk


# ================================ the moment screen on S ================================
DEGEN_TOL = 1e-8      # two subspace eigenvalues closer than this are treated as one degenerate level
MAX_MANIFOLD = 6      # cap on the number of degenerate ground states collected by deflation


def _lowest(Hop, n, rng):
    """Lowest eigenpair of a Hermitian operator/LinearOperator of size n (ARPACK, random start)."""
    v0 = rng.standard_normal(n) + 0j
    ev, V = eigsh(Hop, k=1, which='SA', v0=v0, tol=0, maxiter=20000)
    return float(ev[0].real), V[:, 0]


def subspace_ground_manifold(H_S, seed=12345):
    """Ground level of H_S with an explicit degeneracy check.

    Eigenpairs are found one at a time: after each one, the found vectors are shifted up by deflation
    (H_S + sigma * sum_j g_j g_j^H) and the lowest eigenpair of the deflated operator is computed from a
    fresh random start.  This finds a degenerate partner even when a single-vector Lanczos run would
    see only one vector of the level.  Collection stops at the first eigenvalue more than DEGEN_TOL
    above the lowest one; its distance is returned as the gap.
    Returns (evals_in_manifold, G [n x d orthonormal], gap_to_next)."""
    from scipy.sparse.linalg import LinearOperator
    n = H_S.shape[0]
    if n < 40:
        ev, V = np.linalg.eigh(H_S.toarray())
        ev = ev.real
        d = int(np.sum(ev - ev[0] < DEGEN_TOL))
        gap = float(ev[d] - ev[0]) if d < n else float('inf')
        return ev[:d].copy(), V[:, :d].astype(complex), gap
    rng = np.random.default_rng(seed)
    e, g = _lowest(H_S, n, rng)
    evs = [e]; G = [g / np.linalg.norm(g)]
    sigma = 10.0 * (abs(e) + float(abs(H_S).sum(axis=1).max()) + 1.0)
    gap = float('inf')
    while len(G) < min(MAX_MANIFOLD + 1, n - 1):
        Gm = np.column_stack(G)
        op = LinearOperator((n, n), dtype=complex,
                            matvec=lambda x, Gm=Gm: H_S @ x + sigma * (Gm @ (Gm.conj().T @ x)))
        e2, g2 = _lowest(op, n, rng)
        g2 = g2 - Gm @ (Gm.conj().T @ g2); g2 = g2 / np.linalg.norm(g2)
        e2 = float(np.real(np.vdot(g2, H_S @ g2)))
        if e2 - evs[0] < DEGEN_TOL:
            evs.append(e2); G.append(g2)
            if len(G) > MAX_MANIFOLD:
                raise RuntimeError(f"ground manifold larger than {MAX_MANIFOLD}")
        else:
            gap = e2 - evs[0]
            break
    return np.array(evs), np.column_stack(G), float(gap)


def _apply_J_columns(up_S, dn_S, G, L, t):
    """J (full, reaches outside S) applied to each column of G; returned on the common support W
    (union of the column supports) as a |W| x d matrix, with W's (up, dn, keys)."""
    outs = [apply_current_vec(up_S, dn_S, G[:, j].astype(complex), L, t) for j in range(G.shape[1])]
    keysW = np.unique(np.concatenate([o[3] for o in outs]))
    VW = np.zeros((len(keysW), G.shape[1]), dtype=complex)
    for j, (_, _, amp, uk) in enumerate(outs):
        VW[np.searchsorted(keysW, uk), j] = amp
    upW = (keysW >> L).astype(np.int64); dnW = (keysW & ((1 << L) - 1)).astype(np.int64)
    return upW, dnW, keysW, VW


def _herm_range(Q):
    """[min, max] of x^H Q x over unit x (Q Hermitian d x d)."""
    w = np.linalg.eigvalsh(0.5 * (Q + Q.conj().T))
    return [float(w[0]), float(w[-1])]


def _abs_range(r):
    lo, hi = r
    if lo <= 0.0 <= hi:
        return [0.0, max(-lo, hi)]
    return sorted([abs(lo), abs(hi)])


def _relerr_range(r, ref):
    lo, hi = r
    if lo <= ref <= hi:
        return [0.0, max(abs(lo - ref), abs(hi - ref)) / abs(ref)]
    return sorted([abs(lo - ref) / abs(ref), abs(hi - ref) / abs(ref)])


def screen_on_subspace(up_S, dn_S, L, t, U):
    """Given the SAMPLED distinct configs S=(up_S,dn_S), compute the shared-state test (classical
    post-processing of the sampled subspace): subspace ground level of H_S, e0, within-S moments
    m_bar_{0,1}, full-operator moments m_hat_{0,1}, Delta_{0,1}.

    Every moment is a Hermitian quadratic form in the subspace ground state.  If the ground level is
    degenerate (gap < DEGEN_TOL) the point values below refer to one arbitrary (seeded) vector of the
    manifold and are NOT well defined; the '<name>_manifold' entries give [min, max] over all
    normalized states of the ground manifold.  For a non-degenerate level min = max = point value."""
    n = len(up_S)
    keyS = _pack(up_S, dn_S, L)
    o = np.argsort(keyS); keyS_sorted = keyS[o]; orderS = np.arange(n)[o]

    H_S = build_op_on_set(up_S, dn_S, keyS_sorted, orderS, L, t, U, 'hop')
    J_S = build_op_on_set(up_S, dn_S, keyS_sorted, orderS, L, t, 0.0, 'cur')
    H_S = 0.5 * (H_S + H_S.getH()); J_S = 0.5 * (J_S + J_S.getH())     # symmetrize (numerical hygiene)

    evs, G, gap = subspace_ground_manifold(H_S)
    d = G.shape[1]
    e0 = float(evs[0])
    # connected components of the sampled-subspace hopping graph (diagnostic for degeneracy)
    from scipy.sparse.csgraph import connected_components
    ncomp = int(connected_components(abs(H_S) > 0, directed=False)[0])

    # RECONSTRUCTION (within-S), as d x d forms on the manifold
    VS = J_S @ G
    Qbar0 = VS.conj().T @ VS
    Qbar1 = VS.conj().T @ (H_S @ VS) - e0 * Qbar0
    # INDEPENDENT branch (full J,H reach OUTSIDE S) on the SAME states
    upW, dnW, keysW, VW = _apply_J_columns(up_S, dn_S, G, L, t)
    ow = np.argsort(keysW); orderW = np.arange(len(keysW))[ow]   # keysW already sorted
    H_W = build_op_on_set(upW, dnW, keysW[ow], orderW, L, t, U, 'hop')
    H_W = 0.5 * (H_W + H_W.getH())
    Qhat0 = VW.conj().T @ VW
    Qhat1 = VW.conj().T @ (H_W @ VW) - e0 * Qhat0

    # point values on the first manifold vector (the unique ground state when d = 1)
    m_bar_0 = float(Qbar0[0, 0].real); m_bar_1 = float(Qbar1[0, 0].real)
    m_hat_0 = float(Qhat0[0, 0].real); m_hat_1 = float(Qhat1[0, 0].real)
    leak_frac = max(0.0, 1.0 - m_bar_0 / m_hat_0) if m_hat_0 > 0 else 0.0

    # ranges over the ground manifold
    import scipy.linalg as sla
    ratio = sla.eigh(0.5 * (Qbar0 + Qbar0.conj().T), 0.5 * (Qhat0 + Qhat0.conj().T), eigvals_only=True)
    man = dict(m_bar_0=_herm_range(Qbar0), m_hat_0=_herm_range(Qhat0),
               m_bar_1=_herm_range(Qbar1), m_hat_1=_herm_range(Qhat1),
               delta_0=_abs_range(_herm_range(Qhat0 - Qbar0)),
               delta_1=_abs_range(_herm_range(Qhat1 - Qbar1)),
               leak_frac=[max(0.0, 1.0 - float(ratio[-1])), max(0.0, 1.0 - float(ratio[0]))])
    return dict(n=n, e0=e0, m_bar_0=m_bar_0, m_bar_1=m_bar_1, m_hat_0=m_hat_0,
                m_hat_1=m_hat_1, delta_0=abs(m_hat_0 - m_bar_0), delta_1=abs(m_hat_1 - m_bar_1),
                leak_frac=leak_frac, W=len(keysW), manifold_dim=d, gs_gap=gap,
                degenerate=bool(d > 1), manifold_evals=[float(x) for x in evs],
                n_components=ncomp, manifold=man)


# ================================ exact-diagonalization truth (calibration L) ================
def ed_true_moments(L, U, t=1.0):
    """FULL ED true current moments M_0,M_1 on the exact GS (only for ED-feasible calibration L)."""
    def strings(n):
        S = [sum(1 << b for b in c) for c in combinations(range(L), n)]; S.sort(); return S
    nup = nd = L // 2
    Su = np.array(strings(nup), dtype=np.int64); Sd = np.array(strings(nd), dtype=np.int64)
    Du, Dd = len(Su), len(Sd)
    ku = np.argsort(Su); ksu = Su[ku]; ku_ord = np.arange(Du)[ku]
    kd = np.argsort(Sd); ksd = Sd[kd]; kd_ord = np.arange(Dd)[kd]

    def single_hop(S, ks, kord, kind):
        n = len(S); ar = np.arange(n); r = []; c = []; v = []
        for i in range(L - 1):
            j = i + 1
            dirs = ((i, j, -t + 0j), (j, i, -t + 0j)) if kind == 'hop' else ((j, i, -1j * t), (i, j, +1j * t))
            for (p, q, coeff) in dirs:
                betw = between_mask(p, q); bp = 1 << p; bq = 1 << q
                m = (((S >> q) & 1) == 1) & (((S >> p) & 1) == 0); idx = ar[m]
                if idx.size:
                    ns = (S[idx] & ~bq) | bp
                    pos = np.searchsorted(ks, ns); pos = np.clip(pos, 0, n - 1); hit = ks[pos] == ns
                    r.append(kord[pos[hit]]); c.append(idx[hit]); v.append((coeff * _sig(S[idx], betw))[hit])
        if not r:
            return sp.csr_matrix((n, n), dtype=complex)
        return sp.csr_matrix((np.concatenate(v), (np.concatenate(r), np.concatenate(c))), shape=(n, n))

    Tu = single_hop(Su, ksu, ku_ord, 'hop'); Td = single_hop(Sd, ksd, kd_ord, 'hop')
    Ju = single_hop(Su, ksu, ku_ord, 'cur'); Jd = single_hop(Sd, ksd, kd_ord, 'cur')
    upocc = np.array([[(m >> i) & 1 for i in range(L)] for m in Su], float)
    dnocc = np.array([[(m >> i) & 1 for i in range(L)] for m in Sd], float)
    diagU = U * (upocc @ dnocc.T).ravel()
    H = (sp.kron(Tu, sp.identity(Dd)) + sp.kron(sp.identity(Du), Td) + sp.diags(diagU)).tocsr()
    J = (sp.kron(Ju, sp.identity(Dd)) + sp.kron(sp.identity(Du), Jd)).tocsr()
    Hh = 0.5 * (H + H.getH())
    ev, V = eigsh(Hh, k=1, which='SA', maxiter=5000); E0 = float(ev[0].real); psi0 = V[:, 0]
    Jp = J @ psi0
    M0 = float(np.real(np.vdot(Jp, Jp))); M1 = float(np.real(np.vdot(Jp, H @ Jp)) - E0 * M0)
    return E0, M0, M1, Du * Dd, psi0, (Su, Sd, Du, Dd)


# ================================ DMRG (TeNPy): GS, sampling, true moments ================
def dmrg_state_and_truth(L, U, t=1.0, chi=160, svd_min=1e-9, max_sweeps=16):
    from tenpy.models.hubbard import FermiHubbardModel
    from tenpy.models.model import CouplingModel
    from tenpy.networks.mps import MPS
    from tenpy.algorithms import dmrg
    import warnings; warnings.filterwarnings('ignore')
    M = FermiHubbardModel(dict(L=L, t=t, U=U, mu=0.0, bc_MPS='finite',
                               cons_N='N', cons_Sz='Sz', lattice='Chain'))
    prod = [('up' if i % 2 == 0 else 'down') for i in range(L)]
    psi = MPS.from_product_state(M.lat.mps_sites(), prod, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, M, dict(
        trunc_params=dict(chi_max=chi, svd_min=svd_min),
        max_sweeps=max_sweeps, min_sweeps=1, mixer=True, combine=True))
    # sweep-by-sweep with explicit early stopping (avoids running to max_sweeps needlessly)
    E0 = None; nsweeps = 0
    for i in range(max_sweeps):
        eng.options['max_sweeps'] = i + 1; eng.options['min_sweeps'] = i + 1
        E0, psi = eng.run(); nsweeps = i + 1
        Es = eng.sweep_stats['E']
        if i >= 3 and abs(Es[-1] - Es[-2]) < 1e-9:
            break
    psi.canonical_form()
    try:
        trunc = float(eng.sweep_stats['max_trunc_err'][-1])
    except Exception:
        trunc = float('nan')
    maxchi = int(max(psi.chi))
    # current MPO via a CouplingModel on the same lattice
    cm = CouplingModel(M.lat)
    for op_d, op in (('Cdu', 'Cu'), ('Cdd', 'Cd')):
        cm.add_coupling(1j * t, 0, op_d, 0, op, 1, plus_hc=True)
    Jmpo = cm.calc_H_MPO()
    J0 = psi.copy()
    Jmpo.apply(J0, dict(compression_method='SVD',
                        trunc_params=dict(chi_max=max(4 * chi, 800), svd_min=1e-12)))
    M0 = float(J0.norm ** 2)
    M1 = float(M0 * (M.H_MPO.expectation_value(J0) - E0))
    return dict(model=M, psi=psi, E0=float(E0), maxchi=maxchi, trunc=trunc,
                nsweeps=nsweeps, M0=M0, M1=M1)


def sample_configs(psi, model, L, n_shots, seed=0, partial_path=None, flush_every=2000):
    """Perfect-sample n_shots occupation configs ~ |<config|0>|^2; return the ordered list of packed
    keys (a prefix is a valid smaller-budget sample).  Resumable: if partial_path exists, continue
    from there; flush progress every `flush_every` shots so a killed job loses nothing."""
    site0 = model.lat.mps_sites()[0]
    lbl = {v: k for k, v in site0.state_labels.items()}
    keys = np.full(n_shots, -1, dtype=np.int64)
    start = 0
    if partial_path and os.path.exists(partial_path):
        prev = np.load(partial_path)
        m = int((prev >= 0).sum()); m = min(m, n_shots)
        keys[:m] = prev[:m]; start = m
        log(f"[sample] resumed from {start} shots")
    # i.i.d. draws: a fresh rng stream offset by `start` gives valid new samples (no need to replay)
    rng = np.random.default_rng([seed, start])
    for s in range(start, n_shots):
        sig, _ = psi.sample_measurements(rng=rng)
        up = 0; dn = 0
        for i, si in enumerate(sig):
            lab = lbl[int(si)]
            if lab in ('up', 'full'):
                up |= (1 << i)
            if lab in ('down', 'full'):
                dn |= (1 << i)
        keys[s] = (up << L) | dn
        if partial_path and (s + 1) % flush_every == 0:
            np.save(partial_path, keys)
            log(f"[sample]   {s+1}/{n_shots} shots  ({len(np.unique(keys[keys>=0]))} distinct)")
    if partial_path:
        np.save(partial_path, keys)
    return keys


def _pickle_paths(cpath):
    return cpath.replace('.npz', '_psi.pkl'), cpath.replace('.npz', '_keys.partial.npy')


def unpack(keys, L):
    up = (keys >> L).astype(np.int64); dn = (keys & ((1 << L) - 1)).astype(np.int64)
    return up, dn


# ==================================== driver ====================================
def run_L(L, U=8.0, t=1.0, budgets=None, n_max=None, chi=400, seed=0, ed_calib=False):
    from math import comb
    nup = L // 2
    sector = comb(L, nup) ** 2
    log(f"===== L={L}  U/t={U}  half-filling (nup=ndn={nup})  sector dim = C({L},{nup})^2 = {sector:.3e} =====")

    truth = {}
    if ed_calib:
        E0e, M0e, M1e, dim, psi0e, meta = ed_true_moments(L, U, t)
        log(f"[ED] exact  E0={E0e:.6f}  M0={M0e:.6f}  M1={M1e:.6f}  (full dim {dim})")
        truth = dict(E0=E0e, M0=M0e, M1=M1e, source='ED-exact')

    if budgets is None:
        budgets = [500, 1000, 2000, 5000, 10000, 20000, 50000]
    if n_max is None:
        n_max = max(budgets)

    # ---- cache the expensive DMRG + sampling stage (keys sequence + true moments) ----
    os.makedirs(CACHE, exist_ok=True)
    cpath = os.path.join(CACHE, f'beyondED_L{L}_U{U}_chi{chi}_n{n_max}_s{seed}.npz')
    ppkl, ppart = _pickle_paths(cpath)
    if os.path.exists(cpath):
        z = np.load(cpath)
        dmrg = dict(E0=float(z['E0']), maxchi=int(z['maxchi']), trunc=float(z['trunc']),
                    M0=float(z['M0']), M1=float(z['M1']))
        allkeys = z['allkeys']
        log(f"[cache] loaded {cpath}")
    else:
        import pickle
        if os.path.exists(ppkl):
            with open(ppkl, 'rb') as f:
                d = pickle.load(f)
            log(f"[cache] loaded DMRG state {ppkl}")
        else:
            d = dmrg_state_and_truth(L, U, t, chi=chi)
            with open(ppkl, 'wb') as f:
                pickle.dump(dict(psi=d['psi'], model=d['model'], E0=d['E0'], maxchi=d['maxchi'],
                                 trunc=d['trunc'], M0=d['M0'], M1=d['M1']), f)
            log(f"[cache] wrote DMRG state {ppkl}")
        log(f"[DMRG] E0={d['E0']:.6f}  maxchi={d['maxchi']}  trunc_eps={d['trunc']:.2e}  "
            f"M0(true)={d['M0']:.6f}  M1(true)={d['M1']:.6f}")
        log(f"[sample] drawing {n_max} perfect-sampling shots from the DMRG GS ...")
        ts = time.time()
        allkeys = sample_configs(d['psi'], d['model'], L, n_max, seed=seed,
                                 partial_path=ppart, flush_every=2000)
        log(f"[sample] done in {time.time()-ts:.1f}s; "
            f"{len(np.unique(allkeys))} distinct configs at full budget")
        dmrg = dict(E0=d['E0'], maxchi=d['maxchi'], trunc=d['trunc'], M0=d['M0'], M1=d['M1'])
        np.savez(cpath, allkeys=allkeys, E0=dmrg['E0'], maxchi=dmrg['maxchi'],
                 trunc=dmrg['trunc'], M0=dmrg['M0'], M1=dmrg['M1'])
        log(f"[cache] wrote {cpath}")
    d = dmrg
    log(f"[DMRG] E0={d['E0']:.6f}  maxchi={d['maxchi']}  trunc_eps={d['trunc']:.2e}  "
        f"M0(true)={d['M0']:.6f}  M1(true)={d['M1']:.6f}")
    if ed_calib:
        log(f"[calib] |E0_DMRG-E0_ED|={abs(d['E0']-truth['E0']):.2e}  "
            f"|M0-M0|={abs(d['M0']-truth['M0']):.2e}  |M1-M1|={abs(d['M1']-truth['M1']):.2e}")
    M0_true = truth['M0'] if ed_calib else d['M0']
    M1_true = truth['M1'] if ed_calib else d['M1']

    rows = []
    for N in budgets:
        keys = allkeys[:N]
        uk = np.unique(keys)
        up_S, dn_S = unpack(uk, L)
        sc = screen_on_subspace(up_S, dn_S, L, t, U)
        rec_err0 = abs(sc['m_bar_0'] - M0_true) / abs(M0_true)
        rec_err1 = abs(sc['m_bar_1'] - M1_true) / abs(M1_true)
        hat_err0 = abs(sc['m_hat_0'] - M0_true) / abs(M0_true)
        hat_err1 = abs(sc['m_hat_1'] - M1_true) / abs(M1_true)
        man = dict(sc['manifold'])
        man['recon_err_m0'] = _relerr_range(man['m_bar_0'], M0_true)
        man['recon_err_m1'] = _relerr_range(man['m_bar_1'], M1_true)
        man['hat_err_m0'] = _relerr_range(man['m_hat_0'], M0_true)
        man['hat_err_m1'] = _relerr_range(man['m_hat_1'], M1_true)
        row = dict(N_shots=N, S=int(sc['n']), W=int(sc['W']), e0=sc['e0'],
                   degenerate=sc['degenerate'], manifold_dim=sc['manifold_dim'],
                   gs_gap=(sc['gs_gap'] if np.isfinite(sc['gs_gap']) else None),
                   manifold_evals=sc['manifold_evals'], n_components=sc['n_components'],
                   m_bar_0=sc['m_bar_0'], m_hat_0=sc['m_hat_0'], delta_0=sc['delta_0'],
                   m_bar_1=sc['m_bar_1'], m_hat_1=sc['m_hat_1'], delta_1=sc['delta_1'],
                   leak_frac=sc['leak_frac'],
                   recon_err_m0=rec_err0, recon_err_m1=rec_err1,
                   hat_err_m0=hat_err0, hat_err_m1=hat_err1,
                   manifold=man)
        rows.append(row)
        if sc['degenerate']:
            log(f"  N={N:6d} DEGENERATE ground level of H_S: dim={sc['manifold_dim']} "
                f"gap={sc['gs_gap']:.3e} components={sc['n_components']}  "
                f"Delta0 in [{man['delta_0'][0]:.4f},{man['delta_0'][1]:.4f}]  "
                f"m_hat0 in [{man['m_hat_0'][0]:.4f},{man['m_hat_0'][1]:.4f}]  "
                f"m_bar0 in [{man['m_bar_0'][0]:.4f},{man['m_bar_0'][1]:.4f}]")
        log(f"  N={N:6d} |S|={sc['n']:6d} |W|={sc['W']:7d} e0={sc['e0']:8.4f} | "
            f"m_bar0={sc['m_bar_0']:7.3f} m_hat0={sc['m_hat_0']:7.3f} D0={sc['delta_0']:7.3f} "
            f"leak={100*sc['leak_frac']:5.1f}% | m_bar1={sc['m_bar_1']:7.2f} m_hat1={sc['m_hat_1']:7.2f} "
            f"D1={sc['delta_1']:7.2f} | reconErr(m0)={100*rec_err0:5.1f}% (m1)={100*rec_err1:5.1f}%")
    return dict(L=L, U=U, sector_dim=float(sector), ed_calib=ed_calib, chi_max=chi, seed=seed,
                E0_DMRG=d['E0'], maxchi=d['maxchi'],
                trunc=(d['trunc'] if np.isfinite(d['trunc']) else None),
                M0_true=M0_true, M1_true=M1_true, truth_source=('ED' if ed_calib else 'DMRG'),
                E0_ED=(truth['E0'] if ed_calib else None),
                dmrg_M0=d['M0'], dmrg_M1=d['M1'], sweep=rows)


def _fmt_verdict(r, sm):
    """One factual block per run, formatted from the data."""
    s = r['sweep']; nd = [w for w in s if not w['degenerate']]; dg = [w for w in s if w['degenerate']]
    a, b = nd[0], nd[-1]
    txt = (f"L={r['L']} (U/t={r['U']:g}, S_z=0 sector dim {r['sector_dim']:.3e}; "
           f"truth from {r['truth_source']}: m0={r['M0_true']:.2f}, m1={r['M1_true']:.2f}; "
           f"DMRG E0={r['E0_DMRG']:.4f}, max bond dim {r['maxchi']}"
           + (f", truncation error {r['trunc']:.1e}" if r['trunc'] is not None else "") + "). "
           f"Over the {len(nd)} non-degenerate budgets N={a['N_shots']}..{b['N_shots']} "
           f"(|S|={a['S']}..{b['S']}): Delta_0 {a['delta_0']:.2f} -> {b['delta_0']:.2f} "
           f"(monotone decreasing: {sm['delta0_monotonic_decreasing']}); "
           f"relative error of m_bar_0 {a['recon_err_m0']:.2f} -> {b['recon_err_m0']:.2f}; "
           f"relative error of m_hat_0 {a['hat_err_m0']:.3f} -> {b['hat_err_m0']:.3f}; "
           f"leak fraction {a['leak_frac']:.2f} -> {b['leak_frac']:.2f}; "
           f"subspace e0 {a['e0']:.3f} -> {b['e0']:.3f}; "
           + (f"Pearson corr(Delta_0, rel. error of m_bar_0) = {sm['corr_delta0_reconErr']:.4f} "
              f"over {sm['n_rows_nondegenerate']} rows." if sm['corr_delta0_reconErr'] is not None
              else "too few rows for a correlation."))
    for w in dg:
        m = w['manifold']
        txt += (f" At N={w['N_shots']} (|S|={w['S']}) the subspace ground level is {w['manifold_dim']}-fold "
                f"degenerate (next gap {w['gs_gap']:.3g}); over it Delta_0 ranges "
                f"{m['delta_0'][0]:.2f}-{m['delta_0'][1]:.2f}, m_hat_0 {m['m_hat_0'][0]:.2f}-{m['m_hat_0'][1]:.2f}, "
                f"relative error of m_bar_0 {m['recon_err_m0'][0]:.3f}-{m['recon_err_m0'][1]:.3f}.")
    return txt


def main():
    import tenpy, scipy
    out = {'_provenance': {'script': 'src/beyond_ed_dmrg_catch.py', 'date': DATE,
                           'model': '1D Fermi-Hubbard OBC, half filling, U/t=8, probe=current J',
                           'test': 'shared-state moment test (classical post-processing of the sampled '
                                   'subspace): Delta_k=|m_hat_k-m_bar_k| on the subspace ground state |0_S>',
                           'validator': 'DMRG/MPS moments M_k=<0|J(H-E0)^k J|0> (ED at L=8)',
                           'degeneracy_rule': (f'ground level of H_S degenerate if next eigenvalue - E_min < {DEGEN_TOL:g}; '
                                               'point values of a degenerate row refer to one arbitrary manifold vector; '
                                               'the manifold entry gives [min,max] over the ground manifold; '
                                               'monotonicity and correlations use non-degenerate rows only'),
                           'versions': {'python': sys.version.split()[0], 'numpy': np.__version__,
                                        'scipy': scipy.__version__, 'tenpy': tenpy.__version__}},
           'runs': []}
    # (1) L=8 ED calibration: locks conventions in-run.
    out['runs'].append(run_L(8, U=8.0, budgets=[200, 500, 1000, 2000, 4000], n_max=4000,
                             chi=200, ed_calib=True))
    # (2) L=20 beyond ED (sector 3.4e10).
    out['runs'].append(run_L(20, U=8.0, budgets=[500, 1000, 2000, 5000, 10000, 15000],
                             n_max=15000, chi=160))
    # (3) L=24 beyond ED (sector 7.3e12).
    out['runs'].append(run_L(24, U=8.0, budgets=[500, 1000, 2000, 5000, 10000, 15000],
                             n_max=15000, chi=160))

    def _mono_dec(x):
        return all(x[i + 1] <= x[i] + 1e-9 for i in range(len(x) - 1))
    summ = []
    for r in out['runs']:
        s = [w for w in r['sweep'] if not w['degenerate']]
        d0 = [w['delta_0'] for w in s]; leak = [w['leak_frac'] for w in s]
        re0 = [w['recon_err_m0'] for w in s]
        cc = float(np.corrcoef(d0, re0)[0, 1]) if len(d0) > 2 else None
        summ.append(dict(L=r['L'], sector_dim=r['sector_dim'], truth_source=r['truth_source'],
                         n_rows_nondegenerate=len(s),
                         degenerate_budgets=[w['N_shots'] for w in r['sweep'] if w['degenerate']],
                         delta0_monotonic_decreasing=_mono_dec(d0),
                         leak_monotonic_decreasing=_mono_dec(leak),
                         reconErr_monotonic_decreasing=_mono_dec(re0),
                         corr_delta0_reconErr=cc,
                         delta0_range=[d0[-1], d0[0]], leak_range=[leak[-1], leak[0]],
                         reconErr_m0_range=[re0[-1], re0[0]],
                         hatErr_m0_range=[min(w['hat_err_m0'] for w in s), max(w['hat_err_m0'] for w in s)]))
    out['summary'] = summ
    out['verdict'] = ' '.join(_fmt_verdict(r, sm) for r, sm in zip(out['runs'], summ)) + (
        " Scope: the open 1D Hubbard chain is DMRG-solvable, so this checks the shared-state test code at "
        "sector dimensions beyond full ED against DMRG references; it is not a measure of detection power "
        "and not a beyond-classical claim.")

    os.makedirs(RES, exist_ok=True)
    outpath = os.path.join(RES, f'{DATE}_beyond_ed_dmrg_catch.json')
    with open(outpath, 'w') as f:
        json.dump(out, f, indent=2, allow_nan=False)
    log(f"WROTE {outpath}")
    log("[verdict] " + out['verdict'])


if __name__ == '__main__':
    main()
