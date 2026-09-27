# -*- coding: utf-8 -*-
r"""within_sector_lp.py -- R3 + M1 (+ critic item 5): what within_sector_control.py only PRINTS, stored.

System: the L=6, U/t=4 doped current Lehmann measure of within_sector_control.py (same construction,
spectral_lanczos sector basis, degenerate eigenstates aggregated into N distinct poles).

R3  (a) LP min/max of Delta_{K+1} over redistributions delta preserving m_0..m_K exactly with
        w_n + delta_n >= 0, for K = 1, 2, 3 (Eq. wclp of the paper), with the optimal vectors;
    (b) null-space dimension of the order-0..K moment rows (N-K-1) -- the reason guaranteed power is zero;
    (c) sigma_min and kappa of the N x N node Vandermonde on RAW nodes and on nodes rescaled affinely to
        [-1,1] (critic item 5: kappa is unit-dependent; guaranteed power is zero because of the null space);
    (d) three explicit counterexample vectors (maximal ||delta||_1, exact vertex enumeration of the
        polytope, cross-checked by 2^N sign-pattern LPs):
          CE1  m0,m1,m2 preserved exactly;
          CE2  m0,m1 preserved exactly and |Delta_2| <= tau_2;
          CE3  |Delta_k| <= tau_k for k = 0,1,2 (passes the deployed battery);
        with both weight conventions: ||delta||_1/m0 and the transferred weight ||delta||_1/2 (as a
        fraction of m0), the emptied poles and the dominant residue before/after;
    (e) tau_k from the deployed rule tau_k = z sqrt(Var(O_k^loc)/N_s + (b m_k)^2), z=1.96, N_s=5e4, b=2%
        (the rule of interval_battery.py), recomputed here, plus tau_3, tau_4 by the same rule (context only:
        m3, m4 are not deployed; the local-estimator variance is a LOWER bound on the per-shot variance).
M1  shifted-power bound: if delta preserves m_0..m_{s-2} and is supported in a cluster of diameter D and
    midpoint c, then Delta_{s-1} = sum_n delta_n (omega_n - c)^{s-1}, so |Delta_{s-1}| <= rho (D/2)^{s-1},
    rho = ||delta||_1. Checked on CE1, CE2, the LP optimisers, the committed random cases of
    within_sector_control.py (seed 0), 20000 random null-space redistributions on random pole subsets
    (seed 20260927), and synthetic clusters (tightness and D^{s-1} scaling). Also Theorem 2(a):
    max_j |Delta_j| >= sigma_min(V) ||delta||_2 / sqrt(N-1-K).

Writes keys 'R3_within_sector_lp' and 'M1_shifted_power_bound' of data/2026-09-27_theory_numerics.json.
Does NOT touch data/2026-08-24_within_sector_control.json. Run: cd src && python within_sector_lp.py
"""
import os, sys, time, itertools
import numpy as np
import scipy.sparse as sp
from scipy.linalg import null_space
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import spectral_lanczos as sl
from small_checks import record_theory_numerics

L, U = 6, 4.0
NS, Z, BIAS = 50000, 1.96, 0.02
SEED_MC = 20260927


# ------------------------------------------------------------------------------------------------
def l6u4_system():
    """The L=6, U/t=4, (N_up,N_dn)=(2,2) ring sector of within_sector_control.py / interval_battery.py."""
    nup = nd = L // 3
    Tu, Su, iu = sl.hop(L, nup); Td, Sd, idd = sl.hop(L, nd)
    Du, Dd = len(Su), len(Sd)
    upocc = sl.occ_matrix(Su, L); dnocc = sl.occ_matrix(Sd, L)
    Hs = (sp.kron(Tu, sp.identity(Dd)) + sp.kron(sp.identity(Du), Td)
          + sp.diags(U * (upocc @ dnocc.T).ravel())).tocsr()
    Ju = sl.current_string(L, nup); Jd = sl.current_string(L, nd)
    Js = (sp.kron(Ju, sp.identity(Dd)) + sp.kron(sp.identity(Du), Jd)).tocsr()
    Hd = np.real(Hs.toarray())
    E, V = np.linalg.eigh(Hd)
    return dict(H=Hd, Hs=Hs, J=Js.toarray(), Js=Js, E=E, V=V, E0=float(E[0]), psi0=V[:, 0],
                gap=float(E[1] - E[0]), dim=Hd.shape[0])


def lehmann_measure(S):
    """distinct poles and aggregated weights, exactly as within_sector_control.py:42-59."""
    omega = S['E'] - S['E0']
    w = np.abs(S['V'].T @ (S['Js'] @ S['psi0'])) ** 2
    keep = (w > 1e-10) & (omega > 1e-8)
    omega, w = omega[keep], w[keep]
    n_eig = int(len(w))
    order = np.argsort(omega, kind='stable'); omega, w = omega[order], w[order]
    uom, uw = [], []
    for o, wi in zip(omega, w):
        if uom and abs(o - uom[-1]) < 1e-6:
            uw[-1] += wi
        else:
            uom.append(float(o)); uw.append(float(wi))
    return np.array(uom), np.array(uw), n_eig


def deployed_taus(S, kmax=4):
    """tau_k = z sqrt(Var(O_k^loc)/N_s + (b m_k)^2), O_k = J (H-E0)^k J, local estimator over |psi0(x)|^2
    (interval_battery.py:30-43). Returns m_k, var_k, tau_k (2% bias) and tau_k (no bias), k = 0..kmax."""
    psi0 = S['psi0']; Js = S['Js']
    Hc = S['Hs'] - S['E0'] * sp.identity(S['dim'])
    p = np.abs(psi0) ** 2; mask = np.abs(psi0) > 1e-12
    Jpsi = Js @ psi0
    v = Jpsi.copy()
    out = []
    for k in range(kmax + 1):
        Okpsi = Js @ v                                   # O_k psi0 = J (H-E0)^k J psi0
        loc = np.zeros(len(psi0), dtype=complex); loc[mask] = Okpsi[mask] / psi0[mask]
        mean = float(np.sum(p * loc).real); var = float(np.sum(p * np.abs(loc) ** 2) - mean ** 2)
        m_exact = float(np.vdot(Jpsi, v).real)           # <J psi0 | (H-E0)^k | J psi0>
        out.append(dict(k=k, m_exact=m_exact, m_local_mean=mean, var_local=var,
                        tau_bias2pct=Z * np.sqrt(max(var, 0) / NS + (BIAS * abs(mean)) ** 2),
                        tau_nobias=Z * np.sqrt(max(var, 0) / NS)))
        v = Hc @ v
    return out


# ------------------------------------------------------------------------------------------------
def lp_range(om, w, K):
    c = om ** (K + 1)
    Aeq = np.array([om ** j for j in range(K + 1)]); beq = np.zeros(K + 1)
    bnds = [(-wi, None) for wi in w]
    hi = linprog(-c, A_eq=Aeq, b_eq=beq, bounds=bnds, method='highs')
    lo = linprog(c, A_eq=Aeq, b_eq=beq, bounds=bnds, method='highs')
    return float(lo.fun), float(-hi.fun), lo.x, hi.x


def max_l1_vertex(om, w, tol):
    """max ||delta||_1 over the polytope {|sum delta omega^k| <= tol_k (k<len(tol)), delta >= -w}:
    a convex function is maximised at a vertex; enumerate all vertices (ws_vertex.py, exact data)."""
    N = len(w); rows = []
    for k, tk in enumerate(tol):
        a = om ** k
        rows.append((a, float(tk))); rows.append((-a, float(tk)))
    for n in range(N):
        e = np.zeros(N); e[n] = -1.0
        rows.append((e, float(w[n])))
    A = np.array([r[0] for r in rows]); b = np.array([r[1] for r in rows])
    scale = np.maximum(1.0, np.abs(A).max(axis=1))
    best = (0.0, None); nvert = 0
    for Sidx in itertools.combinations(range(len(rows)), N):
        M = A[list(Sidx)]
        if abs(np.linalg.det(M / scale[list(Sidx), None])) < 1e-12:
            continue
        x = np.linalg.solve(M, b[list(Sidx)])
        if np.all(A @ x <= b + 1e-9 * scale):
            nvert += 1
            l1 = float(np.abs(x).sum())
            if l1 > best[0] + 1e-12:
                best = (l1, x)
    return best[0], best[1], nvert


def max_l1_signpattern(om, w, tol):
    """cross-check: max over sign patterns s of the LP max s.delta (lp_check.py 'corrected')."""
    N = len(w); best = 0.0
    Aub = []; bub = []
    for k, tk in enumerate(tol):
        Aub.append(om ** k); bub.append(tk); Aub.append(-om ** k); bub.append(tk)
    Aub = np.array(Aub); bub = np.array(bub)
    for s in itertools.product([1, -1], repeat=N):
        s = np.array(s)
        bnds = [(0, None) if s[i] > 0 else (-w[i], 0) for i in range(N)]
        r = linprog(-s, A_ub=Aub, b_ub=bub, bounds=bnds, method='highs')
        if r.status == 0 and -r.fun > best:
            best = float(-r.fun)
    return best


def describe(om, w, d, taus):
    m0 = float(w.sum()); wc = w + d
    D = [float(np.sum(d * om ** k)) for k in range(6)]
    tau = [t['tau_bias2pct'] for t in taus]
    l1 = float(np.abs(d).sum())
    idom = int(np.argmax(w))
    passes = all(abs(D[k]) <= tau[k] * (1 + 1e-9) + 1e-12 for k in range(3))
    over = [k for k in range(len(tau)) if abs(D[k]) > tau[k] * (1 + 1e-9) + 1e-12]
    return dict(delta=d, w_true=w, w_corrupt=wc, Delta_signed_k0_to_5=D,
                abs_Delta_over_tau_k0_to_4=[abs(D[k]) / tau[k] for k in range(len(tau))],
                l1=l1, l1_over_m0=l1 / m0, transferred_half_l1=l1 / 2, transferred_half_l1_over_m0=l1 / 2 / m0,
                n_poles_emptied=int(np.sum(wc < 1e-9 * m0)), min_corrupt_weight=float(wc.min()),
                dominant_pole_omega=float(om[idom]), dominant_residue_true=float(w[idom]),
                dominant_residue_corrupt=float(wc[idom]),
                passes_deployed_battery_m0_m1_m2=bool(passes),
                orders_exceeding_same_rule_tau_k_le_4=over,
                corrupted_measure_positive=bool(np.all(wc >= -1e-12)))


# ------------------------------------------------------------------------------------------------
def shifted_power_check(om, d, tag):
    """M1: s-2 = highest order such that m_0..m_{s-2} are preserved exactly (relative 1e-9);
    D, c from the support of delta; verify the identity and the bound for Delta_{s-1}."""
    scale = [np.sum(np.abs(d) * np.abs(om) ** k) for k in range(len(om) + 1)]
    D = [float(np.sum(d * om ** k)) for k in range(len(om) + 1)]
    kpres = -1
    for k in range(len(om)):
        if abs(D[k]) <= 1e-9 * max(scale[k], 1e-300):
            kpres = k
        else:
            break
    s = kpres + 2
    supp = np.abs(d) > 1e-12 * np.abs(d).max()
    lo, hi = om[supp].min(), om[supp].max(); Dm = float(hi - lo); c = 0.5 * (lo + hi)
    rho = float(np.abs(d).sum())
    lhs = D[s - 1]
    shifted = float(np.sum(d * (om - c) ** (s - 1)))
    bound = rho * (Dm / 2) ** (s - 1)
    return dict(case=tag, preserved_through_order=kpres, s=s, order_tested=s - 1,
                Delta_s_minus_1=lhs, shifted_form=shifted, identity_abs_err=abs(lhs - shifted),
                identity_rel_err=abs(lhs - shifted) / max(abs(lhs), 1e-300),
                rho_l1=rho, cluster_diameter_D=Dm, midpoint_c=c, n_support=int(supp.sum()),
                bound_rho_Dhalf_pow=bound, ratio_abs_Delta_over_bound=abs(lhs) / bound,
                bound_holds=bool(abs(lhs) <= bound * (1 + 1e-9)))


def theorem_a_check(om, d, K, sigma_min, tag):
    N = len(om)
    V = np.vander(om, N, increasing=True).T
    Dv = V @ d
    lhs = float(np.max(np.abs(Dv[K + 1:])))
    rhs = float(sigma_min * np.linalg.norm(d) / np.sqrt(N - 1 - K))
    return dict(case=tag, K=K, max_j_abs_Delta_j=lhs, sigma_min_l2_over_sqrt=rhs, holds=bool(lhs >= rhs * (1 - 1e-9)),
                l2=float(np.linalg.norm(d)), l1=float(np.abs(d).sum()),
                l1_over_l2=float(np.abs(d).sum() / np.linalg.norm(d)), sqrtN=float(np.sqrt(N)))


def redistribute_preserving(om, w, K, seed=0):
    """verbatim copy of within_sector_control.py:redistribute_preserving (committed random cases)."""
    Vand = np.array([om ** j for j in range(K + 1)])
    ns = null_space(Vand)
    rng = np.random.default_rng(seed)
    dw = ns @ rng.standard_normal(ns.shape[1])
    dw = dw / np.max(np.abs(dw))
    neg = dw < 0
    alpha = 0.8 * (np.min(-w[neg] / dw[neg]) if neg.any() else 1.0)
    return alpha * dw


# ------------------------------------------------------------------------------------------------
def main():
    t0 = time.time()
    S = l6u4_system()
    om, w, n_eig = lehmann_measure(S)
    N = len(w); m0 = float(w.sum())
    mom = [float(np.sum(w * om ** k)) for k in range(6)]
    taus = deployed_taus(S, kmax=4)
    tau = [t['tau_bias2pct'] for t in taus]
    print(f"L={L} U/t={U}: {n_eig} weighted eigenstates -> N={N} distinct poles {np.round(om, 4)}")
    print(f"weights {np.round(w, 5)}  m0..m2 = {[round(x, 5) for x in mom[:3]]}")
    print(f"tau (2% bias, Ns=5e4) k=0..4: {[round(x, 4) for x in tau]}")

    # ---- (a) LP ranges ----
    lp = []
    for K in (1, 2, 3):
        lo, hi, xlo, xhi = lp_range(om, w, K)
        lp.append(dict(K=K, order=K + 1, Delta_min=lo, Delta_max=hi, contains_zero=bool(lo < 0 < hi),
                       tau_same_rule=tau[K + 1], max_over_tau=hi / tau[K + 1],
                       argmin_delta=xlo, argmax_delta=xhi))
        print(f"K={K}: Delta_{K+1} in [{lo:.4f}, {hi:.4f}]  (tau_{K+1} = {tau[K+1]:.4f})")

    # ---- (b) null spaces ----
    nulls = {f'rows_0..{K}': int(null_space(np.array([om ** j for j in range(K + 1)])).shape[1]) for K in range(N)}

    # ---- (c) conditioning ----
    V = np.vander(om, N, increasing=True).T
    sv = np.linalg.svd(V, compute_uv=False)
    a, b = om.min(), om.max(); x = (2 * om - (a + b)) / (b - a)
    Vx = np.vander(x, N, increasing=True).T; svx = np.linalg.svd(Vx, compute_uv=False)
    x2 = om / om.max(); V2 = np.vander(x2, N, increasing=True).T; sv2 = np.linalg.svd(V2, compute_uv=False)
    V3 = np.array([om ** k for k in range(3)]); sv3 = np.linalg.svd(V3, compute_uv=False)
    cond = dict(raw=dict(sigma_min=float(sv[-1]), sigma_max=float(sv[0]), kappa=float(sv[0] / sv[-1])),
                rescaled_minus1_1=dict(sigma_min=float(svx[-1]), sigma_max=float(svx[0]), kappa=float(svx[0] / svx[-1])),
                scaled_by_max=dict(sigma_min=float(sv2[-1]), sigma_max=float(sv2[0]), kappa=float(sv2[0] / sv2[-1])),
                deployed_rows_0_2_singular_values=sv3, min_pole_separation=float(np.min(np.diff(om))),
                committed_json_values={'sigma_min_V': 0.011805717328315082, 'cond_V': 1248380.672393826})
    print(f"raw sigma_min={sv[-1]:.4e} kappa={sv[0]/sv[-1]:.4e}; rescaled[-1,1] sigma_min={svx[-1]:.4f} "
          f"kappa={svx[0]/svx[-1]:.2f}")

    # ---- (d) counterexamples ----
    ces = {}
    specs = [('CE1_m0_m1_m2_exact', [0.0, 0.0, 0.0]),
             ('CE2_m0_m1_exact_absD2_le_tau2', [0.0, 0.0, tau[2]]),
             ('CE3_all_within_tau', [tau[0], tau[1], tau[2]])]
    for tag, tol in specs:
        l1, d, nvert = max_l1_vertex(om, w, tol)
        l1_sp = max_l1_signpattern(om, w, tol)
        rec = describe(om, w, d, taus)
        rec.update(tolerances=tol, n_vertices=nvert, l1_signpattern_crosscheck=l1_sp,
                   crosscheck_abs_diff=abs(l1 - l1_sp))
        ces[tag] = rec
        print(f"{tag}: ||d||_1/m0 = {100*rec['l1_over_m0']:.1f}%  transferred = "
              f"{100*rec['transferred_half_l1_over_m0']:.1f}%  emptied={rec['n_poles_emptied']}  "
              f"dom {rec['dominant_residue_true']:.3f}->{rec['dominant_residue_corrupt']:.3f}  "
              f"passes={rec['passes_deployed_battery_m0_m1_m2']}  (sign-LP {100*l1_sp/m0:.1f}%)")
    # blind-spot size vs exactly preserved order
    exact_pres = []
    for j in range(N):
        l1, d, _ = max_l1_vertex(om, w, [0.0] * (j + 1))
        exact_pres.append(dict(preserved_exactly='m0..m%d' % j, max_l1_over_m0=l1 / m0,
                               max_transferred_over_m0=l1 / 2 / m0))
    # committed K=1 random case of within_sector_control.py (C5 relabel)
    committed = []
    for K in (1, 2, 3):
        d = redistribute_preserving(om, w, K, seed=0)
        committed.append(dict(K=K, l1_over_m0=float(np.abs(d).sum() / m0),
                              transferred_half_l1_over_m0=float(np.abs(d).sum() / 2 / m0),
                              Delta_next=float(np.sum(d * om ** (K + 1))),
                              Delta_next_over_tau=float(abs(np.sum(d * om ** (K + 1))) / tau[K + 1]), delta=d))

    guaranteed = dict(
        statement='min over m0..mK-preserving delta with ||delta||_1 >= rho of |Delta_{K+1}| is 0 whenever the '
                  'order-0..K+1 rows have a nontrivial null space (N-K-2 >= 1), independent of conditioning',
        deployed_K1_zero_power_up_to_l1_over_m0=ces['CE1_m0_m1_m2_exact']['l1_over_m0'],
        reason='CE1 preserves m0,m1,m2 exactly; alpha*CE1 (0<alpha<=1) is feasible, so |Delta_2| = 0 for every '
               'rho up to ||CE1||_1',
        null_space_dim_deployed_rows_0_2=nulls['rows_0..2'])

    R3 = dict(
        system='L=6, U/t=4 doped current (2,2) ring, spectral_lanczos basis (within_sector_control.py)',
        n_weighted_eigenstates=n_eig, N_distinct_poles=N, poles=om, weights=w, moments_m0_to_m5=mom,
        ground_state_sector_gap=S['gap'],
        tau_rule='tau_k = 1.96 sqrt(Var(O_k^loc)/5e4 + (0.02 m_k)^2), O_k = J(H-E0)^k J (interval_battery.py)',
        tau=taus, tau_deployed_k0_to_2=tau[:3],
        lp_ranges=lp, null_space_dims=nulls, conditioning=cond,
        counterexamples=ces, blind_spot_vs_exactly_preserved_order=exact_pres,
        committed_random_cases_seed0=committed, guaranteed_power=guaranteed,
        weight_conventions={'l1_over_m0': '||delta||_1 / m0',
                            'transferred_half_l1_over_m0': '(||delta||_1 / 2) / m0 = weight moved from '
                                                           'depleted to enhanced poles (sum delta = 0)'},
        plan_expectation='LP [-0.9207,7.0913],[-0.1364,5.2638],[-0.3033,0.3917]; CE 143.3/171.8/203.4% of m0 '
                         '(71.7/85.9/101.7% moved); raw sigma_min 1.18e-2, kappa 1.25e6; rescaled kappa 50, '
                         'sigma_min 0.056; null dim 2')

    # ---------------- M1 ----------------
    checks = [shifted_power_check(om, ces['CE1_m0_m1_m2_exact']['delta'], 'CE1'),
              shifted_power_check(om, ces['CE2_m0_m1_exact_absD2_le_tau2']['delta'], 'CE2')]
    for r in lp:
        checks.append(shifted_power_check(om, r['argmax_delta'], f"LP_argmax_K{r['K']}"))
        checks.append(shifted_power_check(om, r['argmin_delta'], f"LP_argmin_K{r['K']}"))
    for c in committed:
        checks.append(shifted_power_check(om, np.array(c['delta']), f"committed_seed0_K{c['K']}"))
    # random null-space redistributions on random pole subsets (clusters)
    rng = np.random.default_rng(SEED_MC)
    worst = dict(ratio=0.0); nmc = 0; nviol = 0; per_s = {}
    for _ in range(20000):
        msub = int(rng.integers(2, N + 1))
        sub = np.sort(rng.choice(N, size=msub, replace=False))
        K = int(rng.integers(0, msub - 1))                 # preserve m0..mK, need null dim msub-K-1 >= 1
        Vand = np.array([om[sub] ** j for j in range(K + 1)])
        ns = null_space(Vand)
        dsub = ns @ rng.standard_normal(ns.shape[1])
        neg = dsub < 0
        alpha = rng.uniform(0.05, 1.0) * np.min(-w[sub][neg] / dsub[neg])
        d = np.zeros(N); d[sub] = alpha * dsub
        r = shifted_power_check(om, d, 'mc')
        nmc += 1
        if not r['bound_holds']:
            nviol += 1
        per_s.setdefault(r['s'], []).append(r['ratio_abs_Delta_over_bound'])
        if r['ratio_abs_Delta_over_bound'] > worst['ratio']:
            worst = dict(ratio=r['ratio_abs_Delta_over_bound'], s=r['s'], n_support=r['n_support'],
                         D=r['cluster_diameter_D'])
    mc = dict(n=nmc, seed=SEED_MC, n_violations=nviol, worst=worst,
              max_ratio_by_s={int(s): float(np.max(v)) for s, v in sorted(per_s.items())},
              count_by_s={int(s): len(v) for s, v in sorted(per_s.items())})
    # synthetic clusters: tightness and D^{s-1} scaling
    synth = []
    for s in (2, 3, 4):
        for Dm in (1.0, 0.1, 0.01):
            nodes = 7.0 + np.linspace(-Dm / 2, Dm / 2, s)
            ns = null_space(np.array([nodes ** j for j in range(s - 1)]))[:, 0]
            dd = ns / np.abs(ns).sum()                          # rho = 1
            Dsm1 = float(np.sum(dd * nodes ** (s - 1)))
            synth.append(dict(s=s, D=Dm, rho=1.0, abs_Delta_s_minus_1=abs(Dsm1),
                              bound=(Dm / 2) ** (s - 1), ratio=abs(Dsm1) / (Dm / 2) ** (s - 1)))
    thA = [theorem_a_check(om, ces['CE1_m0_m1_m2_exact']['delta'], 2, float(sv[-1]), 'CE1 (K=2)'),
           theorem_a_check(om, ces['CE2_m0_m1_exact_absD2_le_tau2']['delta'], 1, float(sv[-1]), 'CE2 (K=1)')]
    for r in lp:
        thA.append(theorem_a_check(om, r['argmax_delta'], r['K'], float(sv[-1]), f"LP_argmax_K{r['K']}"))
    M1 = dict(
        bound='if delta preserves m0..m_{s-2} and supp(delta) lies in a cluster of diameter D, midpoint c: '
              'Delta_{s-1} = sum_n delta_n (omega_n - c)^{s-1}, |Delta_{s-1}| <= rho (D/2)^{s-1}, rho=||delta||_1',
        cases=checks, all_cases_hold=bool(all(c['bound_holds'] for c in checks)),
        max_identity_rel_err=float(max(c['identity_rel_err'] for c in checks)),
        random_null_space_mc=mc, synthetic_clusters=synth,
        synthetic_note='ratio independent of D at fixed s => the bound captures the D^{s-1} scaling; '
                       'ratio = 1 at s=2 (bound tight)',
        theorem2a_sigma_min_check=thA,
        theorem2a_note='Theorem 2(a) uses ||delta||_2 and raw-node sigma_min; ||delta||_2 <= ||delta||_1 <= '
                       'sqrt(N) ||delta||_2 relates it to rho of part (b)')
    for c in checks[:2]:
        print(f"M1 {c['case']}: |Delta_{c['order_tested']}| = {abs(c['Delta_s_minus_1']):.4f} <= "
              f"rho(D/2)^{c['s']-1} = {c['bound_rho_Dhalf_pow']:.4f}  (ratio {c['ratio_abs_Delta_over_bound']:.3f})")
    print(f"M1 MC: {nmc} random cases, violations = {nviol}, worst ratio = {worst['ratio']:.4f}")

    rt = time.time() - t0
    record_theory_numerics('R3_within_sector_lp', R3, 'within_sector_lp.py', rt, seed=None)
    record_theory_numerics('M1_shifted_power_bound', M1, 'within_sector_lp.py', rt, seed=SEED_MC)


if __name__ == '__main__':
    main()
