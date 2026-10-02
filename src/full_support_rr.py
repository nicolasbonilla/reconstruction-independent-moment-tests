# -*- coding: utf-8 -*-
r"""full_support_rr.py -- Rayleigh-Ritz on the full ground-state support of the calibration ring.

System: the doped L=6 Hubbard ring of the deployed battery, (N_up, N_dn) = (2, 2), periodic boundary
conditions, current probe J, in the spectral_lanczos sector basis (the construction of verify.py,
interval_battery.py and within_sector_lp.py), at U/t = 3, 4, 5, 6, 8.

For each U the script takes S = supp(psi0) = {x : |psi0(x)|^2 > 1e-14} (the n_sup of interval_battery.py),
the Rayleigh-Ritz reconstruction on span S (H_S = Pi H Pi, its ground state |0_S>, e0, probe Pi J |0_S>), and
records:
  * moments k = 0..4 of the exact spectrum, m_k = <J psi0|(H-E0)^k|J psi0>, and of the reconstruction,
    mbar_k = <Pi phi_S|(H_S-e0)^k|Pi phi_S>, phi_S = J|0_S>; the signed residual Delta_k = m_k - mbar_k;
  * the same residual in the shared-state (leak) form of Eq. (app-leak),
    Delta_k^leak = <phi_S|(H-e0)^k|phi_S> - <Pi phi_S|(H_S-e0)^k|Pi phi_S>; at full support |0_S> = |0> and
    e0 = E0, so the two forms coincide (checked), and Delta_2 = ||(1-Pi) H phi_S||^2 (checked);
  * the leak score Delta_k/m_k and the ratio |Delta_k|/tau_k, with tau_k from the deployed rule
    tau_k = 1.96 sqrt(Var(O_k^loc)/5e4 + (0.02 m_k)^2) evaluated on that U's exact state (k = 0..4; m3, m4 are
    not in the deployed battery) and the deployed verdict (reject iff |Delta_k| > tau_k for some k <= 2);
  * the relative L1 error of the Lorentzian-broadened reconstruction, int |A_RR - A| / int A at eta = 0.15 t,
    on a wide grid ([-40, 80] t), on the grid of the audit's independent check ([-2, 25] t, 20001 points), and
    in a window-free form (numerator over [-400, 400] t, denominator m0); the windows differ by < 0.01;
  * the distinct poles (degenerate eigenstates aggregated, |omega - omega'| < 1e-6) and weights of both spectra;
  * the Wasserstein-1 distance between the normalized measures (exact, from the CDFs) and the bound
    sqrt(sigma^2 + sigma_RR^2 + (mean - mean_RR)^2) that any coupling (here the independent one) gives;
  * the one-node Krylov deficit m0 Var_A = m2 - m1^2/m0 (relative to m2), for the error class of Sec. blind.

Second part (author decision D4, optional computation of the audit): the deterministic shared-state leak
scores Delta_k^leak / mhat_k (k = 0, 1, 2) of every truncation d = 1..n_sup of the U/t = 4 support (the
deterministic lexsort order of interval_battery.py), next to the reference residuals |m_k - mbar_k|/m_k and the
relative L1 error, with Spearman rank correlations between leak score and L1. Truncations whose subspace
ground state is degenerate (gap < 1e-8) are flagged and left out of the correlations; 'spearman_truncations_only'
also leaves out d = n_sup (the full support). For each d the scan records the noise-free verdict of the default
battery (tau_k of the exact state) with the oracle (exact m_k) and with the shared-state estimator (tau_k as a
convention), and whether the leak score exceeds the numerical tolerance 1e-6.

Output: data/2026-10-02_full_support_rr.json.   Run: cd src && python full_support_rr.py
Dense linear algebra only (sector dimension 225); no random numbers (seed = none); runs in seconds.
"""
import time
T_START = time.time()
import os, sys, json, datetime, platform
import numpy as np
import scipy
import scipy.sparse as sp
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import spectral_lanczos as sl

DATA = os.path.normpath(os.path.join(HERE, '..', 'data'))
OUT_NAME = '2026-10-02_full_support_rr.json'
L = 6
U_LIST = (3.0, 4.0, 5.0, 6.0, 8.0)
NS, Z, BIAS = 50000, 1.96, 0.02
ETA = 0.15
KMAX = 4
SUPP_TOL = 1e-14                 # |psi0|^2 > SUPP_TOL defines supp(psi0) (as n_sup in interval_battery.py)
LOC_MASK = 1e-12                 # |psi0| > LOC_MASK for the local estimator (as interval_battery.py)
POLE_MERGE = 1e-6                # distinct-pole aggregation (as within_sector_lp.py)
W_MIN = 1e-10                    # weights below this are not listed as poles
LEAK_TOL = 1e-6                  # numerical tolerance on the shared-state leak score (Sec. 2.2(a))
ROUND = 12                       # |psi0|^2 key for the deterministic truncation order (interval_battery.py)
GRID_WIDE = np.linspace(-40.0, 80.0, 120001)     # d omega = 1e-3 t
GRID_AUDIT = np.linspace(-2.0, 25.0, 20001)
GRID_LINE = np.linspace(-400.0, 400.0, 800001)   # numerator only; denominator m0 (the whole real line)
_trapz = getattr(np, 'trapezoid', None) or np.trapz


# ------------------------------------------------------------------------------------------------
def ring(U):
    """L=6 doped ring sector (N_up,N_dn)=(2,2): H and current J, exactly as within_sector_lp.l6u4_system."""
    nup = nd = L // 3
    Tu, Su, _ = sl.hop(L, nup); Td, Sd, _ = sl.hop(L, nd)
    Du, Dd = len(Su), len(Sd)
    upocc = sl.occ_matrix(Su, L); dnocc = sl.occ_matrix(Sd, L)
    Hs = (sp.kron(Tu, sp.identity(Dd)) + sp.kron(sp.identity(Du), Td)
          + sp.diags(U * (upocc @ dnocc.T).ravel())).tocsr()
    Ju = sl.current_string(L, nup); Jd = sl.current_string(L, nd)
    Js = (sp.kron(Ju, sp.identity(Dd)) + sp.kron(sp.identity(Du), Jd)).tocsr()
    H = np.real(Hs.toarray()); J = Js.toarray()
    assert np.allclose(H, H.T) and np.allclose(J, J.conj().T)
    return H, J


def moments(H, phi, e, kmax=KMAX):
    """<phi|(H-e)^k|phi>, k=0..kmax, by repeated mat-vecs."""
    out = []; v = phi.copy()
    for _ in range(kmax + 1):
        out.append(float(np.vdot(phi, v).real))
        v = H @ v - e * v
    return out


def spectrum(H, phi, e):
    """Lehmann poles omega = E_a - e and weights |<a|phi>|^2 (all eigenstates)."""
    E, V = np.linalg.eigh(H)
    return E - e, np.abs(V.conj().T @ phi) ** 2


def distinct_poles(om, w):
    keep = w > W_MIN
    om, w = om[keep], w[keep]
    o = np.argsort(om, kind='stable'); om, w = om[o], w[o]
    uo, uw = [], []
    for a, b in zip(om, w):
        if uo and abs(a - uo[-1]) < POLE_MERGE:
            uw[-1] += b
        else:
            uo.append(float(a)); uw.append(float(b))
    return np.array(uo), np.array(uw)


def broadened(om, w, grid):
    return (w[None, :] * (ETA / np.pi) / ((grid[:, None] - om[None, :]) ** 2 + ETA ** 2)).sum(1)


def rel_l1(om_a, w_a, om_b, w_b, grid):
    """int |A_b - A_a| / int A_a on grid (A_a = exact)."""
    A = broadened(om_a, w_a, grid); B = broadened(om_b, w_b, grid)
    return float(_trapz(np.abs(B - A), grid) / _trapz(A, grid)), float(_trapz(A, grid) / w_a.sum())


def rel_l1_line(om_a, w_a, om_b, w_b):
    """window-free form: int |A_b - A_a| over [-400, 400] t (the |A_b - A_a| tail beyond decays as omega^-3)
    divided by m0 = int A_a over the real line."""
    A = broadened(om_a, w_a, GRID_LINE); B = broadened(om_b, w_b, GRID_LINE)
    return float(_trapz(np.abs(B - A), GRID_LINE) / w_a.sum())


def wasserstein1(om_a, w_a, om_b, w_b):
    """W_1 between the normalized discrete measures: int |F_a - F_b| d omega (exact, piecewise constant)."""
    pa = w_a / w_a.sum(); pb = w_b / w_b.sum()
    x = np.unique(np.concatenate([om_a, om_b]))
    Fa = np.array([pa[om_a <= xi].sum() for xi in x]); Fb = np.array([pb[om_b <= xi].sum() for xi in x])
    return float(np.sum(np.abs(Fa - Fb)[:-1] * np.diff(x)))


def deployed_taus(H, J, psi0, E0, kmax=KMAX):
    """tau_k = z sqrt(Var(O_k^loc)/N_s + (b m_k)^2), O_k = J (H-E0)^k J, local estimator over |psi0(x)|^2
    (interval_battery.py / within_sector_lp.deployed_taus)."""
    p = np.abs(psi0) ** 2; mask = np.abs(psi0) > LOC_MASK
    Jpsi = J @ psi0; v = Jpsi.copy(); out = []
    for k in range(kmax + 1):
        Okpsi = J @ v
        loc = np.zeros(len(psi0), dtype=complex); loc[mask] = Okpsi[mask] / psi0[mask]
        mean = float(np.sum(p * loc).real); var = float(np.sum(p * np.abs(loc) ** 2) - mean ** 2)
        out.append(dict(k=k, m_local_mean=mean, var_local=var,
                        tau=Z * np.sqrt(max(var, 0.0) / NS + (BIAS * abs(mean)) ** 2),
                        tau_shot_only=Z * np.sqrt(max(var, 0.0) / NS)))
        v = H @ v - E0 * v
    return out


def subspace_rr(H, J, idx):
    """Rayleigh-Ritz on span{x in idx}: returns e0, embedded |0_S>, gap_S, phi_S = J|0_S>, Pi phi_S."""
    idx = np.sort(np.asarray(idx))
    HS = H[np.ix_(idx, idx)]
    es, vs = np.linalg.eigh(HS)
    g = np.zeros(H.shape[0], dtype=complex); g[idx] = vs[:, 0]
    phi = J @ g
    Pphi = np.zeros_like(phi); Pphi[idx] = phi[idx]
    gap = float(es[1] - es[0]) if len(es) > 1 else float('inf')
    return float(es[0]), g, gap, phi, Pphi, HS, idx


# ------------------------------------------------------------------------------------------------
def one_U(U):
    H, J = ring(U)
    E, V = np.linalg.eigh(H)
    E0 = float(E[0]); psi0 = V[:, 0].astype(complex)
    gs_deg = int(np.sum(np.abs(E - E0) < 1e-8)); gap = float(E[gs_deg] - E0)
    p = np.abs(psi0) ** 2
    supp = np.flatnonzero(p > SUPP_TOL)
    off = np.setdiff1d(np.arange(len(p)), supp)
    phi0 = J @ psi0
    om, w = spectrum(H, phi0, E0)
    m = moments(H, phi0, E0)

    e0, g, gapS, phiS, PphiS, HS, idx = subspace_rr(H, J, supp)
    overlap = float(abs(np.vdot(psi0, g)))
    mhat = moments(H, phiS, e0)                          # full operator on the subspace state
    HSfull = np.zeros_like(H); HSfull[np.ix_(idx, idx)] = HS
    mbar = moments(HSfull, PphiS, e0)                    # reconstruction (confined to S)
    om_r, w_r = spectrum(HS, PphiS[idx], e0)

    taus = deployed_taus(H, J, psi0, E0)
    tau = [t['tau'] for t in taus]
    Delta = [m[k] - mbar[k] for k in range(KMAX + 1)]
    # without a bias budget the k=2 rule rejects once z sqrt(Var(O_2)/N_s) < |Delta_2|
    ns_reject_nobias_k2 = Z ** 2 * taus[2]['var_local'] / Delta[2] ** 2
    Delta_leak = [mhat[k] - mbar[k] for k in range(KMAX + 1)]
    # Delta_2 at unit captured weight: ||(1-Pi)(H-e0) phi_S||^2 = ||(1-Pi) H phi_S||^2 when phi_S in span S
    out_mask = np.ones(len(p), bool); out_mask[idx] = False
    Hphi = H @ phiS - e0 * phiS
    leak2_norm = float(np.vdot(Hphi[out_mask], Hphi[out_mask]).real)
    captured_w = float(np.vdot(PphiS, PphiS).real / np.vdot(phiS, phiS).real)

    hank = np.linalg.eigvalsh(np.array([[mbar[0], mbar[1]], [mbar[1], mbar[2]]]))[0]
    rejects = [k for k in range(3) if abs(Delta[k]) > tau[k]]
    verdict = 'REJECT' if (rejects or hank < -1e-9 or mbar[1] < -1e-9) else 'pass'

    po, pw = distinct_poles(om, w); pro, prw = distinct_poles(om_r, w_r)
    idom = int(np.argmax(pw)); idom_r = int(np.argmax(prw))
    l1w, capw = rel_l1(po, pw, pro, prw, GRID_WIDE)
    l1a, capa = rel_l1(po, pw, pro, prw, GRID_AUDIT)
    l1line = rel_l1_line(po, pw, pro, prw)
    W1 = wasserstein1(po, pw, pro, prw)
    mu = m[1] / m[0]; var = m[2] / m[0] - mu ** 2
    mu_r = mbar[1] / mbar[0]; var_r = mbar[2] / mbar[0] - mu_r ** 2
    w1_bound = float(np.sqrt(var + var_r + (mu - mu_r) ** 2))

    rec = dict(
        U_over_t=U, sector_dim=int(H.shape[0]), E0=E0, ground_degeneracy=gs_deg, ground_gap=gap,
        support_size=int(len(supp)), support_threshold_p=SUPP_TOL,
        min_p_on_support=float(p[supp].min()), max_p_off_support=float(p[off].max()) if len(off) else 0.0,
        rr_ground_energy_e0=e0, e0_minus_E0=e0 - E0, rr_subspace_gap=gapS, overlap_abs_psi0_0S=overlap,
        captured_weight_w=captured_w,
        m_exact_k0_to_4=m, mbar_rr_k0_to_4=mbar, mhat_shared_state_k0_to_4=mhat,
        Delta_signed_k0_to_4=Delta, Delta_leak_form_k0_to_4=Delta_leak,
        max_abs_diff_reference_vs_leak_form=float(max(abs(a - b) for a, b in zip(Delta, Delta_leak))),
        Delta2_as_norm_1minusPi_H_phi_sq=leak2_norm,
        abs_diff_Delta2_vs_norm_identity=abs(Delta_leak[2] - leak2_norm),
        leak_score_Delta_over_m_k0_to_4=[Delta[k] / m[k] for k in range(KMAX + 1)],
        tau_deployed_rule_k0_to_4=tau, tau_detail=taus,
        abs_Delta_over_tau_k0_to_4=[abs(Delta[k]) / tau[k] for k in range(KMAX + 1)],
        Ns_above_which_shot_only_rule_rejects_k2=ns_reject_nobias_k2,
        deployed_battery_m0_m1_m2_hankel=dict(verdict=verdict, orders_rejecting=rejects,
                                              hankel_min_eig_rr=float(hank)),
        orders_exceeding_same_rule_tau_k_le_4=[k for k in range(KMAX + 1) if abs(Delta[k]) > tau[k]],
        rel_L1_eta0p15_wide_grid=l1w, wide_grid=[float(GRID_WIDE[0]), float(GRID_WIDE[-1]), len(GRID_WIDE)],
        exact_weight_fraction_inside_wide_grid=capw,
        rel_L1_eta0p15_audit_grid=l1a, audit_grid=[float(GRID_AUDIT[0]), float(GRID_AUDIT[-1]), len(GRID_AUDIT)],
        exact_weight_fraction_inside_audit_grid=capa,
        rel_L1_eta0p15_full_line=l1line,
        full_line_definition='int_{-400}^{400} |A_RR - A| d omega / m0 (grid step 1e-3 t)',
        poles_exact=dict(omega=po, weight=pw), poles_rr=dict(omega=pro, weight=prw),
        n_distinct_poles_exact=int(len(po)), n_distinct_poles_rr=int(len(pro)),
        dominant_pole_exact=dict(omega=float(po[idom]), weight=float(pw[idom])),
        dominant_pole_rr=dict(omega=float(pro[idom_r]), weight=float(prw[idom_r])),
        lowest_pole_exact=dict(omega=float(po[0]), weight=float(pw[0])),
        lowest_pole_rr=dict(omega=float(pro[0]), weight=float(prw[0])),
        W1_normalized_exact_vs_rr=W1,
        mean_exact=mu, var_exact=var, mean_rr=mu_r, var_rr=var_r,
        W1_bound_independent_coupling=w1_bound,
        W1_bound_equal_first_two_moments_sqrt2_sigma=float(np.sqrt(2 * var)),
        one_node_krylov_m2_deficit=dict(abs=m[2] - m[1] ** 2 / m[0], rel_to_m2=(m[2] - m[1] ** 2 / m[0]) / m[2],
                                        rel_over_bias_budget=(m[2] - m[1] ** 2 / m[0]) / m[2] / BIAS,
                                        tau2_over_m2=tau[2] / m[2],
                                        abs_over_tau2=(m[2] - m[1] ** 2 / m[0]) / tau[2],
                                        passes_tau2=bool((m[2] - m[1] ** 2 / m[0]) <= tau[2]),
                                        note='one Lanczos node = one pole at m1/m0 with weight m0; '
                                             'mbar_2 = m1^2/m0, deficit m0 Var_A. It can exceed the 2% bias '
                                             'budget b and still pass, because tau_2 = 1.96 sqrt(shot^2 + (b m2)^2) '
                                             '>= 3.92% of m2'),
    )
    return rec, (H, J, psi0, E0, p, supp, po, pw, m)


def truncation_scan(H, J, psi0, E0, p, supp, po, pw, m):
    key = np.round(p, ROUND)
    order = np.lexsort((np.arange(len(p)), -key))[:len(supp)]
    tau = [t['tau'] for t in deployed_taus(H, J, psi0, E0, kmax=2)]
    rows = []
    for d in range(1, len(supp) + 1):
        e0, g, gapS, phiS, PphiS, HS, idx = subspace_rr(H, J, order[:d])
        mhat = moments(H, phiS, e0, kmax=2)
        HSfull = np.zeros_like(H); HSfull[np.ix_(idx, idx)] = HS
        mbar = moments(HSfull, PphiS, e0, kmax=2)
        om_r, w_r = spectrum(HS, PphiS[idx], e0)
        pro, prw = distinct_poles(om_r, w_r)
        l1 = rel_l1(po, pw, pro, prw, GRID_AUDIT)[0] if len(pro) else float('nan')
        hank = float(np.linalg.eigvalsh(np.array([[mbar[0], mbar[1]], [mbar[1], mbar[2]]]))[0])
        r_or = [abs(m[k] - mbar[k]) / tau[k] for k in range(3)]          # oracle: exact moments, no noise
        r_sh = [abs(mhat[k] - mbar[k]) / tau[k] for k in range(3)]       # shared state, tau_k as a convention
        lsc = [(mhat[k] - mbar[k]) / mhat[k] if mhat[k] != 0 else float('nan') for k in range(3)]
        rows.append(dict(d=d, coverage_d_over_support=d / len(supp),
                         degenerate_subspace_gs=bool(gapS < 1e-8), subspace_gap=gapS,
                         leak_score=lsc,
                         ref_residual_over_m=[abs(m[k] - mbar[k]) / m[k] for k in range(3)],
                         r_oracle_over_tau=r_or, r_shared_over_tau=r_sh,
                         oracle_passes_tau=bool(max(r_or) <= 1 and hank >= -1e-9 and mbar[1] >= -1e-9),
                         shared_passes_tau=bool(max(r_sh) <= 1 and hank >= -1e-9 and mbar[1] >= -1e-9),
                         leak_flagged_at_tolerance=bool(np.nanmax(np.abs(lsc)) > LEAK_TOL),
                         rel_L1=l1))
    good = [r for r in rows if not r['degenerate_subspace_gs'] and np.isfinite(r['rel_L1'])
            and all(np.isfinite(r['leak_score']))]
    L1v = np.array([r['rel_L1'] for r in good])
    corr = {}
    for k in range(3):
        rho, pv = spearmanr(np.array([r['leak_score'][k] for r in good]), L1v)
        rho2, pv2 = spearmanr(np.array([r['ref_residual_over_m'][k] for r in good]), L1v)
        corr[f'k{k}'] = dict(spearman_leak_score_vs_L1=float(rho), p_value=float(pv),
                             spearman_ref_residual_vs_L1=float(rho2), p_value_ref=float(pv2))
    trunc = [r for r in good if r['d'] < len(supp)]              # truncations only (d = n_sup is full support)
    L1t = np.array([r['rel_L1'] for r in trunc])
    corr_t = {}
    for k in range(3):
        rho, pv = spearmanr(np.array([r['leak_score'][k] for r in trunc]), L1t)
        corr_t[f'k{k}'] = dict(spearman_leak_score_vs_L1=float(rho), p_value=float(pv))

    def passes(key):
        sel = [r for r in rows if r['d'] < len(supp) and r[key]]
        l1 = [r['rel_L1'] for r in sel if r['rel_L1'] is not None and np.isfinite(r['rel_L1'])]
        cov = [r['coverage_d_over_support'] for r in sel]
        return dict(n=len(sel), d=[r['d'] for r in sel],
                    rel_L1_range=[float(min(l1)), float(max(l1))] if l1 else None,
                    coverage_range=[float(min(cov)), float(max(cov))] if cov else None)
    trs = [r for r in rows if r['d'] < len(supp)]
    verdicts = dict(
        tau_k=tau, tau_rule='deployed rule on the exact U/t=4 state (z=1.96, N_s=5e4, 2% bias)',
        n_truncations=len(trs),
        oracle_exact_moments_pass=passes('oracle_passes_tau'),
        shared_state_tau_convention_pass=passes('shared_passes_tau'),
        oracle_pass_rejected_by_shared_state=[r['d'] for r in trs
                                              if r['oracle_passes_tau'] and not r['shared_passes_tau']],
        leak_tolerance=LEAK_TOL,
        n_truncations_flagged_by_leak_score=int(sum(r['leak_flagged_at_tolerance'] for r in trs)),
        min_k2_leak_score_over_truncations=float(min(r['leak_score'][2] for r in trs
                                                      if np.isfinite(r['leak_score'][2]))),
        note='Noise-free. Oracle: |m_k - mbar_k| <= tau_k for k <= 2 with exact m_k (the Appendix-B scan). '
             'Shared state: |mhat_k - mbar_k| <= tau_k, mhat_k = <0_S|O_k|0_S> (tau_k as a convention). '
             'Leak score: |Delta_k/mhat_k| against the numerical tolerance LEAK_TOL.')
    return dict(order='np.lexsort((index, -round(|psi0|^2, 12)))[:n_sup] (interval_battery.py)',
                L1_grid='audit grid [-2,25] t, 20001 points, eta=0.15',
                n_rows=len(rows), n_degenerate_excluded=len(rows) - len(good), spearman=corr,
                n_truncations_in_spearman=len(trunc),
                d_range_in_spearman=[int(min(r['d'] for r in trunc)), int(max(r['d'] for r in trunc))],
                spearman_truncations_only=corr_t, verdicts=verdicts, rows=rows)


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return [_jsonable(v) for v in o.tolist()]
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    return o


def main():
    per_U = {}; keep4 = None
    for U in U_LIST:
        rec, aux = one_U(U)
        per_U[f'U{U:g}'] = rec
        if U == 4.0:
            keep4 = aux
        r = rec
        print(f"U/t={U:g}: E0={r['E0']:.6f} deg={r['ground_degeneracy']} gap={r['ground_gap']:.4f}  "
              f"supp={r['support_size']}/{r['sector_dim']}  |e0-E0|={abs(r['e0_minus_E0']):.1e}  "
              f"w={r['captured_weight_w']:.12f}")
        print("    Delta_k (k=0..4) = " + ", ".join(f"{x:.4g}" for x in r['Delta_signed_k0_to_4'])
              + "   |Delta|/tau = " + ", ".join(f"{x:.3g}" for x in r['abs_Delta_over_tau_k0_to_4']))
        print(f"    tau = {[round(x, 4) for x in r['tau_deployed_rule_k0_to_4']]}  verdict "
              f"{r['deployed_battery_m0_m1_m2_hankel']['verdict']}  relL1 wide={r['rel_L1_eta0p15_wide_grid']:.4f} "
              f"audit={r['rel_L1_eta0p15_audit_grid']:.4f} line={r['rel_L1_eta0p15_full_line']:.4f}  W1={r['W1_normalized_exact_vs_rr']:.4f} "
              f"(bound {r['W1_bound_independent_coupling']:.3f})")
        print(f"    dominant pole {r['dominant_pole_exact']['omega']:.4f} ({r['dominant_pole_exact']['weight']:.4f}) -> "
              f"{r['dominant_pole_rr']['omega']:.4f} ({r['dominant_pole_rr']['weight']:.4f}); lowest "
              f"{r['lowest_pole_exact']['omega']:.4f} -> {r['lowest_pole_rr']['omega']:.4f}; one-node deficit "
              f"{100*r['one_node_krylov_m2_deficit']['rel_to_m2']:.2f}% of m2")
    scan = truncation_scan(*keep4)
    v = scan['verdicts']
    print(f"scan U/t=4: oracle passes {v['oracle_exact_moments_pass']['n']} {v['oracle_exact_moments_pass']['d']}; "
          f"shared-state (tau convention) passes {v['shared_state_tau_convention_pass']['n']} "
          f"{v['shared_state_tau_convention_pass']['d']} L1 {v['shared_state_tau_convention_pass']['rel_L1_range']}; "
          f"flagged at leak tolerance {v['n_truncations_flagged_by_leak_score']}/{v['n_truncations']}; "
          f"oracle L1 {v['oracle_exact_moments_pass']['rel_L1_range']}")
    for k, v in scan['spearman_truncations_only'].items():
        print(f"scan U/t=4 truncations only {k}: Spearman = {v['spearman_leak_score_vs_L1']:+.3f} "
              f"(n={scan['n_truncations_in_spearman']}, d={scan['d_range_in_spearman']})")
    for k, v in scan['spearman'].items():
        print(f"scan U/t=4 {k}: Spearman(leak score, L1) = {v['spearman_leak_score_vs_L1']:+.3f}; "
              f"Spearman(ref residual, L1) = {v['spearman_ref_residual_vs_L1']:+.3f}  "
              f"(n={scan['n_rows'] - scan['n_degenerate_excluded']})")

    summary = dict(
        U_values=list(U_LIST),
        all_pass_deployed_battery=all(per_U[k]['deployed_battery_m0_m1_m2_hankel']['verdict'] == 'pass' for k in per_U),
        rel_L1_wide_by_U={k: per_U[k]['rel_L1_eta0p15_wide_grid'] for k in per_U},
        rel_L1_audit_grid_by_U={k: per_U[k]['rel_L1_eta0p15_audit_grid'] for k in per_U},
        rel_L1_full_line_by_U={k: per_U[k]['rel_L1_eta0p15_full_line'] for k in per_U},
        Delta2_by_U={k: per_U[k]['Delta_signed_k0_to_4'][2] for k in per_U},
        Delta2_over_tau2_by_U={k: per_U[k]['abs_Delta_over_tau_k0_to_4'][2] for k in per_U},
        Delta2_over_m2_by_U={k: per_U[k]['leak_score_Delta_over_m_k0_to_4'][2] for k in per_U},
        max_abs_Delta0_over_m0=max(abs(per_U[k]['leak_score_Delta_over_m_k0_to_4'][0]) for k in per_U),
        max_abs_Delta1_over_m1=max(abs(per_U[k]['leak_score_Delta_over_m_k0_to_4'][1]) for k in per_U),
        one_node_deficit_rel_m2_by_U={k: per_U[k]['one_node_krylov_m2_deficit']['rel_to_m2'] for k in per_U},
        catalog_U_3_to_6_one_node_deficit_range=[
            min(per_U[f'U{u:g}']['one_node_krylov_m2_deficit']['rel_to_m2'] for u in (3, 4, 5, 6)),
            max(per_U[f'U{u:g}']['one_node_krylov_m2_deficit']['rel_to_m2'] for u in (3, 4, 5, 6))],
        one_node_deficit_over_tau2_by_U={k: per_U[k]['one_node_krylov_m2_deficit']['abs_over_tau2'] for k in per_U},
        tau2_over_m2_by_U={k: per_U[k]['one_node_krylov_m2_deficit']['tau2_over_m2'] for k in per_U},
        Ns_shot_only_rule_rejects_k2_by_U={k: per_U[k]['Ns_above_which_shot_only_rule_rejects_k2'] for k in per_U},
    )
    doc = dict(
        _provenance=dict(
            script='src/full_support_rr.py', command='cd src && python full_support_rr.py',
            date=datetime.date.today().isoformat(), seed=None, runtime_s=round(time.time() - T_START, 2),
            python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
            platform=platform.platform(),
            audit_items='A2-01 (M5) full-support Rayleigh-Ritz; B2 one-node deficit on the ring; SOTA1-02 W_1; '
                        'D4 leak scores vs L1 across truncations',
            system='L=6 doped (2,2) Hubbard ring, PBC, current probe J (spectral_lanczos basis); dense ED',
            scope='Exact classical numerics; no sampling, no device data. tau_k is the deployed rule evaluated '
                  'on each U state (local-estimator variance, N_s=5e4, assumed 2% bias, z=1.96); the RR '
                  'residuals are deterministic.'),
        definitions=dict(
            Delta_k='m_k - mbar_k, m_k = <J psi0|(H-E0)^k|J psi0>, mbar_k = <Pi phi_S|(H_S-e0)^k|Pi phi_S>',
            Delta_k_leak='<phi_S|(H-e0)^k|phi_S> - <Pi phi_S|(H_S-e0)^k|Pi phi_S> (Eq. app-leak)',
            leak_score='Delta_k / m_k (deterministic; its tolerance is numerical, not shot noise)',
            rel_L1='int |A_RR - A| d omega / int A d omega, Lorentzian eta = 0.15 t, trapezoid on the stated grid',
            W1='Wasserstein-1 between the m0-normalized exact and RR measures (integral of |F - F_RR|)'),
        summary=summary, by_U=per_U, truncation_scan_U4=scan)
    path = os.path.join(DATA, OUT_NAME)
    tmp = path + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as fh:
        json.dump(_jsonable(doc), fh, indent=1, ensure_ascii=False)
    os.replace(tmp, path)
    print(f"wrote data/{OUT_NAME}  ({time.time() - T_START:.1f}s)")


if __name__ == '__main__':
    main()
