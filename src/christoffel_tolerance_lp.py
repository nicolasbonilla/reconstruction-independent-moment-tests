# -*- coding: utf-8 -*-
r"""christoffel_tolerance_lp.py -- R4: the tolerance-inflated extremal atom (port of p3math/robust.py).

The Christoffel/Chebyshev-Markov-Stieltjes bracket |F_rec(t) - F_true(t)| <= W_n(t) assumes the reconstruction
matches m_0..m_2n EXACTLY. A reconstruction that PASSES the deployed battery only satisfies |Delta_k| <= tau_k.
The relevant bracket width is then the largest atom a positive measure with moments inside the tolerance box
can place at t:
    A_tau(t) = max { mu({t}) : mu >= 0, |m_k(mu) - m_k| <= tau_k, k = 0,1,2 }      (deployed order n = 1)
computed two ways:
  * LP over positive measures on a fine frequency grid [0, Omega] (+ the point t) -- robust.py;
  * closed form over all of R (Hamburger): A_tau(t) = max a such that (M0-a)(M2-a t^2) >= dist(a t, [m1-tau1,
    m1+tau1])^2 with M0 = m0+tau0, M2 = m2+tau2 (the residual measure exists iff its 2x2 Hankel is PSD; the box
    optimum sits at the largest m0', m2' and the m1' nearest a t). With tau = 0 this is W_1(t).
Two measures:
  U4  the L=6, U/t=4 within-sector current measure (poles 2.73..10.66) with the deployed tau_k (rule of
      interval_battery.py; within_sector_lp.deployed_taus) -- the numbers the plan quotes (robust.py);
  U8  the L=6, U/t=8 doped-current measure of the Christoffel/Markov-Krein figure (export_christoffel_dat.py,
      markov_krein_window.current_measure), built here in its (N=4, N_up=2) sector with hubbard_ed so that the
      same rule gives its tau_k (local-estimator variance at N_s=5e4, 2% bias). Its moments are checked against
      current_measure() (normalized). A second tau variant uses the U4 relative tolerances tau_k/m_k.
Writes key 'R4_tolerance_inflated_atom' of data/2026-09-27_theory_numerics.json. Deterministic (no RNG).
Run: cd src && python christoffel_tolerance_lp.py
"""
import os, sys, time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from small_checks import record_theory_numerics
import within_sector_lp as wsl

NS, Z, BIAS = 50000, 1.96, 0.02


def W1(m, t):
    mu = m[1] / m[0]; var = m[2] / m[0] - mu ** 2
    return m[0] / (1.0 + (t - mu) ** 2 / var)


def atom_closed_form(m, tau, t):
    """sup of the atom at t over positive measures on R with |m_k' - m_k| <= tau_k (k=0,1,2)."""
    M0, M2 = m[0] + tau[0], m[2] + tau[2]
    lo1, hi1 = m[1] - tau[1], m[1] + tau[1]

    def feas(a):
        if a > M0 or M2 - a * t * t < 0:
            return False
        at = a * t
        dist = 0.0 if lo1 <= at <= hi1 else min(abs(at - lo1), abs(at - hi1))
        return (M0 - a) * (M2 - a * t * t) >= dist ** 2
    lo, hi = 0.0, M0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if feas(mid):
            lo = mid
        else:
            hi = mid
    return lo


def atom_lp(m, tau, t, omega_max, step=0.01):
    g = np.union1d(np.arange(0.0, omega_max + step / 2, step), [t])
    it = int(np.argmin(np.abs(g - t)))
    A = np.array([g ** k for k in range(3)])
    c = np.zeros(len(g)); c[it] = -1.0
    if np.all(np.asarray(tau) == 0):
        r = linprog(c, A_eq=A, b_eq=np.asarray(m[:3]), bounds=[(0, None)] * len(g), method='highs')
    else:
        Aub = np.vstack([A, -A]); bub = np.concatenate([np.asarray(m[:3]) + tau, -(np.asarray(m[:3]) - tau)])
        r = linprog(c, A_ub=Aub, b_ub=bub, bounds=[(0, None)] * len(g), method='highs')
    return float(-r.fun) if r.status == 0 else float('nan'), len(g)


def rows_for(m, tau, tlist, omega_max):
    rows = []
    for lab, t in tlist:
        w1 = W1(m, t)
        lp_exact, ng = atom_lp(m, [0, 0, 0], t, omega_max)
        lp_tau, _ = atom_lp(m, tau, t, omega_max)
        cf_tau = atom_closed_form(m, tau, t)
        rows.append(dict(label=lab, t_omega=float(t), W1_closed=w1, W1_over_m0=w1 / m[0],
                         lp_exact_moment_atom=lp_exact, inflated_atom_lp=lp_tau, inflated_atom_closed_R=cf_tau,
                         inflated_atom_over_m0=cf_tau / m[0], ratio_inflated_over_W1=cf_tau / w1,
                         lp_vs_closed_abs_diff=abs(lp_tau - cf_tau), grid_points=ng, grid=f'[0,{omega_max}] step 0.01'))
    return rows


def u8_sector_measure_and_taus():
    """L=6, U/t=8 doped current, (N=4, N_up=2) sector of the hubbard_ed Fock operators used by
    markov_krein_window.current_measure (the triplet ground state's M=0 member)."""
    import hubbard_ed as H
    from run_sumrule_falsifier import current_operator
    L, U = 6, 8.0
    c = H.build_operators(L); cd = [ci.getH() for ci in c]
    Ham = H.hubbard_hamiltonian(L, t=1.0, U=U, pbc=True, c=c).tocsr()
    Nd = np.real(sum((cd[q] @ c[q]).diagonal() for q in range(2 * L)))
    Nu = np.real(sum((cd[q] @ c[q]).diagonal() for q in range(0, 2 * L, 2)))
    J = current_operator(L, c).tocsr()
    idx = np.where((np.abs(Nd - 4) < 1e-9) & (np.abs(Nu - 2) < 1e-9))[0]
    Hs = np.real(Ham[idx][:, idx].toarray()); Js = J[idx][:, idx].toarray()
    E, V = np.linalg.eigh(Hs); E0 = E[0]; psi0 = V[:, 0]
    om = E - E0; w = np.abs(V.conj().T @ (Js @ psi0)) ** 2
    k = (om > 1e-6) & (w > 1e-12); om, w = om[k], w[k]
    o = np.argsort(om, kind='stable'); om, w = om[o], w[o]
    uo, uw = [], []
    for a, b in zip(om, w):
        if uo and abs(a - uo[-1]) < 1e-6:
            uw[-1] += b
        else:
            uo.append(float(a)); uw.append(float(b))
    om, w = np.array(uo), np.array(uw)
    # rule taus on this sector state
    p = np.abs(psi0) ** 2; mask = np.abs(psi0) > 1e-12
    Hc = Hs - E0 * np.eye(len(Hs)); Jpsi = Js @ psi0; v = Jpsi.copy(); taus = []
    for kk in range(3):
        Ok = Js @ v
        loc = np.zeros(len(psi0), dtype=complex); loc[mask] = Ok[mask] / psi0[mask]
        mean = float(np.sum(p * loc).real); var = float(np.sum(p * np.abs(loc) ** 2) - mean ** 2)
        taus.append(dict(k=kk, m_exact=float(np.vdot(Jpsi, v).real), var_local=var,
                         tau_bias2pct=Z * np.sqrt(max(var, 0) / NS + (BIAS * abs(mean)) ** 2)))
        v = Hc @ v
    return om, w, taus, dict(E0=float(E0), gap=float(E[1] - E[0]), sector_dim=len(idx))


def main():
    t0 = time.time()
    # ---------------- U/t = 4 within-sector measure ----------------
    S = wsl.l6u4_system(); om4, w4, _ = wsl.lehmann_measure(S)
    m4 = [float(np.sum(w4 * om4 ** k)) for k in range(3)]
    tau4 = [t['tau_bias2pct'] for t in wsl.deployed_taus(S, kmax=2)]
    mu4 = m4[1] / m4[0]
    t4 = [('t=10', 10.0), ('t=12', 12.0), ('t=14', 14.0), ('t=15', 15.0), ('t=20', 20.0), ('t=25', 25.0),
          ('mean', mu4), ('dominant_pole', float(om4[np.argmax(w4)]))]
    rows4 = rows_for(m4, tau4, t4, omega_max=40.0)
    for r in rows4:
        print(f"U4 {r['label']:>13s} t={r['t_omega']:7.3f}: W1={r['W1_closed']:.4f}  LP(exact)={r['lp_exact_moment_atom']:.4f}  "
              f"inflated LP={r['inflated_atom_lp']:.4f} closed={r['inflated_atom_closed_R']:.4f}  "
              f"x{r['ratio_inflated_over_W1']:.2f}")
    # ---------------- U/t = 8 Christoffel measure ----------------
    om8, w8, taus8, info8 = u8_sector_measure_and_taus()
    m8 = [float(np.sum(w8 * om8 ** k)) for k in range(3)]
    tau8 = [t['tau_bias2pct'] for t in taus8]
    import markov_krein_window as mkw
    uom, wn = mkw.current_measure()
    mcm = [float(np.sum(wn * uom ** k)) for k in range(3)]
    check = [abs(m8[k] / m8[0] - mcm[k] / mcm[0]) for k in range(3)]
    lo_e, hi_e = float(uom.min()), float(uom.max())              # the figure's full-Fock rescaling range
    tmap = lambda x: lo_e + (x + 1) * (hi_e - lo_e) / 2
    mu8 = m8[1] / m8[0]
    t8 = [('t=12', 12.0), ('t=13', 13.0), ('t=14', 14.0), ('t=15', 15.0), ('t=20', 20.0),
          ('rescaled_0.3', tmap(0.3)), ('rescaled_0.6', tmap(0.6)), ('rescaled_0.9', tmap(0.9)),
          ('mean', mu8), ('dominant_pole', float(om8[np.argmax(w8)]))]
    rows8 = rows_for(m8, tau8, t8, omega_max=60.0)
    rel4 = [tau4[k] / m4[k] for k in range(3)]
    tau8_rel = [rel4[k] * m8[k] for k in range(3)]
    rows8b = rows_for(m8, tau8_rel, t8, omega_max=60.0)
    for r in rows8:
        print(f"U8 {r['label']:>13s} t={r['t_omega']:7.3f}: W1/m0={r['W1_over_m0']:.4f}  inflated/m0={r['inflated_atom_over_m0']:.4f}  "
              f"x{r['ratio_inflated_over_W1']:.2f}  (LP-closed {r['lp_vs_closed_abs_diff']:.1e})")
    rt = time.time() - t0
    R4 = dict(
        definition='A_tau(t) = max mu({t}) over positive measures with |m_k(mu)-m_k| <= tau_k, k=0,1,2 '
                   '(deployed order n=1); with tau=0 it equals the Christoffel W_1(t)',
        methods=['LP on a frequency grid [0,Omega] step 0.01 plus t (support omega >= 0)',
                 'closed form over R by bisection (Hamburger; 2x2 Hankel of the residual measure PSD)'],
        U4_within_sector=dict(measure='L=6, U/t=4 doped current (within_sector_lp.lehmann_measure), unnormalized',
                              moments_m0_m2=m4, tau=tau4, tau_rule='deployed (interval_battery.py) 2% bias, Ns=5e4',
                              rows=rows4),
        U8_christoffel_measure=dict(
            measure='L=6, U/t=8 doped current, (N=4,N_up=2) sector of hubbard_ed (Christoffel/Markov-Krein figure '
                    'measure); unnormalized; divide by m0 for the figure normalization',
            sector_info=info8, n_weighted_atoms=int(len(w8)), poles=om8, weights=w8, moments_m0_m2=m8,
            normalized_moment_match_vs_current_measure_abs=check,
            tau_rule_values=tau8, tau_rule='same rule on this state: 1.96 sqrt(Var(O_k^loc)/5e4 + (0.02 m_k)^2)',
            taus_rule_detail=taus8, rows_rule_tau=rows8,
            tau_relative_variant=tau8_rel, tau_relative_variant_def='tau_k/m_k copied from the U4 deployed battery',
            rows_relative_tau=rows8b,
            figure_rescaling_range_omega=[lo_e, hi_e]),
        note_mean_rows='at t = mean (and near it) W_1 = m0 is a supremum approached only by sending a vanishing mass '
                       'to infinity, so the bounded-grid LP [0,Omega] stays below the closed form over R there; away '
                       'from the mean the two agree to <= 1.4e-7',
        plan_expectation='U4: 0.196/0.098/0.042 at t=12/15/20, i.e. 3.6-5.7x W1 (robust.py, rounded taus)')
    record_theory_numerics('R4_tolerance_inflated_atom', R4, 'christoffel_tolerance_lp.py', rt, seed=None)


if __name__ == '__main__':
    main()
