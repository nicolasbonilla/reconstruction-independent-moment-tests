# -*- coding: utf-8 -*-
r"""estimator_form_offeigenstate.py -- R6 / M3: the first-moment estimator forms coincide only on an exact
eigenstate (port of p3math/forms_check.py, 2026-09-26).

System: L=6, U/t=4, (N_up,N_dn)=(2,2) doped ring, current probe J (spectral_lanczos basis; dim 225; the
ground state has 186 nonzero amplitudes). A truncated reconstruction keeps the d most probable determinants
S, and |x> is the ground state of H restricted to S (energy e0 = <x|H|x>). Compared on |x>:
    A  = <x| J (H - e0) J |x>                  (the form fixed as the deployed estimator)
    B  = 1/2 <x| [J,[H,J]] |x>                 (double commutator / f-sum form)
    C  = Re <x| J [H,J] |x>                    (adjoint form <J ad_H(J)>)
    mb = <x| J P_S (H - e0) P_S J |x>          (the reconstruction's own read-off m_bar_1)
and the exact m1 = <0|J(H-E0)J|0>. On an eigenstate A = B = C; off it they differ.
Truncation order: DETERMINISTIC np.lexsort((index, -round(p,12))) (R7 convention), and, for the record,
the np.argsort(-p) order used by the 2026-09-26 scratch (forms_check.py), which breaks ties in |psi|^2
arbitrarily. Degenerate |psi|^2 shells cut by d are reported.
Verdicts use the deployed tau_1 (within_sector_lp.deployed_taus; 2% bias, Ns=5e4): a form 'fires' if
|form - m_bar_1| > tau_1.
Writes key 'R6_estimator_forms_off_eigenstate' of data/2026-09-27_theory_numerics.json. Deterministic.
Run: cd src && python estimator_form_offeigenstate.py
"""
import os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from small_checks import record_theory_numerics
import within_sector_lp as wsl

DS = (200, 186, 150, 100, 75, 50, 30)


def forms(H, J, S, D):
    Hs = H[np.ix_(S, S)]
    es, vs = np.linalg.eigh(Hs)
    x = np.zeros(D, dtype=complex); x[S] = vs[:, 0]; e0 = float(es[0])
    I = np.eye(D)
    HJ = H @ J - J @ H                           # [H,J]
    A = float(np.real(x.conj() @ J.conj().T @ (H - e0 * I) @ J @ x))
    B = 0.5 * float(np.real(x.conj() @ (J @ HJ - HJ @ J) @ x))
    C = float(np.real(x.conj() @ (J @ HJ) @ x))
    P = np.zeros((D, D)); P[S, S] = 1.0
    mb = float(np.real(x.conj() @ J.conj().T @ P @ (H - e0 * I) @ P @ J @ x))
    e_x = float(np.real(x.conj() @ H @ x))
    return A, B, C, mb, e0, e_x


def main():
    t0 = time.time()
    Sy = wsl.l6u4_system()
    H = Sy['H']; J = Sy['J']; psi = Sy['psi0']; E0 = Sy['E0']; D = Sy['dim']
    tau1 = wsl.deployed_taus(Sy, kmax=1)[1]['tau_bias2pct']
    m1_exact = float(np.real(psi.conj() @ J.conj().T @ (H - E0 * np.eye(D)) @ J @ psi))
    p = np.abs(psi) ** 2; nsup = int(np.sum(p > 1e-14))
    orders = {'lexsort_deterministic': np.lexsort((np.arange(D), -np.round(p, 12))),
              'argsort_scratch_2026-09-26': np.argsort(-p)}
    res = {}
    for oname, order in orders.items():
        rows = []
        for d in DS:
            S = order[:d]
            A, B, C, mb, e0, ex = forms(H, J, S, D)
            pc = np.round(p[order[d - 1]], 12); pn = np.round(p[order[d]], 12) if d < D else -1.0
            shell = int(np.sum(np.round(p, 12) == pc)) if pc == pn else 0
            rows.append(dict(d=d, coverage_d_over_dim=d / D, coverage_d_over_support=min(d, nsup) / nsup,
                             A_J_Hme0_J=A, B_half_double_comm=B, C_J_adH_J=C, m_bar_1=mb, e0=e0,
                             e0_equals_expH=abs(e0 - ex) < 1e-10,
                             rel_diff_A_minus_B_over_A=(A - B) / A,
                             A_minus_mbar=A - mb, B_minus_mbar=B - mb, exact_minus_mbar=m1_exact - mb,
                             A_fires=bool(abs(A - mb) > tau1), B_fires=bool(abs(B - mb) > tau1),
                             exact_fires=bool(abs(m1_exact - mb) > tau1),
                             verdicts_A_B_disagree=bool((abs(A - mb) > tau1) != (abs(B - mb) > tau1)),
                             cut_inside_degenerate_shell=bool(pc == pn), shell_size=shell))
        res[oname] = rows
    same = all(abs(r1['A_J_Hme0_J'] - r2['A_J_Hme0_J']) < 1e-9 and abs(r1['B_half_double_comm'] - r2['B_half_double_comm']) < 1e-9
               for r1, r2 in zip(res['lexsort_deterministic'], res['argsort_scratch_2026-09-26']))
    rt = time.time() - t0
    for oname, rows in res.items():
        print(f"--- order: {oname} ---")
        print(f"{'d':>4} {'cov%':>6} {'A=J(H-e0)J':>11} {'B=1/2[J,[H,J]]':>15} {'mbar1':>8} {'A fires':>8} {'B fires':>8} shell")
        for r in rows:
            print(f"{r['d']:>4} {100*r['coverage_d_over_dim']:6.1f} {r['A_J_Hme0_J']:11.4f} {r['B_half_double_comm']:15.4f} "
                  f"{r['m_bar_1']:8.4f} {str(r['A_fires']):>8} {str(r['B_fires']):>8} {r['shell_size']}")
    print(f"exact m1 = {m1_exact:.6f}; tau1 = {tau1:.4f}; orders agree: {same}")
    R6 = dict(system='L=6, U/t=4 (2,2) ring current, dim 225, 186 nonzero GS amplitudes',
              m1_exact=m1_exact, tau1_deployed=tau1, n_support=nsup, dim=D,
              forms={'A': '<x|J(H-e0)J|x>, e0=<x|H|x>', 'B': '1/2<x|[J,[H,J]]|x>', 'C': 'Re<x|J[H,J]|x>',
                     'm_bar_1': '<x|J P_S (H-e0) P_S J|x>'},
              rows_by_order=res, orders_give_identical_numbers=bool(same),
              plan_expectation='coverage 88.9/66.7/44.4/33.3/22.2/13.3%: A 4.227/4.082/4.605/5.909/7.859/12.62; '
                               'B 4.227/4.048/3.707/4.336/5.492/8.000; m_bar_1 4.227/3.392/3.810/2.651/1.153/0')
    record_theory_numerics('R6_estimator_forms_off_eigenstate', R6, 'estimator_form_offeigenstate.py', rt, seed=None)


if __name__ == '__main__':
    main()
