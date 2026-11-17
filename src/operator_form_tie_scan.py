# -*- coding: utf-8 -*-
"""Operator-form example (main text Sec. II C): robustness of the d=100 verdict pair to the tie order.

The d = 100 truncation of the doped L=6, U/t=4 ring (current probe) falls inside an exactly degenerate shell of
|psi0|^2. This script keeps every determinant above the shell and draws, 300 times (seed 20261108), which members of
the shell are kept, then evaluates on the subspace ground state the O_1 form <J(H-e0)J>, the double-commutator form
(1/2)<[J,[H,J]]> and the reconstruction moment mbar_1, with the same functions as estimator_form_offeigenstate.py,
against the exact-state default threshold tau_1 (within_sector_lp.deployed_taus). Writes
data/2026-11-08_operator_form_tie_scan.json. Post hoc (2026-11-08). Runs in seconds (5.9 s in the stored record).
"""
import os, sys, json, time, platform
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import within_sector_lp as wsl                      # noqa: E402
from estimator_form_offeigenstate import forms      # noqa: E402

D_CUT, NDRAW, SEED = 100, 300, 20261108


def main():
    t0 = time.time()
    Sy = wsl.l6u4_system()
    H, J, psi, E0, D = Sy['H'], Sy['J'], Sy['psi0'], Sy['E0'], Sy['dim']
    taus = wsl.deployed_taus(Sy, kmax=1)
    tau0, tau1 = taus[0]['tau_bias2pct'], taus[1]['tau_bias2pct']
    m1_exact = float(np.real(psi.conj() @ J.conj().T @ (H - E0 * np.eye(D)) @ J @ psi))
    p = np.round(np.abs(psi) ** 2, 12)
    order = np.lexsort((np.arange(D), -p))
    pc = p[order[D_CUT - 1]]
    above = np.where(p > pc)[0]; shell = np.where(p == pc)[0]
    n_keep = D_CUT - len(above)
    rng = np.random.default_rng(SEED)
    rows = []
    for i in range(NDRAW):
        S = np.sort(np.concatenate([above, rng.choice(shell, n_keep, replace=False)]))
        A, B, C, mb, e0, ex, mh0, mb0 = forms(H, J, S, D)
        rows.append(dict(A=A, B=B, mbar1=mb, deficit=(m1_exact - mb) / m1_exact,
                         A_rejects=bool(abs(A - mb) > tau1), B_rejects=bool(abs(B - mb) > tau1),
                         D0_over_tau0=abs(mh0 - mb0) / tau0))
    S = order[:D_CUT]
    A, B, C, mb, e0, ex, mh0, mb0 = forms(H, J, S, D)
    det = dict(A=A, B=B, mbar1=mb, leak_A_minus_mbar=A - mb, error_m1_minus_mbar=m1_exact - mb,
               deficit=(m1_exact - mb) / m1_exact, D0_over_tau0=abs(mh0 - mb0) / tau0)
    dfc = np.array([r['deficit'] for r in rows])
    out = dict(_what='tie-order robustness of the d=100 operator-form verdicts', _date='2026-11-08',
               _script='src/operator_form_tie_scan.py', python=platform.python_version(),
               m1_exact=m1_exact, tau1=tau1, tau0=tau0, d=D_CUT, n_above_shell=int(len(above)),
               shell_size=int(len(shell)), shell_members_kept=int(n_keep), n_draws=NDRAW, seed=SEED,
               deterministic_lexsort=det,
               draws_A_rejects=int(sum(r['A_rejects'] for r in rows)),
               draws_B_rejects=int(sum(r['B_rejects'] for r in rows)),
               draws_pair_A_rejects_B_passes=int(sum(r['A_rejects'] and not r['B_rejects'] for r in rows)),
               deficit_min=float(dfc.min()), deficit_max=float(dfc.max()),
               D0_over_tau0_min=float(min(r['D0_over_tau0'] for r in rows)),
               runtime_s=round(time.time() - t0, 1), draws=rows)
    path = os.path.normpath(os.path.join(HERE, '..', 'data', '2026-11-08_operator_form_tie_scan.json'))
    json.dump(out, open(path, 'w'), indent=1)
    print({k: v for k, v in out.items() if k != 'draws'})


if __name__ == '__main__':
    main()
