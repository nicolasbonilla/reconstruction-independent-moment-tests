# -*- coding: utf-8 -*-
"""Generator of paper/figs/gauss_state.dat (fig_gausslaw), plan items R10/M15 (2026-09-27).

What the figure plots. The exact state-level Gauss-law violation <psi~| sum_x G_x^2 |psi~> of the U(1) quantum
link model (lattice Schwinger model: N=4 staggered sites, 3 spin-1/2 links, OBC, w=1, m=0.6; H, G_n from
run_gausslaw_falsifier.build) for the physical ground state psi0 mixed with one unphysical state phi_u,
    psi~ = normalize((1 - eps) psi0 + eps phi_u),   eps in {0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30}.
eps is an AMPLITUDE mixing parameter: the unphysical weight of psi~ is w(eps) = eps^2 / ((1-eps)^2 + eps^2)
(0.155 at eps = 0.3). Because psi0 is orthogonal to phi_u and sum G^2 psi0 = 0,
    <sum G^2>(eps) = g_u * w(eps),   g_u = <phi_u| sum G^2 |phi_u>,
which rises monotonically from 0 for 0 <= eps <= 1/2 and depends on phi_u ONLY through g_u.

Provenance and the selection rule. The committed .dat is the state-level column
(screen_vs_gaugebreaking[].gauss_violation) of the 2026-08-25 run of work/03_src/_superseded/
run_gauss_spectral_screen.py:36-114 (JSON in the work tree, not in this repository; its circular moment-screen
columns are not reproduced here). That script chose phi_u as "the unphysical eigenstate of H with the largest
|<k|M|psi0>|", M = sum_n (-1)^n n_n. M is gauge-invariant, so M psi0 lies entirely in the physical subspace and
EVERY such overlap is zero up to rounding (<= ~2e-15): the rule picks by floating-point noise among unphysical
eigenstates whose g_u ranges from ~1 to 10. On the machine that made the figure it picked a state with g_u = 2.
This generator therefore uses an explicit, platform-independent rule that yields the same curve:
    phi_u = lowest-energy eigenstate of H restricted to the charge-neutral sector (Q_tot = 0, the sector of psi0)
            with the MINIMAL Gauss-law violation sum G^2 = 2 (one site at G=+1, one at G=-1).
Every vector of that subspace has sum G^2 = 2 exactly (the G_n are diagonal in the computational basis), so the
curve 2 w(eps) does not depend on which vector is taken. The ported noise rule is still evaluated and recorded as a
diagnostic (its pick, g_u and overlaps), together with the spread of g_u over all unphysical eigenstates.

Outputs: data/2026-09-27_gauss_state.json and paper/figs/gauss_state.dat (same two-column format and precision
as before). No random numbers are drawn."""
import time
T_START = time.time()                  # runtime recorded in the provenance includes the imports
import os, sys, json, platform, collections
import numpy as np
import scipy, scipy.sparse as sp
from run_gausslaw_falsifier import build, Q, N, NL, op_on, W, M as MASS

EPS = [0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30]
G_TARGET = 2.0                         # minimal charge-conserving Gauss-law violation
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, '..', 'data'))
FIGS = os.path.normpath(os.path.join(HERE, '..', 'paper', 'figs'))
OUT_JSON = '2026-09-27_gauss_state.json'
OUT_DAT = 'gauss_state.dat'


def sweep(psi0, phi, G2):
    rows = []
    for eps in EPS:
        tilde = (1 - eps) * psi0 + eps * phi
        tilde = tilde / np.linalg.norm(tilde)
        gviol = float(np.real(np.vdot(tilde, G2 @ tilde)))
        rows.append({'eps': eps, 'gauss_violation': gviol,
                     'unphysical_weight': eps ** 2 / ((1 - eps) ** 2 + eps ** 2)})
    return rows


def main():
    H, G, q = build()
    Hd = H.toarray()
    Z = sp.csr_matrix((2 ** Q, 2 ** Q), dtype=complex)
    G2 = sum((g @ g for g in G), Z).tocsr()
    Qtot = sum(q, Z).tocsr()
    comm = max(float(sp.linalg.norm(H @ g - g @ H)) for g in G)            # gauge invariance [H, G_n] = 0
    g2_diag, q_diag = np.real(G2.diagonal()), np.real(Qtot.diagonal())
    offdiag = float(sp.linalg.norm(G2 - sp.diags(G2.diagonal())))         # G2 is diagonal in the basis

    E, V = np.linalg.eigh(Hd)
    g2_all = np.array([float(np.real(np.vdot(V[:, k], G2 @ V[:, k]))) for k in range(len(E))])
    phys = np.where(g2_all < 1e-6)[0]
    i0 = int(phys[0]); E0, psi0 = float(E[i0]), V[:, i0]
    q0 = float(np.real(np.vdot(psi0, Qtot @ psi0)))

    # --- explicit rule: lowest-energy state of H in the (Q_tot = 0, sum G^2 = 2) subspace ---
    sub = np.where((np.abs(g2_diag - G_TARGET) < 1e-9) & (np.abs(q_diag - round(q0)) < 1e-9))[0]
    Hsub = Hd[np.ix_(sub, sub)]
    es, vs = np.linalg.eigh(Hsub)
    phi_u = np.zeros(2 ** Q, dtype=complex); phi_u[sub] = vs[:, 0]
    g_u = float(np.real(np.vdot(phi_u, G2 @ phi_u)))
    leak = float(np.linalg.norm(Hd @ phi_u - es[0] * phi_u))               # is it an eigenstate of the full H?
    rows = sweep(psi0, phi_u, G2)
    for r in rows:
        r['closed_form'] = g_u * r['unphysical_weight']
        r['abs_diff_closed_form'] = abs(r['gauss_violation'] - r['closed_form'])

    # --- diagnostic: the ported rule (largest |<k|M|psi0>| over unphysical eigenstates) ---
    NUMq = sp.csr_matrix(np.array([[0, 0], [0, 1]], complex))
    Mop = sum((((-1) ** n) * op_on(n, NUMq) for n in range(N)), Z).tocsr()
    unphys = np.where(g2_all > 0.5)[0]
    ov = np.abs(V[:, unphys].conj().T @ (Mop @ psi0))
    k_noise = int(unphys[np.argmax(ov)])
    Mpsi = Mop @ psi0
    frac_M_physical = float(np.linalg.norm(V[:, phys].conj().T @ Mpsi) / np.linalg.norm(Mpsi))
    rows_noise = sweep(psi0, V[:, k_noise], G2)
    g_spread = collections.Counter(np.round(g2_all[unphys], 3).tolist())

    # --- comparison with the committed .dat (read before it is rewritten) ---
    dat_path = os.path.join(FIGS, OUT_DAT)
    old_text = open(dat_path).read() if os.path.isfile(dat_path) else None
    committed = ([{'eps': float(a), 'g2': float(b)} for a, b in
                  (l.split() for l in old_text.strip().splitlines()[1:])] if old_text else None)
    new_text = 'eps g2\n' + ''.join(f"{r['eps']:.4f} {r['gauss_violation']:.5f}\n" for r in rows)
    identical = (old_text is not None) and (old_text.replace('\r\n', '\n') == new_text)
    max_dev = (max(abs(c['g2'] - r['gauss_violation']) for c, r in zip(committed, rows))
               if committed and len(committed) == len(rows) else None)
    mono = all(b['gauss_violation'] > a['gauss_violation'] for a, b in zip(rows, rows[1:]))

    print("=== fig_gausslaw generator: exact <sum_x G_x^2> of a physical state mixed with an unphysical one ===")
    print(f"U(1) link model N={N}, {NL} links, OBC, w={W}, m={MASS};  max||[H,G_n]|| = {comm:.1e}; "
          f"||offdiag(sum G^2)|| = {offdiag:.1e}")
    print(f"physical GS: index {i0}, E0 = {E0:.6f}, <sum G^2> = {g2_all[i0]:.1e}, Q_tot = {q0:.1e}")
    print(f"explicit phi_u: lowest state of H in (Q_tot=0, sum G^2=2) [{len(sub)} basis states], E = {es[0]:.6f}, "
          f"g_u = {g_u:.12f}, ||H phi_u - E phi_u|| = {leak:.1e}")
    print(f"ported rule: max |<k|M|psi0>| over {len(unphys)} unphysical eigenstates = {ov.max():.1e} "
          f"(M psi0 is {100*frac_M_physical:.10f}% physical) -> noise pick index {k_noise}, g_u = {g2_all[k_noise]:.6f}")
    print(f"g_u over unphysical eigenstates (value: count): {dict(sorted(g_spread.items()))}")
    print(" eps   weight w   <sum G^2>   closed form 2w   |diff|    (ported-rule pick)")
    for r, rn in zip(rows, rows_noise):
        print(f" {r['eps']:.2f}  {r['unphysical_weight']:.5f}   {r['gauss_violation']:.8f}   {r['closed_form']:.8f}   "
              f"{r['abs_diff_closed_form']:.1e}   {rn['gauss_violation']:.8f}")
    print(f"monotone increasing in eps: {mono}")
    print(f"committed paper/figs/{OUT_DAT}: identical after regeneration (content, line endings aside): {identical}; "
          f"max |dev| from the 5-decimal file = {max_dev:.1e}")

    runtime = time.time() - T_START
    out = {'_provenance': {'script': 'src/export_gauss_state.py', 'date': '2026-09-27', 'plan_item': 'R10/M15',
                           'seed': 'none (no random numbers; exact diagonalization with numpy.linalg.eigh)',
                           'runtime_s': round(runtime, 3), 'python': sys.version.split()[0],
                           'numpy': np.__version__, 'scipy': scipy.__version__, 'platform': platform.platform(),
                           'ported_from': 'work/03_src/_superseded/run_gauss_spectral_screen.py:36-114 '
                                          '(state-level <sum G^2> column only; the circular moment columns are not '
                                          'reproduced)',
                           'figure': 'paper/figs/fig_gausslaw.tex plots paper/figs/gauss_state.dat (eps, g2)'},
           'model': {'name': 'U(1) quantum link model (lattice Schwinger)', 'N_sites': N, 'N_links': NL, 'bc': 'OBC',
                     'w': W, 'm': MASS, 'max_comm_H_G': comm, 'sumG2_offdiag_norm': offdiag},
           'physical_ground_state': {'index': i0, 'E0': E0, 'sum_G2': float(g2_all[i0]), 'Q_tot': q0},
           'eps_meaning': 'amplitude mixing parameter; unphysical weight w = eps^2/((1-eps)^2+eps^2)',
           'unphysical_admixture': {'rule': 'lowest-energy eigenstate of H restricted to the (Q_tot = 0, sum G^2 = 2) '
                                            'subspace: the minimal charge-conserving Gauss-law violation',
                                    'subspace_dim': int(len(sub)), 'energy': float(es[0]), 'g_u': g_u,
                                    'eigenstate_residual_full_H': leak},
           'closed_form': 'g_u * eps^2 / ((1-eps)^2 + eps^2)',
           'monotone_increasing': mono,
           'rows': rows,
           'ported_rule_diagnostic': {'rule': 'unphysical eigenstate of H with the largest |<k|M|psi0>|, '
                                              'M = sum_n (-1)^n n_n',
                                      'max_overlap': float(ov.max()), 'median_overlap': float(np.median(ov)),
                                      'fraction_of_M_psi0_in_physical_subspace': frac_M_physical,
                                      'verdict': 'ill-posed: all overlaps vanish by gauge invariance of M; '
                                                 'the pick is set by rounding noise',
                                      'pick_index': k_noise, 'pick_energy': float(E[k_noise]),
                                      'pick_g_u': float(g2_all[k_noise]),
                                      'rows': rows_noise,
                                      'g_u_over_unphysical_eigenstates': {str(k): v for k, v in sorted(g_spread.items())}},
           'committed_dat': committed, 'committed_dat_identical_content': identical,
           'committed_dat_max_abs_dev': max_dev}
    json.dump(out, open(os.path.join(DATA, OUT_JSON), 'w'), indent=2)
    with open(dat_path, 'w', newline='\r\n' if (old_text and '\r\n' in open(dat_path, newline='').read()) else '\n') as f:
        f.write(new_text)
    print(f"wrote data/{OUT_JSON} and paper/figs/{OUT_DAT}  (runtime {runtime:.2f} s)")


if __name__ == '__main__':
    main()
