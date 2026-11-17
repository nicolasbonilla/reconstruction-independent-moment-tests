# -*- coding: utf-8 -*-
"""Export the PGFPlots data of Fig. krylov from 2026-11-08_krylov_detection_order.json.
Writes ../paper/figs/krylov_order_n{1..4}.dat (k >= 2n, relative residual |m_k - mbar_k|/m_k at U/t = 4),
../paper/figs/krylov_order_bound.dat (the Chebyshev-Jensen bound on the relative residual at k = 2n, U/t = 4),
../paper/figs/krylov_order_tau.dat (tau_k/m_k, default rule and shot-only, k = 1..10) and prints
U/t = 3..6 ranges quoted in the text and in Table krylov."""
import os, json
HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.normpath(os.path.join(HERE, '..', 'paper', 'figs'))
d = json.load(open(os.path.join(HERE, '..', 'data', '2026-11-08_krylov_detection_order.json')))
rows = d['rows']
U0 = 4.0
for n in range(1, 5):
    r = next(x for x in rows if x['U'] == U0 and x['n_nodes'] == n)
    with open(os.path.join(FIGS, f'krylov_order_n{n}.dat'), 'w') as f:
        f.write('k rel\n')
        for k in range(2 * n, len(r['rel_residual_by_k'])):
            f.write(f"{k} {r['rel_residual_by_k'][k]:.6e}\n")
with open(os.path.join(FIGS, 'krylov_order_bound.dat'), 'w') as f:
    f.write('k bound\n')
    for n in range(1, 5):
        r = next(x for x in rows if x['U'] == U0 and x['n_nodes'] == n)
        f.write(f"{2 * n} {r['cheb_jensen_bound_rel']:.6e}\n")
r = next(x for x in rows if x['U'] == U0 and x['n_nodes'] == 1)
with open(os.path.join(FIGS, 'krylov_order_tau.dat'), 'w') as f:
    f.write('k taudef taushot\n')
    for k in range(1, len(r['tau_rel_default_by_k'])):
        f.write(f"{k} {r['tau_rel_default_by_k'][k]:.6e} {r['tau_rel_shot_by_k'][k]:.6e}\n")
for n in range(1, 5):
    sel = [x for x in rows if x['n_nodes'] == n]
    def rng(key, fmt):
        v = [x[key] for x in sel if x[key] is not None]
        return (fmt % min(v)) + '--' + (fmt % max(v)) if v else 'none'
    print(f"n={n}: rel D_2n {rng('rel_residual_2n','%.2e')}  D_2n/tau {rng('ratio_default','%.3f')}  "
          f"L1 {rng('rel_L1_eta0p15_wide_grid','%.2f')}  Ns_expose(b=0) {rng('Ns_expose_nobias','%.1e')}  "
          f"Ns_80%(b=0) {rng('Ns_power80_nobias','%.1e')}  bias_max {rng('bias_frac_max_for_exposure','%.2e')}  "
          f"first exposed default {[x['first_order_exposed_default'] for x in sel]}  shot-only "
          f"{[x['first_order_exposed_shot_only'] for x in sel]}  gauss relerr max "
          f"{max(x['gauss_identity_relerr'] for x in sel):.1e}  below-2n max {max(x['max_rel_residual_below_2n'] for x in sel):.1e}")
    if n == 1:
        print('   n=1 D_3/tau_3 by U:', [round(x['ratio_default_by_k'][3], 2) for x in sel],
              ' rel D_3:', [round(x['rel_residual_by_k'][3], 4) for x in sel])
    if n == 2:
        print('   n=2 D_10/tau_10 by U:', [round(x['ratio_default_by_k'][10], 2) for x in sel],
              ' D_8/tau_8:', [round(x['ratio_default_by_k'][8], 2) for x in sel])
print('tau_k/m_k default (U=4):', [round(x, 4) for x in r['tau_rel_default_by_k']])
print('tau_k/m_k shot-only (U=4):', [round(x, 4) for x in r['tau_rel_shot_by_k']])
print('L1 at U=4 by n:', [round(next(x for x in rows if x['U'] == 4.0 and x['n_nodes'] == n)['rel_L1_eta0p15_wide_grid'], 3) for n in range(1, 5)])
