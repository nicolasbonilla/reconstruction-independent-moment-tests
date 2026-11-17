# -*- coding: utf-8 -*-
"""Write ../paper/figs/krylov_table.tex (rows of Table S-krylovfull) from 2026-11-08_krylov_detection_order.json."""
import os, json
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, '..', 'data', '2026-11-08_krylov_detection_order.json')))


def sci(x, nd=1):
    if x is None:
        return '---'
    m, e = f"{x:.{nd}e}".split('e')
    return f"${m}\\times10^{{{int(e)}}}$"


def ratio(x):
    return f"${x:.3f}$" if x >= 1e-3 else sci(x)


lines = []
for r in sorted(d['rows'], key=lambda r: (r['n_nodes'], r['U'])):
    k_def = r['first_order_exposed_default']; k_shot = r['first_order_exposed_shot_only']
    lines.append(' & '.join([
        f"${r['n_nodes']}$", f"${r['U']:g}$", f"${r['rel_L1_eta0p15_wide_grid']:.2f}$",
        sci(r['rel_residual_2n'], 2), sci(r['cheb_jensen_bound_rel'], 2), ratio(r['ratio_default']),
        sci(r['Ns_expose_nobias']), sci(r['Ns_power80_nobias']),
        sci(r['bias_frac_max_for_exposure'], 2),
        (f"${k_def}$" if k_def is not None else '---') + ' / ' + (f"${k_shot}$" if k_shot is not None else '---'),
    ]) + r' \\')
open(os.path.join(HERE, '..', 'paper', 'figs', 'krylov_table.tex'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
