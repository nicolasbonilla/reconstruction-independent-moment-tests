# -*- coding: utf-8 -*-
"""Check that numbers the present manuscript prints from the two 2026-11-08 analyses match their JSON records.
Records (deposit, from version 2.0.3): ../data/2026-11-08_krylov_detection_order.json and
../data/2026-11-08_operator_form_tie_scan.json. Text: main.tex + supplement.tex of this folder.
42 of its 43 checks rebuild a printed string from one of the two records and require it verbatim in the text; the
remaining one reads no text: it checks that the tie-scan record counts the same two verdicts (the O_1 form rejects,
the commutator form does not) in all of its 300 draws. Runs in seconds."""
import os, json, sys, math
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, '..', 'data'))
PAPER = HERE
K = json.load(open(os.path.join(DATA, '2026-11-08_krylov_detection_order.json')))
T = json.load(open(os.path.join(DATA, '2026-11-08_operator_form_tie_scan.json')))
tex = ''.join(open(os.path.join(PAPER, f), encoding='utf-8').read() for f in ('main.tex', 'supplement.tex'))
if os.environ.get('CHECK_TEX_OVERRIDE'):            # optional: check the text of the file named there instead
    tex = open(os.environ['CHECK_TEX_OVERRIDE'], encoding='utf-8').read()
R = K['rows']
by = lambda n: [r for r in R if r['n_nodes'] == n]
fails = 0


def need(printed, why):
    """PASS only if the sentence fragment appears and EVERY occurrence of the same fragment with any digits carries
    exactly these digits (a copy with a changed number elsewhere in the text fails the check)."""
    global fails
    import re
    # generalize the printed values only: not subscripts, exponents or constants after '_', '^', '{', '-', '(' or a
    # letter, and not a leading table-row label such as '$1$ & '
    esc = re.escape(printed)
    lead = re.match(r'(\\\$\d\\\$\\ \\&\\ ){1,2}', esc)
    keep = lead.end() if lead else 0
    gen = esc[:keep] + re.sub(r'(?<![_^{\-(A-Za-z\\])\d+(?:\\\.\d+)?', lambda m: r'\d+(?:\.\d+)?', esc[keep:])
    found = re.findall(gen, tex)
    ok = printed in tex and all(f == printed for f in found)
    fails += not ok
    bad = [f for f in found if f != printed]
    print(f"[{'PASS' if ok else 'FAIL'}] {why}: {printed}" + (f"   (other occurrences: {bad[:3]})" if bad else ''))


def rng(v, f):
    return f % min(v), f % max(v)


a, b = rng([100 * r['rel_residual_2n'] for r in by(1)], '%.1f'); need(f"${a}$--${b}\\%$ of $m_2$", 'n=1 m2 residual')
a, b = rng([r['ratio_default'] for r in by(1)], '%.2f'); need(f"(${a}$--${b}\\,\\tau_2$)", 'n=1 D2/tau2')
a, b = rng([100 * r['rel_residual_by_k'][3] for r in by(1)], '%.1f'); need(f"${a}$--${b}\\%$, exceeds", 'n=1 m3 residual')
a, b = rng([r['ratio_default_by_k'][3] for r in by(1)], '%.2f'); need(f"(${a}$--${b}\\,\\tau_3$)", 'n=1 D3/tau3')
a, b = rng([100 * r['rel_residual_2n'] for r in by(2)], '%.2f'); need(f"${a}$--${b}\\%$ of $m_4$", 'n=2 m4 residual')
a, b = rng([r['ratio_default'] for r in by(2)], '%.3f'); need(f"(${a}$--${b}\\,\\tau_4$)", 'n=2 D4/tau4')
ns = [r['Ns_expose_nobias'] for r in by(2)]
need(f"${min(ns)/1e7:.1f}\\times10^7$ to ${max(ns)/1e8:.1f}\\times10^8$ shots", 'n=2 shots')
for n, k, sc in ((3, 6, 1e-5), (4, 8, 1e-6)):
    a, b = rng([r['rel_residual_2n'] / sc for r in by(n)], '%.1f')
    need(f"${a}$--${b}\\times10^{{{int(round(math.log10(sc)))}}}$ of $m_{k}$", f'n={n} residual')
taus = [x for r in R for x in r['tau_rel_default_by_k'][1:]]
# "for every k" encloses all 160 values, so its ends round outward: 4.852 -> 4.8, 7.223 -> 7.3
lo_, hi_ = math.floor(1000 * min(taus) + 1e-9) / 10, math.ceil(1000 * max(taus) - 1e-9) / 10
need(f"between ${lo_:.1f}\\%$ and ${hi_:.1f}\\%$ for every", 'threshold range (SM, outward)')
# upper bounds ('at most', 'to X relative') are rounded UP: 2.40003 -> 2.5, 1.329 -> 1.4
ceil1 = lambda x: math.ceil(10 * x - 1e-9) / 10
need(f"at most ${ceil1(max(r['max_rel_residual_below_2n'] for r in R)/1e-15):.1f}\\times10^{{-15}}$", 'below-2n')
need(f"to ${ceil1(max(r['gauss_identity_relerr'] for r in R)/1e-9):.1f}\\times10^{{-9}}$ relative", 'Gauss identity')
# Corollary 2 (Chebyshev--Jensen bound)
assert all(r['bound_holds'] for r in R), 'bound violated'
need(f"so $D/(4\\bar\\omega)={sum(r['D_over_4wbar'] for r in R)/len(R):.2f}$", 'D/(4 wbar) on the ring')
f = [r['bound_over_actual'] for r in R]; need(f"loose by factors of ${min(f):.1f}$ to ${max(f):.0f}$", 'bound looseness')
need(f"all ${len(R)}$ $(U,n)$ pairs", 'number of pairs')
assert all(r['bias_floor_excludes_order_2n'] == (r['n_nodes'] >= 2) for r in R), 'bias-floor claim'
RE = {1: -2, 2: -3, 3: -5, 4: -6}                     # exponent of the r_2n cell of Table III


def rcell(n):
    a, b = rng([r['rel_residual_2n'] / 10.0 ** RE[n] for r in by(n)], '%.1f')
    return f"$({a}$--${b})\\times10^{{{RE[n]}}}$"


for n in (1, 2):
    a, b = rng([r['cheb_jensen_bound_rel'] for r in by(n)], '%.2f' if n == 1 else '%.3f')
    need(f"{rcell(n)} & ${a}$--${b}$", f'Table III r_2n and bound n={n}')
for n, sc in ((3, 1e-3), (4, 1e-5)):
    a, b = rng([r['cheb_jensen_bound_rel'] / sc for r in by(n)], '%.1f')
    need(f"{rcell(n)} & $({a}$--${b})\\times10^{{{int(round(math.log10(sc)))}}}$", f'Table III r_2n and bound n={n}')
for n in range(1, 5):
    a, b = rng([r['rel_L1_eta0p15_wide_grid'] for r in by(n)], '%.2f')
    need(f"${n}$ & ${a}$--${b}$", f'Table III row {n} L1')
# Table SII (supplement): every row, formatted exactly as export_krylov_table_tex.py formats it


def _sci(x, nd=1):
    if x is None:
        return '---'
    m, e = f"{x:.{nd}e}".split('e')
    return f"${m}\\times10^{{{int(e)}}}$"


def _ratio(x):
    return f"${x:.3f}$" if x >= 1e-3 else _sci(x)


for r in sorted(R, key=lambda r: (r['n_nodes'], r['U'])):
    kd, ks = r['first_order_exposed_default'], r['first_order_exposed_shot_only']
    row = ' & '.join([f"${r['n_nodes']}$", f"${r['U']:g}$", f"${r['rel_L1_eta0p15_wide_grid']:.2f}$",
                      _sci(r['rel_residual_2n'], 2), _sci(r['cheb_jensen_bound_rel'], 2), _ratio(r['ratio_default']),
                      _sci(r['Ns_expose_nobias']), _sci(r['Ns_power80_nobias']), _sci(r['bias_frac_max_for_exposure'], 2),
                      (f"${kd}$" if kd is not None else '---') + ' / ' + (f"${ks}$" if ks is not None else '---')]) + r' \\'
    need(row, f"Table SII row n={r['n_nodes']} U={r['U']:g}")
# operator-form tie scan
d = T['deterministic_lexsort']
need(f"its residual, ${d['leak_A_minus_mbar']:.2f}$, is the leak and exceeds the actual first-moment error, ${d['error_m1_minus_mbar']:.2f}$", 'leak vs error')
need(f"$\\bar m_1$ is ${100*T['deficit_min']:.1f}$--${100*T['deficit_max']:.1f}\\%$ low", 'tie-scan deficit range')
ok = T['draws_pair_A_rejects_B_passes'] == T['n_draws'] == 300
print(f"[{'PASS' if ok else 'FAIL'}] tie scan: O1 rejects and commutator passes in {T['draws_pair_A_rejects_B_passes']}/{T['n_draws']}")
fails += not ok
need(f"$\\bar m_1$ is ${100*d['deficit']:.0f}\\%$ below the exact value", 'deterministic deficit 11%')
print('RESULT:', 'ALL PASS' if fails == 0 else f'{fails} FAIL')
sys.exit(1 if fails else 0)
