# -*- coding: utf-8 -*-
"""JOINT INTERVAL-MOMENT BATTERY (m0,m1,m2 + Hankel/Stieltjes), $0 -- the remedy for a lone-moment miss.
Honest reframe baked in: this is a SHOT-BUDGET-INFORMED FORECAST computed on the exact state (device counts
unsaved), NOT a load-bearing hardware test; shot term = ground-state local-estimator variance (a LOWER BOUND
on true per-shot uncertainty); bias term = an ASSUMED residual (swept). REJECT = any |m_hat_k - m_bar_k| > delta_k
OR interval-Hankel/Stieltjes infeasibility. Prints, per focus d, the lone-m1 verdict and the battery verdict.

Truncation order (2026-09-27, plan item R7). The top-d determinants are ranked by |psi0|^2 with a DETERMINISTIC
tie-break, np.lexsort((index, -round(|psi0|^2, 12))): descending probability, then ascending sector index. The
ground-state probabilities fall in exactly degenerate shells (6/12/24-fold); the earlier np.argsort(prob)[::-1]
broke ties by floating-point noise, which is why the committed run (interval_battery_showcase.json, |g1| = 0.2250
at d=98, an m1-alone miss at 2% bias) and fresh clones (|g1| = 0.254, a reject) disagreed. The lexsort choice is
reproducible but is still a convention whenever the cut falls inside a shell, so this script also
  (i)   scans EVERY truncation d = 1..n_sup-1 (not only the focus list) for lone-m1 misses and says which the
        battery catches, and repeats verify.py's auto-scan rule on the deterministic order;
  (ii)  records the degenerate shell at every cut and the convention-free shell-boundary cuts;
  (iii) measures the tie-break sensitivity at the focus cuts, at d=97 and at every cut the whole battery passes:
        R random choices (seeded per cut) of which shell members
        are kept, with the fraction of choices on which m1 alone rejects.
Every verdict printed is computed in this run. Output: data/2026-09-27_interval_battery.json (the committed
interval_battery_showcase.json is left untouched as the record of the argsort run)."""
import time
T_START = time.time()                   # runtime recorded in the provenance includes the imports
import numpy as np, json, os, sys, platform
import scipy, scipy.sparse as sp
import spectral_lanczos as sl

L, U, Ns, z = 6, 4.0, 50000, 1.96
ROUND = 12                               # decimals used to identify exactly degenerate |psi0|^2 shells
TIE_SEED, TIE_R = 20260927, 300          # tie-break sensitivity: seed and number of random shell choices
BIASES = (0.0, 0.02, 0.04)
DEPLOYED_BIAS = 0.02
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, '..', 'data'))
OUT_NAME = '2026-09-27_interval_battery.json'
COMMITTED = 'interval_battery_showcase.json'   # argsort-ordered record (read only, for the comparison)

nup = nd = L // 3
Tu, Su, iu = sl.hop(L, nup); Td, Sd, idd = sl.hop(L, nd)
Du, Dd = len(Su), len(Sd)
upocc = sl.occ_matrix(Su, L); dnocc = sl.occ_matrix(Sd, L)
Hs = (sp.kron(Tu, sp.identity(Dd)) + sp.kron(sp.identity(Du), Td) + sp.diags(U*(upocc@dnocc.T).ravel())).tocsr()
Ju = sl.current_string(L, nup); Jd = sl.current_string(L, nd)
Js = (sp.kron(Ju, sp.identity(Dd)) + sp.kron(sp.identity(Du), Jd)).tocsr()
E0, psi0 = sl._sub_gs(Hs)
Jp = Js @ psi0
Hc = (Hs - E0*sp.identity(Hs.shape[0])).tocsr()          # H - E0

# independent moments m_k = <Jp|(H-E0)^k|Jp>, k=0,1,2  and their LOCAL-estimator variances (shot term, LOWER bound)
phi = Jp.copy()
Hphi = Hc @ phi
H2phi = Hc @ Hphi
m0_op = float(np.vdot(phi, phi).real)
m1_op = float(np.vdot(phi, Hphi).real)
m2_op = float(np.vdot(phi, H2phi).real)
p = np.abs(psi0)**2; mask = np.abs(psi0) > 1e-12


def locvar(Ok_psi):                                       # O_k psi0 -> local estimator O_k^loc(x), var over |c_x|^2
    loc = np.zeros(len(psi0), dtype=np.result_type(Ok_psi, psi0)); loc[mask] = Ok_psi[mask]/psi0[mask]
    loc = np.real(loc)                                    # the imaginary part is identically zero here
    mean = float(np.sum(p*loc)); return mean, float(np.sum(p*loc**2)-mean**2)


# O_0=J^2, O_1=J(H-E0)J, O_2=J(H-E0)^2 J
m0_loc, var0 = locvar(Js@(Jp))
m1_loc, var1 = locvar(Js@(Hc@Jp))
m2_loc, var2 = locvar(Js@(H2phi))


def delta(var, m, bias_frac):
    return z*np.sqrt(max(var, 0)/Ns + (bias_frac*abs(m))**2)


def hankel_ok(m0, m1, m2):
    # T=0 support omega>=0: Hankel H1=[[m0,m1],[m1,m2]]>=0  AND shifted Stieltjes [m1]>=0 (m1>=0) etc.
    H1 = np.array([[m0, m1], [m1, m2]])
    return (np.linalg.eigvalsh(H1)[0] >= -1e-9) and (m1 >= -1e-9)


DELTAS = {b: (delta(var0, m0_op, b), delta(var1, m1_op, b), delta(var2, m2_op, b)) for b in BIASES}

# --- deterministic truncation order (R7) ---
prob = np.abs(psi0)**2
prob_key = np.round(prob, ROUND)
order_all = np.lexsort((np.arange(len(prob)), -prob_key))
nsup = int((prob > 1e-14).sum()); order = order_all[:nsup]
# for the comparison only: the old, platform-dependent argsort order on THIS machine
order_argsort = np.argsort(prob)[::-1][:nsup]

# robustness of the key itself: an independent dense eigensolver must give the SAME lexsort order, and no
# |psi0|^2 may sit near a rounding boundary of the 12-decimal key (else solver noise could flip a shell)
_ed, _Vd = np.linalg.eigh(Hs.toarray())
prob_dense = np.abs(_Vd[:, 0])**2
order_dense = np.lexsort((np.arange(len(prob_dense)), -np.round(prob_dense, ROUND)))[:nsup]
_frac = prob[prob > 1e-14] * 10**ROUND
ORDER_CHECK = {'max_abs_p_eigsh_minus_dense': float(np.max(np.abs(prob - prob_dense))),
               'lexsort_order_identical_with_dense_eigh': bool(np.array_equal(order, order_dense)),
               'min_distance_to_rounding_boundary': float(np.min(np.abs(_frac - np.floor(_frac) - 0.5))) / 10**ROUND,
               'ground_state_gap': float(_ed[1] - _ed[0])}
print(f"order check: eigsh vs dense max|dp| = {ORDER_CHECK['max_abs_p_eigsh_minus_dense']:.1e}; lexsort order "
      f"identical: {ORDER_CHECK['lexsort_order_identical_with_dense_eigh']}; closest |psi0|^2 to a rounding boundary "
      f"of the key: {ORDER_CHECK['min_distance_to_rounding_boundary']:.1e}; GS gap {ORDER_CHECK['ground_state_gap']:.3f}")

# degenerate shells of the support, in descending |psi0|^2
shell_vals = np.unique(prob_key[order])[::-1]
shell_sizes = [int(np.sum(prob_key == v)) for v in shell_vals]
shell_bounds = [int(x) for x in np.cumsum(shell_sizes)]            # convention-free cuts


def shell_at_cut(d):
    pc = prob_key[order[d - 1]]
    size = int(np.sum(prob_key == pc)); above = int(np.sum(prob_key > pc))
    return {'p_cut': float(pc), 'shell_size': size, 'n_above_shell': above, 'n_kept_from_shell': int(d - above),
            'cut_inside_shell': bool(0 < d - above < size)}


def trunc_moments(idx):
    idx = np.sort(np.asarray(idx)); d = len(idx)
    Hsub = Hs[idx][:, idx]; Jsub = Js[idx][:, idx]
    e0, g = sl._sub_gs(Hsub); Hcs = (Hsub - e0*sp.identity(d)).tocsr(); Jg = Jsub @ g
    return float(np.vdot(Jg, Jg).real), float(np.vdot(Jg, Hcs @ Jg).real), float(np.vdot(Jg, Hcs @ (Hcs @ Jg)).real)


def verdicts(b0, b1, b2, bias_frac):
    d0, d1, d2 = DELTAS[bias_frac]
    g0, g1, g2 = abs(b0 - m0_op), abs(b1 - m1_op), abs(b2 - m2_op)
    m1v = 'REJECT' if g1 > d1 else 'corrob'
    bv = 'REJECT' if ((g0 > d0) or (g1 > d1) or (g2 > d2) or (not hankel_ok(b0, b1, b2))) else 'corrob'
    return g0, g1, g2, m1v, bv


focus_ds = [186, 167, 98, 88, 64, 42]     # the committed focus list (d=98: the committed m1-only miss; 186: full)

print("=== JOINT INTERVAL BATTERY (m0,m1,m2 + Hankel), honest reframe; deterministic lexsort truncation order ===")
print(f"independent: m0={m0_op:.4f} m1={m1_op:.4f} m2={m2_op:.4f}")
print(f"local-est var (shot LOWER bound): var0={var0:.3g} var1={var1:.3g} var2={var2:.3g}")
print(f"|psi0|^2 shells on the support (sizes, descending p): {shell_sizes}; convention-free cuts d = {shell_bounds}")
rows = []
for bias_frac in BIASES:
    d0, d1, d2 = DELTAS[bias_frac]
    print(f"\n--- bias={100*bias_frac:.0f}% : delta0={d0:.3f} delta1={d1:.3f} delta2={d2:.3f} ---")
    print(" d     cov%  |g0|    |g1|    |g2|    m1-verdict  BATTERY(m0,m1,m2)+Hankel   cut")
    for d in focus_ds:
        d = int(d)
        b0, b1, b2 = trunc_moments(order[:d])
        g0, g1, g2, m1v, bv = verdicts(b0, b1, b2, bias_frac)
        sh = shell_at_cut(d)
        tag = (f"inside {sh['shell_size']}-fold shell ({sh['n_kept_from_shell']} kept)" if sh['cut_inside_shell']
               else 'shell boundary')
        print(f" {d:5d} {100*d/nsup:5.1f}  {g0:6.3f}  {g1:6.3f}  {g2:6.3f}   {m1v:7s}     {bv:7s}                  {tag}")
        rows.append(dict(bias_frac=bias_frac, d=d, cov=d/nsup, g0=g0, g1=g1, g2=g2,
                         delta0=d0, delta1=d1, delta2=d2, m1_verdict=m1v, battery_verdict=bv, **sh))


# --- (i) exhaustive scan over every truncation d = 1..nsup-1 ---
def full_scan(order_used):
    out = []
    for d in range(nsup - 1, 0, -1):
        b0, b1, b2 = trunc_moments(order_used[:d])
        rec = {'d': d, 'm0_trunc': b0, 'm1_trunc': b1, 'm2_trunc': b2}
        for bf in BIASES:
            g0, g1, g2, m1v, bv = verdicts(b0, b1, b2, bf)
            rec[f'bias{int(round(100*bf))}'] = {'m1_verdict': m1v, 'battery_verdict': bv}
        rec.update({'g0': abs(b0 - m0_op), 'g1': abs(b1 - m1_op), 'g2': abs(b2 - m2_op)})
        out.append(rec)
    return out


def summarize(scan, bias_frac):
    key = f'bias{int(round(100*bias_frac))}'
    d0, d1, d2 = DELTAS[bias_frac]
    m1_miss = [r['d'] for r in scan if r[key]['m1_verdict'] != 'REJECT']
    caught_m2 = [r['d'] for r in scan if r[key]['m1_verdict'] != 'REJECT' and r['g2'] > d2]
    battery_pass = [r['d'] for r in scan if r[key]['battery_verdict'] != 'REJECT']
    # verify.py's auto-scan rule: first d from nsup-1 downward to nsup//3+1 with g1 <= delta1 and g2 > delta2
    vpick = next((r['d'] for r in scan if max(2, nsup // 3) < r['d'] <= nsup - 1
                  and r['g1'] <= d1 and r['g2'] > d2), None)
    return {'bias_frac': bias_frac, 'm1_alone_misses_d': m1_miss, 'm1_misses_caught_by_m2_d': caught_m2,
            'battery_passes_d': battery_pass, 'verify_py_rule_pick_d': vpick}


scan_lex = full_scan(order)
scan_arg = full_scan(order_argsort)
summ_lex = {f'bias{int(round(100*b))}': summarize(scan_lex, b) for b in BIASES}
summ_arg = {f'bias{int(round(100*b))}': summarize(scan_arg, b) for b in BIASES}
print("\n--- (i) every truncation d = 1..%d (deterministic order) ---" % (nsup - 1))
for b in BIASES:
    s = summ_lex[f'bias{int(round(100*b))}']
    print(f" bias {100*b:.0f}%: m1 alone misses at d = {s['m1_alone_misses_d'] or 'none'}; of these m2 catches "
          f"d = {s['m1_misses_caught_by_m2_d'] or 'none'}; whole battery passes d = {s['battery_passes_d'] or 'none'}; "
          f"verify.py rule picks d = {s['verify_py_rule_pick_d']}")
sa = summ_arg['bias2']
print(f" (old argsort order on this machine, 2% bias: m1 alone misses d = {sa['m1_alone_misses_d'] or 'none'}; "
      f"verify.py rule picks d = {sa['verify_py_rule_pick_d']})")

# --- (ii) convention-free shell-boundary cuts ---
bound_rows = []
for d in shell_bounds[:-1]:
    r = next(x for x in scan_lex if x['d'] == d)
    bound_rows.append({'d': d, 'g0': r['g0'], 'g1': r['g1'], 'g2': r['g2'],
                       'm1_verdict_2pct': r['bias2']['m1_verdict'], 'battery_verdict_2pct': r['bias2']['battery_verdict']})
print("\n--- (ii) convention-free shell-boundary cuts (2% bias) ---")
for r in bound_rows:
    print(f" d={r['d']:4d}: |g1|={r['g1']:.4f} -> m1 {r['m1_verdict_2pct']:7s} | |g2|={r['g2']:.3f} -> battery "
          f"{r['battery_verdict_2pct']}")

# --- (iii) tie-break sensitivity at the focus cuts (and verify.py's cut) ---
tie_rows = []
d1_dep, d2_dep = DELTAS[DEPLOYED_BIAS][1], DELTAS[DEPLOYED_BIAS][2]
# focus cuts, verify.py's cut, and every truncation the whole battery passes at the deployed bias
sens_ds = sorted(set([d for d in focus_ds if d < nsup] + [97] +
                     summ_lex[f'bias{int(round(100*DEPLOYED_BIAS))}']['battery_passes_d']), reverse=True)
print(f"\n--- (iii) tie-break sensitivity: {TIE_R} random choices of the kept shell members (seed [{TIE_SEED}, d] per cut), 2% bias ---")
for d in sens_ds:
    sh = shell_at_cut(d)
    if not sh['cut_inside_shell']:
        tie_rows.append(dict(d=d, **sh, note='cut at a shell boundary: no tie to break'))
        print(f" d={d:4d}: shell boundary, verdict convention-free")
        continue
    above = [i for i in order if prob_key[i] > sh['p_cut']]
    members = np.array([i for i in order if prob_key[i] == sh['p_cut']])
    rng = np.random.default_rng([TIE_SEED, d])          # per-cut stream: independent of which cuts are listed
    g1s, g2s, m1rej, batrej = [], [], 0, 0
    for _ in range(TIE_R):
        pick = rng.choice(members, size=sh['n_kept_from_shell'], replace=False)
        b0, b1, b2 = trunc_moments(np.concatenate([above, pick]))
        g0, g1, g2, m1v, bv = verdicts(b0, b1, b2, DEPLOYED_BIAS)
        g1s.append(g1); g2s.append(g2); m1rej += (m1v == 'REJECT'); batrej += (bv == 'REJECT')
    g1s = np.array(g1s); g2s = np.array(g2s)
    rec = dict(d=d, **sh, R=TIE_R, g1_min=float(g1s.min()), g1_median=float(np.median(g1s)), g1_max=float(g1s.max()),
               g2_min=float(g2s.min()), g2_max=float(g2s.max()), frac_m1_reject=m1rej/TIE_R,
               frac_battery_reject=batrej/TIE_R, delta1=d1_dep, delta2=d2_dep)
    tie_rows.append(rec)
    print(f" d={d:4d} ({sh['n_kept_from_shell']} of a {sh['shell_size']}-fold shell): |g1| in [{rec['g1_min']:.3f}, "
          f"{rec['g1_max']:.3f}] (median {rec['g1_median']:.3f}) vs delta1={d1_dep:.3f}; m1 alone rejects on "
          f"{100*rec['frac_m1_reject']:.0f}% of tie-breaks, battery on {100*rec['frac_battery_reject']:.0f}%")

# --- comparison with the committed argsort record (read only) ---
comparison = None
_cpath = os.path.join(DATA, COMMITTED)
if os.path.isfile(_cpath):
    _c = json.load(open(_cpath))['rows']
    comparison = []
    for r in rows:
        c = next((x for x in _c if int(x['d']) == r['d'] and abs(x['bias_frac'] - r['bias_frac']) < 1e-12), None)
        if c is None:
            continue
        comparison.append({'bias_frac': r['bias_frac'], 'd': r['d'], 'g1_committed': c['g1'], 'g1_lexsort': r['g1'],
                           'g2_committed': c['g2'], 'g2_lexsort': r['g2'],
                           'm1_verdict_committed': c['m1_verdict'], 'm1_verdict_lexsort': r['m1_verdict'],
                           'battery_verdict_committed': c['battery_verdict'],
                           'battery_verdict_lexsort': r['battery_verdict'],
                           'm1_verdict_changed': c['m1_verdict'] != r['m1_verdict'],
                           'battery_verdict_changed': c['battery_verdict'] != r['battery_verdict']})

_r98 = next(r for r in rows if r['d'] == 98 and r['bias_frac'] == DEPLOYED_BIAS)
print(f"\nKEY (computed in this run), d=98 at 2% bias: m1 alone -> {_r98['m1_verdict']} "
      f"(|g1| = {_r98['g1']:.4f} vs delta1 = {_r98['delta1']:.4f}); joint battery -> {_r98['battery_verdict']} "
      f"(|g2| = {_r98['g2']:.3f} vs delta2 = {_r98['delta2']:.3f}); the cut keeps {_r98['n_kept_from_shell']} of a "
      f"{_r98['shell_size']}-fold degenerate |psi0|^2 shell.")
if comparison:
    ch = [(c['d'], c['bias_frac'], c['m1_verdict_committed'], c['m1_verdict_lexsort']) for c in comparison
          if c['m1_verdict_changed'] or c['battery_verdict_changed']]
    print(f"verdict changes vs committed {COMMITTED} (d, bias, m1 committed -> lexsort): {ch if ch else 'none'}")
print("A lone moment can miss a truncation that the battery catches; that is a catch, not a guarantee")
print("(within_sector_control.py has redistributions the whole battery passes).")

runtime = time.time() - T_START
out = {'_provenance': {'script': 'src/interval_battery.py', 'date': '2026-09-27', 'plan_item': 'R7/M7',
                       'truncation_order': f'np.lexsort((index, -np.round(|psi0|^2, {ROUND}))): descending probability, '
                                           'ties by ascending sector index (deterministic convention)',
                       'seed': {'tie_break_sensitivity': f'np.random.default_rng([{TIE_SEED}, d]) per cut d',
                                'eigsh': 'spectral_lanczos._sub_gs seeds the start vector with default_rng(0)'},
                       'runtime_s': round(runtime, 2), 'python': sys.version.split()[0], 'numpy': np.__version__,
                       'scipy': scipy.__version__, 'platform': platform.platform(),
                       'supersedes_for_ordering': f'data/{COMMITTED} (argsort, platform-dependent ties; kept unchanged)',
                       'scope': 'shot-budget-informed forecast on the exact state; shot term is a lower bound; '
                                'bias is assumed, not measured'},
       'L': L, 'U': U, 'Ns': Ns, 'z': z, 'seeded_eigsh': True,
       'independent_moments': {'m0': m0_op, 'm1': m1_op, 'm2': m2_op},
       'local_estimator_variances': {'var0': var0, 'var1': var1, 'var2': var2},
       'deltas': {f'bias{int(round(100*b))}': dict(zip(('delta0', 'delta1', 'delta2'), DELTAS[b])) for b in BIASES},
       'order_robustness_check': ORDER_CHECK,
       'shells': {'sizes_desc_p': shell_sizes, 'p_values_desc': [float(v) for v in shell_vals],
                  'convention_free_cuts_d': shell_bounds},
       'focus_ds': focus_ds, 'nsup': nsup, 'rows': rows,
       'full_scan_summary_lexsort': summ_lex,
       'full_scan_summary_argsort_this_machine': summ_arg,
       'full_scan_lexsort': scan_lex,
       'shell_boundary_cuts': bound_rows,
       'tie_break_sensitivity': tie_rows,
       'comparison_with_committed_argsort_run': comparison}
json.dump(out, open(os.path.join(DATA, OUT_NAME), 'w'), indent=1)
print(f"wrote data/{OUT_NAME}  (runtime {runtime:.1f} s)")
