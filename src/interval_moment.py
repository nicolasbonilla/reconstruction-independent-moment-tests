# -*- coding: utf-8 -*-
"""INTERVAL-MOMENT CLOSURE (Route A make-or-break, $0): upgrade the exact-moment idealization to the ROBUST
interval-moment test the paper repeatedly promises but never shows. The falsification threshold is set by the
REAL IBM Heron shot budget (N_s=50000, ibm_fez, arXiv:2608.16436), NOT an arbitrary 5%. Produces a genuine
CORROBORATE (full reconstruction) vs REJECT (truncated) VERDICT from the interval test.
Honest scope: current-probe instance at the Heron demonstration scale/budget (L=6, U=4); the device is NOT
load-bearing at this classically-reproducible scale (full-sector recovery) -- this closes the interval machinery,
it does not manufacture device-essentiality (that is Route B / gate design).

Truncation order (2026-09-27, plan item R7): the top-d determinants are ranked by |psi0|^2 with a DETERMINISTIC
tie-break, np.lexsort((index, -round(|psi0|^2, 12))): descending probability, then ascending sector index. The
ground-state probabilities come in exactly degenerate shells (6/12/24-fold, by lattice symmetry), and the earlier
np.argsort(prob)[::-1] broke those ties by floating-point noise, so which members of the shell at the cut were kept
depended on the platform (the committed 2026-08-22 d=98 row, 3.9646, is one such draw; fresh clones gave 3.9728).
The lexsort choice inside a shell is reproducible but still a CONVENTION: the shell at each cut is recorded in the
output, and interval_battery.py measures how much the verdict depends on it. All verdicts below are computed in
this run, never asserted. Output: data/2026-09-27_interval_moment_closure.json (the committed 2026-08-22 file is
left untouched as the record of the argsort run)."""
import time
T_START = time.time()                  # runtime recorded in the provenance includes the imports
import numpy as np, json, os, sys, platform
import scipy
import spectral_lanczos as sl

L, U, Ns = 6, 4.0, 50000            # Heron demonstration scale + REAL ibm_fez shot budget
z = 1.96                            # 95% interval
b_frac = 0.02                       # ASSUMED residual-bias budget as a fraction of m1 (not measured). The ibm_fez
                                    # runs used measurement twirling, Pauli gate twirling and dynamical decoupling;
                                    # no readout-error mitigation was applied ("TREX" in the draft is a misnomer).
ROUND = 12                          # decimals used to identify exactly degenerate |psi0|^2 shells
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, '..', 'data'))
OUT_NAME = '2026-09-27_interval_moment_closure.json'
COMMITTED = '2026-08-22_interval_moment_closure.json'    # argsort-ordered record (read only, for the comparison)

# --- build the sector, ground state, independent estimator m1_op, and its LOCAL-estimator shot variance ---
import scipy.sparse as sp
nup = nd = L // 3
Tu, Su, iu = sl.hop(L, nup); Td, Sd, idd = sl.hop(L, nd)
Du, Dd = len(Su), len(Sd)
upocc = sl.occ_matrix(Su, L); dnocc = sl.occ_matrix(Sd, L)
diagU = U * (upocc @ dnocc.T).ravel()
Hs = (sp.kron(Tu, sp.identity(Dd)) + sp.kron(sp.identity(Du), Td) + sp.diags(diagU)).tocsr()
Ju = sl.current_string(L, nup); Jd = sl.current_string(L, nd)
Js = (sp.kron(Ju, sp.identity(Dd)) + sp.kron(sp.identity(Du), Jd)).tocsr()
E0, psi0 = sl._sub_gs(Hs)
Jp = Js @ psi0
m1_op = float(np.real(np.vdot(Jp, Hs @ Jp)) - E0 * np.real(np.vdot(Jp, Jp)))   # independent estimator (exact mean)
m0 = float(np.real(np.vdot(Jp, Jp)))                                            # <J^2>

# LOCAL-estimator shot variance: O1|0> = J(H-E0)J|0>; O1_loc(x) = (O1 psi0)_x / psi0_x ; sample dist p_x=|c_x|^2
O1psi = Js @ (Hs @ Jp - E0 * Jp)
p = np.abs(psi0) ** 2
mask = np.abs(psi0) > 1e-12
O1loc = np.zeros(len(psi0), dtype=np.result_type(O1psi, psi0))
O1loc[mask] = O1psi[mask] / psi0[mask]
O1loc = np.real(O1loc)
mean_loc = float(np.sum(p * O1loc))                     # == m1_op (check)
var_loc = float(np.sum(p * O1loc ** 2) - mean_loc ** 2) # local-estimator variance
sigma_shot = np.sqrt(max(var_loc, 0.0) / Ns)
bias = b_frac * abs(m1_op)
delta1 = z * np.sqrt(sigma_shot ** 2 + bias ** 2)       # the robust interval half-width (REAL 50k-shot budget)

# --- deterministic truncation order: descending |psi0|^2, ties by ascending sector index (R7) ---
prob = np.abs(psi0) ** 2
prob_key = np.round(prob, ROUND)
order = np.lexsort((np.arange(len(prob)), -prob_key))
nsup = int((prob > 1e-14).sum()); order = order[:nsup]


def shell_at_cut(d):
    """the exactly degenerate |psi0|^2 shell containing the d-th kept determinant (1-based)."""
    pc = prob_key[order[d - 1]]
    size = int(np.sum(prob_key == pc)); above = int(np.sum(prob_key > pc))
    return {'p_cut': float(pc), 'shell_size': size, 'n_above_shell': above, 'n_kept_from_shell': int(d - above),
            'cut_inside_shell': bool(0 < d - above < size)}


ds = sl._logspace_d(nsup, 22, max(20, int(0.05 * nsup)))
rows = []
for d in ds:
    d = int(d); idx = np.sort(order[:d])
    Hsub = Hs[idx][:, idx]; Jsub = Js[idx][:, idx]
    e0, g = sl._sub_gs(Hsub); Jg = Jsub @ g
    m1_trunc = float(np.real(np.vdot(Jg, Hsub @ Jg)) - e0 * np.real(np.vdot(Jg, Jg)))
    gap = abs(m1_trunc - m1_op)
    verdict = 'REJECT' if gap > delta1 else 'corroborate'
    rows.append((d, d / nsup, m1_trunc, gap, verdict))
# largest d that still REJECTs
d_reject_max = max([r[0] for r in rows if r[4] == 'REJECT'], default=None)
misses = [r[0] for r in rows if r[0] < nsup and r[4] != 'REJECT']     # truncated, yet m1 alone corroborates

# --- comparison with the committed argsort-ordered record (read only) ---
comparison = None
_cpath = os.path.join(DATA, COMMITTED)
if os.path.isfile(_cpath):
    _c = {int(r['d']): r for r in json.load(open(_cpath))['sweep']}
    comparison = []
    for d, cov, m1t, gap, v in rows:
        c = _c.get(d)
        if c is None:
            continue
        comparison.append({'d': d, 'm1_trunc_committed': c['m1_trunc'], 'gap_committed': c['gap'],
                           'verdict_committed': c['verdict'], 'm1_trunc_lexsort': m1t, 'gap_lexsort': gap,
                           'verdict_lexsort': v, 'm1_trunc_changed': bool(abs(c['m1_trunc'] - m1t) > 1e-9),
                           'verdict_changed': bool(c['verdict'] != v)})

print("=== INTERVAL-MOMENT CLOSURE (real 50k-shot budget, ibm_fez) -- deterministic lexsort truncation order ===")
print(f"L={L} U={U}  sector={Du*Dd}  n_support={nsup}  E0={E0:.4f}")
print(f"independent m1_op = {m1_op:.5f}   m0=<J^2> = {m0:.4f}   (local-est check {mean_loc:.5f})")
print(f"local-estimator var = {var_loc:.4f}  -> shot sigma(50k) = {sigma_shot:.5f}")
print(f"mitigation-bias budget b = {bias:.5f} ({100*b_frac:.0f}% of m1)")
print(f"ROBUST INTERVAL half-width delta_1 = {delta1:.5f}  (= {100*delta1/abs(m1_op):.2f}% of m1)")
print(f"\n FULL reconstruction (d={nsup}): gap={abs(rows[0][3]):.2e} < delta -> {rows[0][4].upper()}")
print(" d       cov%    m1_trunc   |gap|     verdict   cut shell (size, kept)")
for d, cov, m1t, gap, v in rows:
    sh = shell_at_cut(d)
    tag = f"inside {sh['shell_size']}-fold shell, {sh['n_kept_from_shell']} kept" if sh['cut_inside_shell'] else 'shell boundary'
    print(f" {d:6d}  {100*cov:5.1f}   {m1t:8.4f}  {gap:8.4f}   {v:11s} {tag}")
print(f"\nHONEST FRAMING (per adversarial verify wf_e83094d3-0bd): this is a SHOT-BUDGET-INFORMED FORECAST of the")
print(f"robust interval-moment screen at the ibm_fez scale/budget, computed as a MODEL on the EXACT state (device")
print(f"counts unsaved). shot term = ground-state local-estimator variance (a LOWER BOUND on true per-shot")
print(f"uncertainty); bias term = an ASSUMED 2% residual (swept, not measured); the interval is {100*delta1/abs(m1_op):.1f}% of m1.")
print(f"A WIDER real interval only makes the refutation-only screen MORE conservative (fewer, safer rejects). The gap")
print(f"here is vs the EXACT oracle m1_op (a calibration at accessible scale), NOT the deployable same-sample falsifier.")
_r98 = [r for r in rows if r[0] == 98]
if _r98:   # computed in this run
    sh = shell_at_cut(98)
    print(f"d=98 in this run: |gap| = {_r98[0][3]:.4f} vs delta_1 = {delta1:.4f} -> {_r98[0][4]} by m1 alone "
          f"(cut inside a {sh['shell_size']}-fold |psi0|^2 shell: {sh['n_above_shell']} above it, "
          f"{sh['n_kept_from_shell']} of its {sh['shell_size']} kept).")
print(f"truncated grid points NOT rejected by m1 alone (2% bias): {misses if misses else 'none'}")
if comparison:
    ch = [c['d'] for c in comparison if c['verdict_changed']]
    mv = [c['d'] for c in comparison if c['m1_trunc_changed']]
    print(f"vs committed {COMMITTED}: m1_trunc changed at d={mv}; verdict changed at d={ch if ch else 'none'}")
print(f"A lone moment can miss a truncation (interval_battery.py scans every d); the joint (m0,m1,m2)+Hankel battery")
print(f"is the remedy, not a guarantee. Device NOT load-bearing at this full-sector scale.")

runtime = time.time() - T_START
out = {'_provenance': {'script': 'src/interval_moment.py', 'date': '2026-09-27', 'plan_item': 'R7/M7',
       'scale': 'L=6 U=4 current-probe, Heron demonstration scale',
       'shot_budget': Ns, 'backend': 'ibm_fez (arXiv:2608.16436)',
       'truncation_order': f'np.lexsort((index, -np.round(|psi0|^2, {ROUND}))): descending probability, '
                           'ties by ascending sector index (deterministic convention)',
       'seed': 'none drawn; spectral_lanczos._sub_gs seeds the eigsh start vector with default_rng(0)',
       'runtime_s': round(runtime, 2), 'python': sys.version.split()[0], 'numpy': np.__version__,
       'scipy': scipy.__version__, 'platform': platform.platform(),
       'supersedes_for_ordering': f'data/{COMMITTED} (argsort, platform-dependent ties; kept unchanged)',
       'honest_scope': 'closes the robust interval-moment gap; device NOT load-bearing at this full-sector-recovery scale'},
       'm1_op': m1_op, 'm0': m0, 'var_loc': var_loc, 'sigma_shot': sigma_shot, 'bias': bias, 'b_frac': b_frac,
       'z': z, 'delta1': delta1, 'delta1_pct': 100*delta1/abs(m1_op), 'nsup': nsup, 'd_reject_max': d_reject_max,
       'truncated_grid_points_not_rejected_by_m1': misses,
       'sweep': [dict({'d': d, 'cov': cov, 'm1_trunc': m1t, 'gap': gap, 'verdict': v}, **shell_at_cut(d))
                 for d, cov, m1t, gap, v in rows],
       'comparison_with_committed_argsort_run': comparison}
json.dump(out, open(os.path.join(DATA, OUT_NAME), 'w'), indent=2)
print(f"\nwrote data/{OUT_NAME}  (runtime {runtime:.1f} s)")
