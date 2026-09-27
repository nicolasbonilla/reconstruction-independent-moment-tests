# -*- coding: utf-8 -*-
"""DIRECT MONTE-CARLO DEMONSTRATION of the interval-moment test at the real ibm_fez
shot budget (converts the App-B *forecast* into a *demonstrated-in-simulation* result).
SIM-ONLY, zero QPU: shots are drawn from the exact-state distribution p_x=|<x|0>|^2 --
device bitstring counts were not retained -- so this CALIBRATES the analytic interval at a
classically reproducible scale; it does not certify a device line shape.

Reuses the exact sector/ground-state/local-estimator block of interval_moment.py, then:
  (i)   validates the analytic Gaussian interval (realized std vs analytic; skew/kurtosis),
  (ii)  reports the realized false-positive rate at full reconstruction vs the nominal alpha,
  (iii) measures detection power vs determinant coverage on synthetic truncated reconstructions,
        on the 22-point log grid of the committed figure AND on every truncation d = 1..n_sup-1,
  (iv)  gives a ROC at d=98 (the committed "hard case") and at the lowest-power truncations found.
Honest scope unchanged: refutation-only, necessary-not-sufficient, one probe (m1); the
(m0,m1,m2)+Hankel battery (interval_battery.py) is a partial remedy for a lone-moment miss, not a guarantee.

Truncation order (2026-09-27, plan item R7): np.lexsort((index, -round(|psi0|^2, 12))), i.e. descending
probability with ties broken by ascending sector index. The earlier argsort broke the exactly degenerate
|psi0|^2 shells by floating-point noise (committed d=98: m1_trunc 4.0024, power 0.80; fresh clones: 3.9728,
0.885). The shot replicas (seed 1) do not depend on the order, so the histogram and the scalars are unchanged.
Output: data/2026-09-27_interval_moment_mc.json (the committed 2026-08-24 file is kept as the argsort record);
paper/figs/momentmc_{hist,power,scalars}.dat (grid, same format as before) and momentmc_power_scan.dat (every d).
"""
import time
T_START = time.time()                  # runtime recorded in the provenance includes the imports
import os, sys, json, platform
import numpy as np
import scipy, scipy.sparse as sp
import spectral_lanczos as sl

L, U, Ns = 6, 4.0, 50000          # Heron demonstration scale + REAL ibm_fez shot budget
z = 1.96                          # 95% two-sided
M = 4000                          # Monte-Carlo replicas of the N_s-shot estimator
SEED = 1
ROUND = 12                        # decimals used to identify exactly degenerate |psi0|^2 shells
B_FRAC = 0.02                     # the deployed (assumed) bias budget, used only for the secondary power column
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, '..', 'data'))
OUT_NAME = '2026-09-27_interval_moment_mc.json'
COMMITTED = '2026-08-24_interval_moment_mc.json'   # argsort-ordered record (read only, for the comparison)

# --- exact sector, GS, independent estimator, local-estimator variance (== interval_moment.py) ---
nup = nd = L // 3
Tu, Su, iu = sl.hop(L, nup); Td, Sd, idd = sl.hop(L, nd)
Du, Dd = len(Su), len(Sd)
upocc = sl.occ_matrix(Su, L); dnocc = sl.occ_matrix(Sd, L)
Hs = (sp.kron(Tu, sp.identity(Dd)) + sp.kron(sp.identity(Du), Td)
      + sp.diags(U * (upocc @ dnocc.T).ravel())).tocsr()
Ju = sl.current_string(L, nup); Jd = sl.current_string(L, nd)
Js = (sp.kron(Ju, sp.identity(Dd)) + sp.kron(sp.identity(Du), Jd)).tocsr()
E0, psi0 = sl._sub_gs(Hs)
Jp = Js @ psi0
m1_op = float(np.real(np.vdot(Jp, Hs @ Jp)) - E0 * np.real(np.vdot(Jp, Jp)))
O1psi = np.real(Js @ (Hs @ Jp - E0 * Jp))
psi0r = np.real(psi0)
p = np.abs(psi0) ** 2; p = p / p.sum()
mask = np.abs(psi0) > 1e-12
O1loc = np.zeros(len(psi0)); O1loc[mask] = O1psi[mask] / psi0r[mask]
var_loc = float(np.sum(p * O1loc ** 2) - np.sum(p * O1loc) ** 2)
sigma_shot = np.sqrt(max(var_loc, 0.0) / Ns)
thr = z * sigma_shot              # shot-only 95% threshold (the interval the MC validates)
thr_dep = z * np.sqrt(sigma_shot ** 2 + (B_FRAC * abs(m1_op)) ** 2)   # bias-inclusive deployment threshold

# --- deterministic truncation order (R7): descending |psi0|^2, ties by ascending sector index ---
p_key = np.round(p, ROUND)
order = np.lexsort((np.arange(len(p)), -p_key)); nsup = int((p > 1e-14).sum()); order = order[:nsup]


def shell_at_cut(d):
    pc = p_key[order[d - 1]]
    size = int(np.sum(p_key == pc)); above = int(np.sum(p_key > pc))
    return {'shell_size': size, 'n_kept_from_shell': int(d - above), 'cut_inside_shell': bool(0 < d - above < size)}


def m1_trunc(d):
    idx = np.sort(order[:int(d)])
    Hsub = Hs[idx][:, idx]; Jsub = Js[idx][:, idx]
    e0, g = sl._sub_gs(Hsub); Jg = Jsub @ g
    return float(np.real(np.vdot(Jg, Hsub @ Jg)) - e0 * np.real(np.vdot(Jg, Jg)))


# --- draw M replicas of the N_s-shot independent estimator (once; estimator is d-independent) ---
rng = np.random.default_rng(SEED)
cdf = np.cumsum(p)
mhat = np.empty(M)
for k in range(M):
    idx = np.searchsorted(cdf, rng.random(Ns))
    mhat[k] = O1loc[idx].mean()


def wilson(k, n, zc=1.96):
    ph = k / n; d = 1 + zc**2 / n
    c = (ph + zc**2 / (2*n)) / d
    hw = zc * np.sqrt(ph*(1-ph)/n + zc**2/(4*n**2)) / d
    return max(0.0, c-hw), min(1.0, c+hw)


# --- (i)-(ii) faithfulness + false-positive at full reconstruction (target = m1_op) ---
std_emp = float(mhat.std(ddof=1))
mu = mhat.mean(); sd = std_emp
skew = float(np.mean(((mhat-mu)/sd)**3)); exkurt = float(np.mean(((mhat-mu)/sd)**4) - 3)
fp_rej = int(np.sum(np.abs(mhat - m1_op) > thr))
fp = fp_rej / M; fp_lo, fp_hi = wilson(fp_rej, M)

# --- (iii) power vs determinant coverage: every d (exhaustive) and the committed 22-point log grid ---
scan = {}
for d in range(nsup, 0, -1):
    mt = m1_trunc(d)
    rej = int(np.sum(np.abs(mhat - mt) > thr)); lo, hi = wilson(rej, M)
    rej_dep = int(np.sum(np.abs(mhat - mt) > thr_dep))
    scan[d] = dict({'d': d, 'cov': d / nsup, 'm1_trunc': mt, 'gap': abs(mt - m1_op), 'power': rej / M,
                    'wilson95': [lo, hi], 'power_thr_2pct_bias': rej_dep / M}, **shell_at_cut(d))
ds = sl._logspace_d(nsup, 22, max(20, int(0.05 * nsup)))
power_rows = [(int(d), int(d)/nsup, scan[int(d)]['m1_trunc'], scan[int(d)]['power'],
               scan[int(d)]['wilson95'][0], scan[int(d)]['wilson95'][1]) for d in ds]
trunc = [scan[d] for d in range(nsup - 1, 0, -1)]
low_power = [r['d'] for r in trunc if r['power'] < 0.5]                  # truncations the shot-only test mostly misses
below_full = [r['d'] for r in trunc if r['power'] < 0.95]
grid_min = min((r for r in power_rows if r[0] < nsup), key=lambda r: r[3])
scan_min = min(trunc, key=lambda r: r['power'])


def roc_at(d):
    mt = scan[d]['m1_trunc']
    pts = []
    for zz in np.linspace(0.2, 4.0, 40):
        t = zz * sigma_shot
        pts.append((float(zz), float(np.mean(np.abs(mhat - m1_op) > t)), float(np.mean(np.abs(mhat - mt) > t))))
    fprs = np.array([r[1] for r in pts]); tprs = np.array([r[2] for r in pts])
    o = np.argsort(fprs, kind='stable')
    return {'d': d, 'm1_trunc': mt, 'auc': float(np.trapz(tprs[o], fprs[o])), 'power_z196': scan[d]['power'],
            'points': [{'z': zz, 'fpr': f, 'tpr': t} for zz, f, t in pts]}


roc98 = roc_at(98)
roc_gridmin = roc_at(grid_min[0])
roc_scanmin = roc_at(scan_min['d'])

# --- comparison with the committed argsort-ordered record (read only) ---
comparison = None
_cpath = os.path.join(DATA, COMMITTED)
if os.path.isfile(_cpath):
    _c = {int(r['d']): r for r in json.load(open(_cpath))['power']}
    comparison = [{'d': d, 'm1_trunc_committed': _c[d]['m1_trunc'], 'power_committed': _c[d]['power'],
                   'm1_trunc_lexsort': mt, 'power_lexsort': pw}
                  for d, cov, mt, pw, lo, hi in power_rows if d in _c]

print("=" * 74)
print("INTERVAL-MOMENT MONTE-CARLO (sim-only, ibm_fez budget) -- deterministic lexsort truncation order")
print(f"L={L} U/t={U:.0f}  N_s={Ns}  replicas M={M}  seed={SEED}")
print("=" * 74)
print(f"independent m1_op = {m1_op:.4f}   n_support = {nsup}")
print(f"analytic shot sigma = sqrt(var_loc/N_s) = {sigma_shot:.4f}   (var_loc={var_loc:.1f})")
print(f"realized estimator std = {std_emp:.4f}   ratio realized/analytic = {std_emp/sigma_shot:.3f}")
print(f"sampling distribution: skew = {skew:+.3f}   excess kurtosis = {exkurt:+.3f}  (Gaussian ~ 0,0)")
print(f"realized FALSE-POSITIVE rate at full recon = {100*fp:.2f}%  "
      f"[Wilson95 {100*fp_lo:.1f}-{100*fp_hi:.1f}%]  vs nominal {100*(1-0.95):.0f}%")
print(f"\n power vs coverage on the committed 22-point grid (shot-only threshold {thr:.4f}):")
print(f" {'d':>5} {'cov%':>6} {'m1_trunc':>9} {'power':>7}  [Wilson95]")
for d, cov, mt, pw, lo, hi in power_rows:
    print(f" {d:>5} {100*cov:>6.1f} {mt:>9.3f} {100*pw:>6.1f}%  [{100*lo:.0f}-{100*hi:.0f}]")
print(f"\n lowest-power truncation on the grid: d={grid_min[0]} ({100*grid_min[1]:.1f}%), power {100*grid_min[3]:.1f}%")
print(f" d=98 (committed hard case): m1_trunc={roc98['m1_trunc']:.4f}, power(z=1.96)={100*roc98['power_z196']:.1f}%, "
      f"ROC AUC={roc98['auc']:.3f}")
print(f"\n EXHAUSTIVE scan d=1..{nsup-1}: power < 50% at d = {low_power if low_power else 'none'}")
print(f"                            power < 95% at d = {below_full if below_full else 'none'}")
print(f" lowest power over all truncations: d={scan_min['d']} ({100*scan_min['cov']:.1f}%), "
      f"power {100*scan_min['power']:.1f}% (FP anchor {100*fp:.1f}%)")
if comparison:
    print("\n vs committed " + COMMITTED + " (argsort): changed grid points (d: power committed -> lexsort):")
    for c in comparison:
        if abs(c['m1_trunc_committed'] - c['m1_trunc_lexsort']) > 1e-9:
            print(f"   d={c['d']:4d}: m1_trunc {c['m1_trunc_committed']:.4f} -> {c['m1_trunc_lexsort']:.4f}, "
                  f"power {c['power_committed']:.4f} -> {c['power_lexsort']:.4f}")
print(" -> a lone moment misses whole windows of truncations; see interval_battery.py for what the battery catches")

runtime = time.time() - T_START
out = {'_provenance': {'script': 'src/interval_moment_mc.py', 'date': '2026-09-27', 'plan_item': 'R7/M7',
        'sim_only': True, 'seed': SEED, 'scale': f'L={L} U/t={U} current probe', 'budget': Ns, 'replicas': M,
        'truncation_order': f'np.lexsort((index, -np.round(|psi0|^2, {ROUND}))): descending probability, '
                            'ties by ascending sector index (deterministic convention)',
        'runtime_s': round(runtime, 2), 'python': sys.version.split()[0], 'numpy': np.__version__,
        'scipy': scipy.__version__, 'platform': platform.platform(),
        'supersedes_for_ordering': f'data/{COMMITTED} (argsort, platform-dependent ties; kept unchanged)',
        'note': 'shots drawn from exact-state |<x|0>|^2; device counts not retained; power uses the shot-only '
                'threshold z*sigma_shot (power_thr_2pct_bias uses z*sqrt(sigma^2+(0.02 m1)^2) with the same '
                'unbiased replicas)'},
       'm1_op': m1_op, 'sigma_shot': sigma_shot, 'thr_shot_only': thr, 'thr_2pct_bias': thr_dep,
       'std_emp': std_emp, 'std_ratio': std_emp/sigma_shot, 'skew': skew, 'excess_kurtosis': exkurt,
       'fp_rate': fp, 'fp_wilson95': [fp_lo, fp_hi], 'nsup': nsup,
       'power': [{'d': d, 'cov': cov, 'm1_trunc': mt, 'power': pw, 'wilson95': [lo, hi]}
                 for d, cov, mt, pw, lo, hi in power_rows],
       'power_scan_all_d': [scan[d] for d in range(nsup, 0, -1)],
       'truncations_power_below_50pct': low_power, 'truncations_power_below_95pct': below_full,
       'lowest_power_grid': {'d': grid_min[0], 'cov': grid_min[1], 'power': grid_min[3]},
       'lowest_power_all_d': {'d': scan_min['d'], 'cov': scan_min['cov'], 'power': scan_min['power']},
       'roc_d98': roc98, 'roc_lowest_power_grid': roc_gridmin, 'roc_lowest_power_all_d': roc_scanmin,
       'mhat_hist': np.histogram(mhat, bins=40)[0].tolist(),
       'mhat_edges': np.histogram(mhat, bins=40)[1].tolist(),
       'comparison_with_committed_argsort_run': comparison}
json.dump(out, open(os.path.join(DATA, OUT_NAME), 'w'), indent=2)
print(f"\nwrote data/{OUT_NAME}  (runtime {runtime:.1f} s)")

# --- native-figure .dat mirrors (PGFPlots) into the paper's figs/ dir ---
cnt, edges = np.histogram(mhat, bins=40)
ctr = 0.5 * (edges[:-1] + edges[1:]); bw = edges[1] - edges[0]
gauss = M * bw / (sd * np.sqrt(2*np.pi)) * np.exp(-0.5*((ctr - mu)/sd)**2)
for figdir in (os.path.normpath(os.path.join(HERE, '..', 'paper', 'figs')),):
    if not os.path.isdir(figdir):
        continue
    with open(os.path.join(figdir, 'momentmc_hist.dat'), 'w') as f:
        f.write('x count gauss\n')
        for x, c, g in zip(ctr, cnt, gauss):
            f.write(f'{x:.6g} {c:d} {g:.6g}\n')
    with open(os.path.join(figdir, 'momentmc_power.dat'), 'w') as f:
        f.write('cov power lo hi\n')
        for d, cov, mt, pw, lo, hi in power_rows:
            f.write(f'{100*cov:.4g} {pw:.4g} {lo:.4g} {hi:.4g}\n')
    with open(os.path.join(figdir, 'momentmc_power_scan.dat'), 'w') as f:
        f.write('d cov power lo hi\n')
        for d in range(nsup, 0, -1):
            r = scan[d]
            f.write(f"{d} {100*r['cov']:.4g} {r['power']:.4g} {r['wilson95'][0]:.4g} {r['wilson95'][1]:.4g}\n")
    with open(os.path.join(figdir, 'momentmc_scalars.dat'), 'w') as f:
        f.write(f'm1op {m1_op:.6g}\nsigma {sigma_shot:.6g}\nthr {thr:.6g}\n'
                f'fp {fp:.6g}\nmuhat {mu:.6g}\n')
    print('wrote', figdir + '/momentmc_{hist,power,power_scan,scalars}.dat')
