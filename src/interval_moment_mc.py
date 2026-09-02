# -*- coding: utf-8 -*-
"""DIRECT MONTE-CARLO DEMONSTRATION of the interval-moment test at the real ibm_fez
shot budget (converts the App-B *forecast* into a *demonstrated-in-simulation* result).
SIM-ONLY, zero QPU: shots are drawn from the exact-state distribution p_x=|<x|0>|^2 --
device bitstring counts were not retained -- so this CALIBRATES the analytic interval at a
classically reproducible scale; it does not certify a device line shape.

Reuses the exact sector/ground-state/local-estimator block of interval_moment.py, then:
  (i)  validates the analytic Gaussian interval (realized std vs analytic; skew/kurtosis),
  (ii) reports the realized false-positive rate at full reconstruction vs the nominal alpha,
  (iii) measures detection power vs determinant coverage on synthetic truncated reconstructions,
  (iv) gives a ROC at the hard d=98 case (accidental first-moment match).
Honest scope unchanged: refutation-only, necessary-not-sufficient, one probe (m1); the
(m0,m1,m2)+Hankel battery closes the d=98 blind spot exhibited here (interval_battery.py).
"""
import os, json
import numpy as np
import scipy.sparse as sp
import spectral_lanczos as sl

L, U, Ns = 6, 4.0, 50000          # Heron demonstration scale + REAL ibm_fez shot budget
z = 1.96                          # 95% two-sided
M = 4000                          # Monte-Carlo replicas of the N_s-shot estimator
SEED = 1

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

# --- reconstruction moment m1_trunc(d) on determinant-truncated subspaces ---
order = np.argsort(p)[::-1]; nsup = int((p > 1e-14).sum()); order = order[:nsup]
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

# --- (iii) power vs determinant coverage ---
ds = sl._logspace_d(nsup, 22, max(20, int(0.05 * nsup)))
power_rows = []
for d in ds:
    mt = m1_trunc(d)
    rej = int(np.sum(np.abs(mhat - mt) > thr)); pw = rej / M
    lo, hi = wilson(rej, M)
    power_rows.append((int(d), int(d)/nsup, mt, pw, lo, hi))

# --- (iv) ROC at the hard d=98 (accidental first-moment match) ---
d_hard = 98; mt98 = m1_trunc(d_hard)
zs = np.linspace(0.2, 4.0, 40)
roc = []
for zz in zs:
    t = zz * sigma_shot
    tpr = float(np.mean(np.abs(mhat - mt98) > t))   # reject a WRONG (truncated) reconstruction
    fpr = float(np.mean(np.abs(mhat - m1_op) > t))   # reject the TRUE reconstruction
    roc.append((float(zz), fpr, tpr))
# trapezoid AUC over fpr (sorted)
fprs = np.array([r[1] for r in roc]); tprs = np.array([r[2] for r in roc])
o = np.argsort(fprs); auc98 = float(np.trapz(tprs[o], fprs[o]))
pw98 = float(np.mean(np.abs(mhat - mt98) > thr))

print("=" * 74)
print("INTERVAL-MOMENT MONTE-CARLO (sim-only, ibm_fez budget)")
print(f"L={L} U/t={U:.0f}  N_s={Ns}  replicas M={M}  seed={SEED}")
print("=" * 74)
print(f"independent m1_op = {m1_op:.4f}   n_support = {nsup}")
print(f"analytic shot sigma = sqrt(var_loc/N_s) = {sigma_shot:.4f}   (var_loc={var_loc:.1f})")
print(f"realized estimator std = {std_emp:.4f}   ratio realized/analytic = {std_emp/sigma_shot:.3f}")
print(f"sampling distribution: skew = {skew:+.3f}   excess kurtosis = {exkurt:+.3f}  (Gaussian ~ 0,0)")
print(f"realized FALSE-POSITIVE rate at full recon = {100*fp:.1f}%  "
      f"[Wilson95 {100*fp_lo:.1f}-{100*fp_hi:.1f}%]  vs nominal {100*(1-0.95):.0f}%")
print(f"\n power vs coverage (reject a truncated reconstruction):")
print(f" {'d':>5} {'cov%':>6} {'m1_trunc':>9} {'power':>7}  [Wilson95]")
for d, cov, mt, pw, lo, hi in power_rows[::3]:
    print(f" {d:>5} {100*cov:>6.1f} {mt:>9.3f} {100*pw:>6.1f}%  [{100*lo:.0f}-{100*hi:.0f}]")
print(f"\n hard case d={d_hard} ({100*d_hard/nsup:.0f}% coverage, m1_trunc={mt98:.3f}): "
      f"power(z=1.96)={100*pw98:.0f}%, ROC AUC={auc98:.3f}")
print(f" -> closed by the (m0,m1,m2)+Hankel battery (interval_battery.py)")

out = {'_provenance': {'script': 'interval_moment_mc.py', 'sim_only': True, 'seed': SEED,
        'scale': f'L={L} U/t={U} current probe', 'budget': Ns, 'replicas': M,
        'note': 'shots drawn from exact-state |<x|0>|^2; device counts not retained'},
       'm1_op': m1_op, 'sigma_shot': sigma_shot, 'std_emp': std_emp,
       'std_ratio': std_emp/sigma_shot, 'skew': skew, 'excess_kurtosis': exkurt,
       'fp_rate': fp, 'fp_wilson95': [fp_lo, fp_hi], 'nsup': nsup,
       'power': [{'d': d, 'cov': cov, 'm1_trunc': mt, 'power': pw, 'wilson95': [lo, hi]}
                 for d, cov, mt, pw, lo, hi in power_rows],
       'roc_d98': {'d': d_hard, 'm1_trunc': mt98, 'auc': auc98, 'power_z196': pw98,
                   'points': [{'z': zz, 'fpr': f, 'tpr': t} for zz, f, t in roc]},
       'mhat_hist': np.histogram(mhat, bins=40)[0].tolist(),
       'mhat_edges': np.histogram(mhat, bins=40)[1].tolist()}
os.makedirs('../06_results', exist_ok=True)
json.dump(out, open('../06_results/2026-08-24_interval_moment_mc.json', 'w'), indent=2)
print("\nwrote 06_results/2026-08-24_interval_moment_mc.json")

# --- native-figure .dat mirrors (PGFPlots) into the paper's figs/ dirs ---
cnt, edges = np.histogram(mhat, bins=40)
ctr = 0.5 * (edges[:-1] + edges[1:]); bw = edges[1] - edges[0]
gauss = M * bw / (sd * np.sqrt(2*np.pi)) * np.exp(-0.5*((ctr - mu)/sd)**2)
for figdir in ('../paper/arxiv-submission/figs', '../paper/figs'):
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
    with open(os.path.join(figdir, 'momentmc_scalars.dat'), 'w') as f:
        f.write(f'm1op {m1_op:.6g}\nsigma {sigma_shot:.6g}\nthr {thr:.6g}\n'
                f'fp {fp:.6g}\nmuhat {mu:.6g}\n')
    print('wrote', figdir + '/momentmc_{hist,power,scalars}.dat')
