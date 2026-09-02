# -*- coding: utf-8 -*-
"""LEVER 1 (crux experiment): the SAME-SAMPLE COVARIANCE of the falsifier Delta_k = m_hat_k - m_bar_k.

The beyond-SOTA gap (SOTA panel, verified 2024-2026): the independent operator-moment estimator
m_hat_k and the reconstruction's read-off moment m_bar_k are computed from the SAME finite device
sample set, so Cov(m_hat_k, m_bar_k) != 0 and Var(Delta_k) = Var(m_hat_k)+Var(m_bar_k)-2Cov.
Wang-Acin/Mortimer (finite-stat SDP) and the super-resolution CRB literature both assume the moment
estimate is an INDEPENDENT exogenous input -> they use Var_indep = Var(m_hat)+Var(m_bar).
This script MEASURES the covariance honestly by a joint Monte-Carlo: per replica, draw N_s shots ONCE
and compute BOTH branches from those SAME shots. It reports rho_k, Var_corr/Var_indep, and the
decision-band ratio sqrt(Var_corr)/sqrt(Var_indep). The whole beyond-SOTA case is CONDITIONAL on
this covariance being load-bearing (rho_k appreciable, band ratio != 1); if rho_k ~ 0 the correlation
is cosmetic and we say so plainly (honest fallback). SIM-ONLY, $0. L=6 U/t=8 doped current measure.
"""
import os, sys, json, numpy as np, scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import spectral_lanczos as sl

L, U = 6, 4.0
SEED = 20260901
M = 3000                    # MC replicas
rng = np.random.default_rng(SEED)

# --- exact sector / GS / independent local estimators (m1,m2), same block as interval_moment_mc ---
nup = nd = L // 3
Tu, Su, iu = sl.hop(L, nup); Td, Sd, idd = sl.hop(L, nd)
Du, Dd = len(Su), len(Sd)
upocc = sl.occ_matrix(Su, L); dnocc = sl.occ_matrix(Sd, L)
Hs = (sp.kron(Tu, sp.identity(Dd)) + sp.kron(sp.identity(Du), Td)
      + sp.diags(U * (upocc @ dnocc.T).ravel())).tocsr()
Ju = sl.current_string(L, nup); Jd = sl.current_string(L, nd)
Js = (sp.kron(Ju, sp.identity(Dd)) + sp.kron(sp.identity(Du), Jd)).tocsr()
E0, psi0 = sl._sub_gs(Hs); psi0 = np.real(psi0)
Dfull = len(psi0)
p = np.abs(psi0) ** 2; p = p / p.sum(); cdf = np.cumsum(p)
Jp = Js @ psi0
# exact true moments (full operator)
m1_true = float(np.real(np.vdot(Jp, Hs @ Jp)) - E0 * np.real(np.vdot(Jp, Jp)))
m0_true = float(np.real(np.vdot(Jp, Jp)))
# local estimators O_k^loc(x): E_p[O_k^loc] = m_k_true  (m1 and m2)
HmE = lambda v: Hs @ v - E0 * v
O1psi = np.real(Js @ HmE(Jp))                       # J (H-E0) J psi
O2psi = np.real(Js @ HmE(HmE(Jp)))                  # J (H-E0)^2 J psi
mask = np.abs(psi0) > 1e-12
O1loc = np.zeros(Dfull); O1loc[mask] = O1psi[mask] / psi0[mask]
O2loc = np.zeros(Dfull); O2loc[mask] = O2psi[mask] / psi0[mask]
nsup = int(mask.sum())
print(f"L={L} U/t={U:.0f}  sector dim={Dfull}  support nsup={nsup}  m0_true={m0_true:.4f} m1_true={m1_true:.4f}")

# --- per-replica reconstruction moment m_bar_k on the SAMPLED support S_r (SQD: rediagonalize on S) ---
def recon_moments(support_idx):
    idx = np.sort(np.unique(support_idx))
    Hsub = Hs[idx][:, idx]; Jsub = Js[idx][:, idx]
    e0, g = sl._sub_gs(Hsub); g = np.real(g)
    Jg = Jsub @ g
    m0b = float(np.vdot(Jg, Jg).real)
    m1b = float(np.vdot(Jg, (Hsub @ Jg - e0 * Jg)).real)
    m2b = float(np.vdot(Jg, (Hsub @ (Hsub @ Jg - e0 * Jg) - e0 * (Hsub @ Jg - e0 * Jg))).real)
    return m0b, m1b, m2b, len(idx)

def run_at(Ns):
    mh1 = np.empty(M); mh2 = np.empty(M); mb1 = np.empty(M); mb2 = np.empty(M); cov = np.empty(M)
    for r in range(M):
        shots = np.searchsorted(cdf, rng.random(Ns))     # ONE shot draw
        mh1[r] = O1loc[shots].mean()                       # independent estimator branch (frequencies)
        mh2[r] = O2loc[shots].mean()
        m0b, m1b, m2b, ns = recon_moments(shots)           # reconstruction branch (support), SAME shots
        mb1[r] = m1b; mb2[r] = m2b; cov[r] = ns / nsup
    out = {}
    for name, mh, mb in [('m1', mh1, mb1), ('m2', mh2, mb2)]:
        d = mh - mb
        vh, vb = mh.var(ddof=1), mb.var(ddof=1)
        cvhb = np.cov(mh, mb, ddof=1)[0, 1]
        rho = cvhb / np.sqrt(max(vh * vb, 1e-30))
        var_indep = vh + vb                                # Mortimer-style (independent)
        var_corr = vh + vb - 2 * cvhb                       # correlation-aware (this work)
        out[name] = {'mean_mhat': float(mh.mean()), 'mean_mbar': float(mb.mean()),
                     'mean_delta': float(d.mean()), 'var_mhat': float(vh), 'var_mbar': float(vb),
                     'cov': float(cvhb), 'rho': float(rho),
                     'var_indep': float(var_indep), 'var_corr': float(var_corr),
                     'band_ratio_corr_over_indep': float(np.sqrt(max(var_corr,0)/var_indep)),
                     'sd_indep': float(np.sqrt(var_indep)), 'sd_corr': float(np.sqrt(max(var_corr,0)))}
    out['mean_coverage'] = float(cov.mean())
    return out

print(f"\n{'Ns':>7} {'cov':>6} | {'moment':>6} {'rho':>7} {'sd_indep':>9} {'sd_corr':>9} {'band_ratio':>10} {'mean_dbar':>10}")
results = {}
for Ns in [2000, 5000, 12000, 30000, 60000]:
    r = run_at(Ns); results[str(Ns)] = r
    for mm in ('m1', 'm2'):
        x = r[mm]
        print(f"{Ns:>7} {r['mean_coverage']:>6.2f} | {mm:>6} {x['rho']:>7.3f} {x['sd_indep']:>9.4f} "
              f"{x['sd_corr']:>9.4f} {x['band_ratio_corr_over_indep']:>10.3f} {x['mean_delta']:>10.4f}")

# --- honest verdict on load-bearingness ---
verdict = {}
for Ns, r in results.items():
    for mm in ('m1', 'm2'):
        x = r[mm]
        # load-bearing if the correlation-aware band differs materially from independent (>~10%) with definite sign
        lb = abs(1 - x['band_ratio_corr_over_indep']) > 0.10
        verdict[f'{mm}@{Ns}'] = {'rho': x['rho'], 'band_ratio': x['band_ratio_corr_over_indep'], 'load_bearing_candidate': bool(lb)}
any_lb = any(v['load_bearing_candidate'] for v in verdict.values())
print(f"\nLOAD-BEARING CANDIDATE (|1-band_ratio|>0.10 somewhere)? {any_lb}")
print("  (this is necessary but NOT sufficient for a verdict flip; a flip needs a realistic corruption")
print("   landing between the two bands. If rho~0 everywhere, the covariance is cosmetic -> honest fallback.)")

RES = os.path.normpath(os.path.join(HERE, '..', '06_results'))
json.dump({'_provenance': {'script': 'joint_covariance_composition.py', 'sim_only': True, 'seed': SEED,
           'model': f'L={L} U/t={U} doped current Lehmann measure', 'M': M,
           'gap': 'same-sample Cov(m_hat_k,m_bar_k); Var(Delta)=Var(mhat)+Var(mbar)-2Cov'},
           'm0_true': m0_true, 'm1_true': m1_true, 'nsup': nsup, 'results': results, 'verdict': verdict,
           'any_load_bearing_candidate': any_lb},
          open(os.path.join(RES, '2026-09-01_joint_covariance.json'), 'w'), indent=2)
print("\nwrote 06_results/2026-09-01_joint_covariance.json")
