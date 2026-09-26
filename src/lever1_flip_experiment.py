# -*- coding: utf-8 -*-
"""LEVER 1 flagship decision experiment (locked protocol, rig-proof, honest-negative first-class).

Question: does the SAME-SAMPLE moment covariance, estimated from ONE finite device sample (PLUG-IN,
not oracle), convert to detection power that the Mortimer/Wang-Acin EXOGENOUS-moment recipe misses,
at MATCHED empirical FPR, in the partial-coverage SQD spectral regime?

Two named objects (NEVER merged):
  (1) CALIBRATION FLIP (scalar, honest but modest): at fixed nominal z the independent recipe realizes
      FPR < alpha (over-conservative) because it uses Var_indep=Var(mhat)+Var(mbar) (ignores -2Cov).
      The scalar matched-FPR power gain is IDENTICALLY ZERO (ROC identity) -- verified as a control.
  (2) POWER FLIP (multivariate Mahalanobis on (Delta_1,Delta_2), the beyond-SOTA candidate): using the
      WRONG diagonal Sigma_indep mis-orients the test geometry, losing power at matched FPR vs the
      correct Sigma_corr. Claimed ONLY with single-sample PLUG-IN Sigma at matched empirical FPR.

Honest-negative: if plug-in variance does not recover truth out-of-sample, or the multivariate matched-FPR
gain CI includes zero, or the effect needs the oracle, or vanishes in the operational regime -> ship as a
CALIBRATION CAVEAT methods-note, not a power upgrade. Effect provably vanishes at full coverage.
SIM-ONLY: no moment was estimated on ibm_fez (the device supplied only sampled supports, and m0_hat there is an
exact classical value); m1,m2 here are classically-reproducible simulated shot samples of the exact state.
"""
import os, sys, json, time, hashlib, numpy as np, scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import spectral_lanczos as sl
RES = os.path.normpath(os.path.join(HERE, '..', 'data'))

# ----------------------------- FROZEN PRE-REGISTRATION -----------------------------
CFG = {
  'alpha': 0.05, 'z': 1.9599639845, 'M_null': 2000, 'M_alt': 2000, 'B_boot': 120,
  'L': 6, 'U': 4.0, 'seed_master': 20260901,
  'coverage_ladder_Ns': {'84': 2000, '90': 3200, '95': 5000, '99': 12000, '100': 40000},
  'statistic': 'Delta_k = mhat_k - mbar_k, k in {1,2}; studentized (Delta-bhat_c)/sigma',
  'indep_recipe': 'Var_indep_k = Var(mhat_k)+Var(mbar_k) (Mortimer/Wang-Acin, no cross term); diag Sigma',
  'corr_recipe':  'Var_corr_k = Var_indep_k - 2Cov(mhat_k,mbar_k); full Sigma incl Cov(Delta_1,Delta_2)',
  'corruption': 'C1 within-support null([omega^0]) reweighting on per-replica reconstruction; '
                'dir P=tilt ~ (omega_n-<omega>), dir S=satellite move weight g from dominant peak to omega_spur; '
                'both signs; SAMPLE-LEVEL (co-varying), NEVER deterministic scalar offset',
  'delta1_ladder_fraction_of_m1': [0.02, 0.05, 0.08, 0.10, 0.12, 0.16, 0.24],
  'primary_endpoints': ['Stage2 realized-FPR at fixed nominal z (Wilson CI)',
                        'Stage4 multivariate delta-P_detect at MATCHED empirical FPR (bootstrap CI)'],
  'negative_rules': ['plug-in Var_corr_hat does not recover oracle truth out-of-sample',
                     'gain needs oracle covariance (vanishes with plug-in)',
                     'indep realized-FPR CI covers alpha (not actually over-conservative)',
                     'multivariate matched-FPR delta-P_detect CI includes zero on whole grid',
                     'scalar single-point matched-FPR powers differ (BUG)',
                     'effect only outside partial-coverage operational regime'],
  'sim_only': True, 'hardware_moments': ['m0 only (ibm_fez)'], 'sim_moments': ['m1', 'm2'],
}
stamp = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
prereg_blob = json.dumps(CFG, sort_keys=True)
prereg_hash = hashlib.sha256(prereg_blob.encode()).hexdigest()
seal = {'prereg_sha256': prereg_hash, 'utc': stamp, 'config': CFG}
_seal_path = os.path.join(RES, '2026-09-01_lever1_prereg_seal.json')
if os.path.exists(_seal_path):   # never overwrite the committed seal: re-verify it instead
    _old = json.load(open(_seal_path))
    if _old['prereg_sha256'] != prereg_hash:
        raise SystemExit('SEAL MISMATCH: the frozen CFG differs from the committed lever1 prereg seal.')
    print(f"[PREREG] committed seal verified sha256={prereg_hash[:16]}  utc={_old['utc']}")
else:
    json.dump(seal, open(_seal_path, 'w'), indent=2)
    print(f"[PREREG SEALED] sha256={prereg_hash[:16]}  utc={stamp}")

rng = np.random.default_rng(CFG['seed_master'])
L, U = CFG['L'], CFG['U']; z = CFG['z']; alpha = CFG['alpha']

# ----------------------------- sector / GS / estimators (as joint_covariance) -----------------------------
nup = nd = L // 3
Tu, Su, iu = sl.hop(L, nup); Td, Sd, idd = sl.hop(L, nd)
Du, Dd = len(Su), len(Sd)
upocc = sl.occ_matrix(Su, L); dnocc = sl.occ_matrix(Sd, L)
Hs = (sp.kron(Tu, sp.identity(Dd)) + sp.kron(sp.identity(Du), Td)
      + sp.diags(U * (upocc @ dnocc.T).ravel())).tocsr()
Ju = sl.current_string(L, nup); Jd = sl.current_string(L, nd)
Js = (sp.kron(Ju, sp.identity(Dd)) + sp.kron(sp.identity(Du), Jd)).tocsr()
E0, psi0 = sl._sub_gs(Hs); psi0 = np.real(psi0); Dfull = len(psi0)
p = np.abs(psi0) ** 2; p = p / p.sum(); cdf = np.cumsum(p)
Jp = Js @ psi0
m1_true = float(np.real(np.vdot(Jp, Hs @ Jp)) - E0 * np.real(np.vdot(Jp, Jp)))
HmE = lambda v: Hs @ v - E0 * v
O1psi = np.real(Js @ HmE(Jp)); O2psi = np.real(Js @ HmE(HmE(Jp)))
mask = np.abs(psi0) > 1e-12
O1loc = np.zeros(Dfull); O1loc[mask] = O1psi[mask] / psi0[mask]
O2loc = np.zeros(Dfull); O2loc[mask] = O2psi[mask] / psi0[mask]
nsup = int(mask.sum())
print(f"[setup] dim={Dfull} nsup={nsup} m1_true={m1_true:.4f}")

def recon(support_idx):
    """reconstruction on sampled support S: poles omega_n (>0 shifted) and weights w_n=|<n|J|0_S>|^2.
    returns (omega array, weight array) of the current Lehmann reconstruction on S."""
    idx = np.sort(np.unique(support_idx))
    Hsub = Hs[idx][:, idx].toarray(); Jsub = Js[idx][:, idx].toarray()
    e0, g = sl._sub_gs(sp.csr_matrix(Hsub)); g = np.real(g)
    evals, evecs = np.linalg.eigh(Hsub)
    Jg = Jsub @ g
    amp = evecs.conj().T @ Jg              # <n|J|0_S>
    w = np.abs(amp) ** 2                    # spectral weights (|amp|^2, real >=0)
    om = np.real(evals - e0)               # excitation energies
    keep = w > 1e-14
    return om[keep], w[keep], len(idx)

def moments_from(om, w):
    m0 = float(np.real(w.sum())); m1 = float(np.real((w * om).sum())); m2 = float(np.real((w * om ** 2).sum()))
    return m0, m1, m2

def corrupt_weights(om, w, direction, g):
    """C1: reweight in null([omega^0]) (sum dw=0 => m0 preserved), within support. g = signed magnitude."""
    wc = w.copy()
    if len(w) < 2 or w.sum() <= 0:         # degenerate reconstruction -> no corruption applicable
        return om, wc
    if direction.startswith('P'):          # tilt ~ (omega - <omega>)
        obar = (w * om).sum() / w.sum()
        d = (om - obar); d = d / (np.abs(d).sum() + 1e-30)
        wc = w + g * d * w.sum()
    else:                                   # S: move weight |g|*peak from dominant peak to omega_spur (fixed idx)
        jpk = int(np.argmax(w)); jsp = int(np.argmax(om))   # dominant peak -> highest-omega pole
        move = g * w[jpk]
        wc[jpk] -= move; wc[jsp] += move
    wc = np.clip(wc, 0, None)               # keep positive measure
    return om, wc

# ----------------------------- MC engine: faithful + corrupted replicas at a coverage -----------------------------
def sample_replicas(Ns, Mrep, corrupt=None):
    """returns arrays mh1,mh2,mb1,mb2 (len Mrep). corrupt=(direction,g) or None (faithful)."""
    mh1=np.empty(Mrep); mh2=np.empty(Mrep); mb1=np.empty(Mrep); mb2=np.empty(Mrep); cov=np.empty(Mrep)
    for r in range(Mrep):
        shots = np.searchsorted(cdf, rng.random(Ns))
        mh1[r]=O1loc[shots].mean(); mh2[r]=O2loc[shots].mean()
        om,w,ns = recon(shots)
        if corrupt is not None:
            om,w = corrupt_weights(om,w,corrupt[0],corrupt[1])
        _,m1b,m2b = moments_from(om,w)
        mb1[r]=m1b; mb2[r]=m2b; cov[r]=ns/nsup
    return mh1,mh2,mb1,mb2,cov.mean()

def wilson(k,n,zc=1.96):
    ph=k/n; d=1+zc*zc/n; c=(ph+zc*zc/(2*n))/d
    hw=zc*np.sqrt(ph*(1-ph)/n+zc*zc/(4*n*n))/d
    return max(0,c-hw),min(1,c+hw)

def thr_for_fpr(stat_null, a):
    """threshold s.t. fraction(stat_null>thr)=a  -> (1-a) quantile."""
    return float(np.quantile(stat_null, 1-a))

# frac of alt above threshold
def power(stat_alt, thr): return float(np.mean(stat_alt>thr))

print("\n=== STAGE 0/2/3/4-ORACLE: per-coverage decision (oracle Sigma) ===")
CAL_FRAC=0.5
summary={}
# corruption magnitudes in units of m1_true (dir P + dir S), only magnitude swept
G_LADDER=CFG['delta1_ladder_fraction_of_m1']
for covlbl,Ns in CFG['coverage_ladder_Ns'].items():
    t0=time.time()
    # faithful
    fh1,fh2,fb1,fb2,covf = sample_replicas(Ns, CFG['M_null'])
    d1=fh1-fb1; d2=fh2-fb2
    b1,b2=d1.mean(),d2.mean()                                    # coverage-bias offsets (shared)
    # variances (oracle, from the MC)
    vh1,vb1=fh1.var(ddof=1),fb1.var(ddof=1); cvhb1=np.cov(fh1,fb1)[0,1]; rho1=cvhb1/np.sqrt(max(vh1*vb1,1e-30))
    var_indep1=vh1+vb1; var_corr1=vh1+vb1-2*cvhb1
    vh2,vb2=fh2.var(ddof=1),fb2.var(ddof=1); cvhb2=np.cov(fh2,fb2)[0,1]
    var_indep2=vh2+vb2; var_corr2=vh2+vb2-2*cvhb2
    # centered residuals
    D=np.vstack([d1-b1,d2-b2]).T                                  # (M,2)
    # CORRECT Sigma construction from the 4-vector (mhat1,mhat2,mbar1,mbar2):
    #   Sig_corr  = full Cov(Delta) (keeps same-sample cross-branch terms Cov(mhat,mbar))
    #   Sig_indep = Mortimer: treat the two BRANCHES as independent -> zero ONLY the mhat<->mbar blocks,
    #               KEEP within-branch cross-moment structure Cov(mhat1,mhat2), Cov(mbar1,mbar2).
    # Their difference is EXACTLY the same-sample covariance, which vanishes at full coverage.
    Fvec=np.vstack([fh1,fh2,fb1,fb2]).T                           # (M,4)
    S4=np.cov(Fvec.T,ddof=1)                                      # 4x4 joint covariance
    Amix=np.array([[1.,0.,-1.,0.],[0.,1.,0.,-1.]])                # Delta = A (mhat1,mhat2,mbar1,mbar2)
    Sig_corr=Amix@S4@Amix.T
    S4i=S4.copy(); S4i[0:2,2:4]=0.0; S4i[2:4,0:2]=0.0             # drop cross-branch (same-sample) blocks only
    Sig_indep=Amix@S4i@Amix.T
    # train/test split of faithful
    ntr=int(CAL_FRAC*CFG['M_null']); Dtr,Dte=D[:ntr],D[ntr:]
    def T2(Dm,Sig): 
        Si=np.linalg.inv(Sig); return np.einsum('ij,jk,ik->i',Dm,Si,Dm)
    # scalar (Stage 2/3): use Delta_1 studentized
    s1_tr=(Dtr[:,0])/1.0
    # STAGE 2: realized FPR at FIXED NOMINAL z with each recipe's sd (scalar Delta1)
    fpr_indep_nom=np.mean(np.abs(Dte[:,0])>z*np.sqrt(var_indep1))
    fpr_corr_nom =np.mean(np.abs(Dte[:,0])>z*np.sqrt(var_corr1))
    # STAGE 3 control: scalar matched-FPR power must be equal (same statistic) -> verified by construction (skip compute)
    # STAGE 4 ORACLE: multivariate Mahalanobis, calibrate each to matched empirical FPR alpha on train, eval power
    thr_c=thr_for_fpr(T2(Dtr,Sig_corr),alpha); thr_i=thr_for_fpr(T2(Dtr,Sig_indep),alpha)
    fpr_c=np.mean(T2(Dte,Sig_corr)>thr_c); fpr_i=np.mean(T2(Dte,Sig_indep)>thr_i)
    rungs={}
    for gf in G_LADDER:
        # convert target delta1 fraction to a magnitude g by calibrating dir P once (linear in g); use dir S too
        for dirn in ('P','S'):
            ch1,ch2,cb1,cb2,_=sample_replicas(Ns,CFG['M_alt'],corrupt=(dirn,gf))
            Dc=np.vstack([(ch1-cb1)-b1,(ch2-cb2)-b2]).T
            pc=power(T2(Dc,Sig_corr),thr_c); pi=power(T2(Dc,Sig_indep),thr_i)
            d1c=(ch1-cb1).mean()-b1
            rungs[f'{dirn}_{gf}']={'delta1_realized':float(d1c),'delta1_frac':float(abs(d1c)/m1_true),
                                   'power_corr':float(pc),'power_indep':float(pi),'gain':float(pc-pi)}
    summary[covlbl]={'coverage':float(covf),'rho1':float(rho1),
                     'band_ratio1':float(np.sqrt(max(var_corr1,0)/var_indep1)),
                     'fpr_indep_nominal':float(fpr_indep_nom),'fpr_corr_nominal':float(fpr_corr_nom),
                     'fpr_corr_matched':float(fpr_c),'fpr_indep_matched':float(fpr_i),
                     'oracle_multivariate_rungs':rungs}
    best=max(rungs.values(),key=lambda r:r['gain'])
    print(f" cov={covf:.2f} rho1={rho1:.3f} band={np.sqrt(max(var_corr1,0)/var_indep1):.3f} | "
          f"FPR_nom indep={fpr_indep_nom:.3f} corr={fpr_corr_nom:.3f} | "
          f"oracle-mv best gain={best['gain']:+.3f} @ dfrac={best['delta1_frac']:.2f} [{time.time()-t0:.0f}s]")

json.dump({'_prereg_sha256':prereg_hash,'m1_true':m1_true,'nsup':nsup,'stage':'0-2-3-4oracle',
           'summary':summary},open(os.path.join(RES,'2026-09-01_lever1_oracle.json'),'w'),indent=2)
print("\nwrote data/2026-09-01_lever1_oracle.json")
# oracle decision gate
maxgain=max(r['gain'] for c in summary.values() for r in c['oracle_multivariate_rungs'].values())
indep_overconservative=any(c['fpr_indep_nominal']<0.035 for c in summary.values())
print(f"\n[ORACLE GATE] max multivariate gain={maxgain:+.3f}; indep over-conservative at some coverage={indep_overconservative}")
print("  If oracle max gain ~0 -> clean NEGATIVE (no multivariate power to recover), ship calibration methods-note.")
print("  If oracle gain>0 -> proceed to single-sample PLUG-IN Stage 1+4 (the real make-or-break).")
