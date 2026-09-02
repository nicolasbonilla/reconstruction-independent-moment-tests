# -*- coding: utf-8 -*-
"""OFF-DIAGONAL MOMENT DISCRIMINANT UNDER ibm_fez-ANCHORED NOISE  (Option A, SIM-ONLY, ZERO QPU).

Converts the paper's *admitted* off-diagonal gap into a falsifiable, quantitative FORECAST:
the density-matrix-EXACT upward bias floor B(eps) of the reconstruction-independent DENSITY-probe
first moment m1 = <drho (H-E0) drho> (PSD operator, biased toward Tr(M1)/dim) versus per-CZ 2q
gate error eps, bracketed [twirled-depolarizing (optimistic), coherent (pessimistic)], with the
crossing eps*(Delta) where a PRE-REGISTERED reconstruction corruption can still be refuted ABOVE
the floor + shot band. Honest expected outcome = NEGATIVE-WITH-A-NUMBER (off-diagonal unusable at
ibm_fez noise; eps* sits ~1 decade below current 2q error). Locked by the 7-agent design panel.

RED LINES (design panel): noise lives on CZ-bearing carriers, NEVER on the 1q measurement rotation;
floor is density-matrix EXACT (not an assumed global-depolarizing scalar); floor delivered as a
BRACKET; corruption size referenced to m1_EXACT (never the biased noisy estimate); m0 is a CONTROL
not a finding; L=6/12q is exactly diagonalizable => calibration/forecast, NOT beyond-classical;
Aer cannot model leakage => reported floor is a LOWER bound on badness, eps* an UPPER bound on
tolerable error.
"""
import os, json, time, hashlib
import numpy as np
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error
import bond_moment_estimator as bme

# ----------------------------- FROZEN MANIFEST (pre-registered) -----------------------------
MANIFEST = {
    'script': 'offdiag_noise_forecast.py', 'sim_only': True, 'qpu_used': False,
    'master_seed': 20260830, 'robustness_seeds': [1, 2, 3, 4, 5],
    'model': 'Hubbard L=6 chain OBC, 12 qubits, U/t=4, density probe rho_q q=pi',
    'estimand': 'm1 = <drho (H-E0) drho>  (PSD, biased toward Tr(M1)/dim under depolarizing)',
    'eps_grid': [0.0, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 2e-3, 2.5e-3, 3e-3, 5e-3, 1e-2, 3e-2, 1e-1],
    'ibm_fez_2q_anchor': 2.5e-3, 'ibm_fez_readout': [0.008, 0.018], 'ibm_fez_1q': 2.5e-4,
    'G_headline': 292, 'G_surface': [20, 60, 100, 200, 292],
    'corruption_coverages': [1.0, 0.9, 0.75, 0.5, 0.25, 0.1],
    'ns_shots': 50000, 'shot_reps': 32, 'z': 1.96,
    'fair_test_margin': 5.0,            # corruption must clear 5x shot-CI noiselessly
    'crossing_def': 'largest eps with Delta - |B(eps)| > z*sigma_eff (upward bias-aligned sign)',
    'pins': {'qiskit': '1.4.6', 'qiskit_aer': '0.17.2'},
    'honest_expectation': 'NEGATIVE-WITH-A-NUMBER: eps* well below ibm_fez 2.5e-3',
}
Z = MANIFEST['z']; NS = MANIFEST['ns_shots']; GH = MANIFEST['G_headline']

# ----------------------------- STEP 2: operators + RE-VERIFY anchors -----------------------------
Mops, E0, psi, m_exact = bme.build_operators(nmom=3)
M0op, M1op, M2op = Mops
NQ = bme.NQ; DIM = 2 ** NQ
gs = np.asarray(psi.data, dtype=complex)
m0e, m1e, m2e = m_exact
M0d = M0op.to_matrix(sparse=True).toarray()
M1d = M1op.to_matrix(sparse=True).toarray()

def tr_over_dim(M):
    for lab, co in M.to_list():
        if set(lab) == {'I'}:
            return float(np.real(co))
    return 0.0
trM0, trM1 = tr_over_dim(M0op), tr_over_dim(M1op)
Bmax = trM1 - m1e                         # full-depolarizing bias asymptote for m1
assert abs(m1e - 10.6085) < 1e-3 and abs(trM1 - 34.2662) < 1e-3, "anchor drift -- ABORT"
print(f"[anchors] E0={E0:.4f}  m1_exact={m1e:.4f}  Tr(M1)/dim={trM1:.4f}  ratio={trM1/m1e:.3f}  Bmax={Bmax:.3f}")
print(f"[anchors] m0_exact={m0e:.4f}  Tr(M0)/dim={trM0:.4f}  ratio={trM0/m0e:.3f}  (diagonal CONTROL)")

# ----------------------------- STEP 3: density-probe corruption ladder -----------------------------
Hmat = bme.build_hubbard().to_matrix(sparse=True).tocsr()
rho_mat = bme.build_rho().to_matrix(sparse=True).tocsr()
import scipy.sparse.linalg as sla
p_gs = np.abs(gs) ** 2
order = np.argsort(p_gs)[::-1]
nsup = int((p_gs > 1e-14).sum())
order = order[:nsup]

def m1_trunc_density(frac):
    """density-probe first moment on a determinant-truncated (top-frac support) reconstruction."""
    d = max(2, int(round(frac * nsup)))
    idx = np.sort(order[:d])
    Hs = Hmat[idx][:, idx]; rs = rho_mat[idx][:, idx]
    e0, V = sla.eigsh(Hs.tocsc(), k=1, which='SA')
    g = V[:, 0]
    rev = float(np.real(g.conj() @ (rs @ g)))
    dg = rs @ g - rev * g                       # drho |g>
    HmE = Hs @ dg - float(e0[0]) * dg
    return float(np.real(dg.conj() @ HmE))

ladder = []
for fr in MANIFEST['corruption_coverages']:
    mt = m1_trunc_density(fr)
    ladder.append({'coverage': fr, 'm1_trunc': mt, 'Delta': abs(mt - m1e), 'signed': mt - m1e})
print("[ladder] density-probe truncation corruptions (Delta = |m1_trunc - m1_exact|):")
for L in ladder:
    print(f"   cov={L['coverage']:.2f}  m1_trunc={L['m1_trunc']:.4f}  Delta={L['Delta']:.4f} (signed {L['signed']:+.4f})")

# ----------------------------- shot sigma (noiseless, for the fair-test + band) -----------------------------
est, ng = bme.estimate_moments(psi, [M0op, M1op], ns=NS, seed=MANIFEST['master_seed'], reps=MANIFEST['shot_reps'])
sigma_shot = float(est[:, 1].std(ddof=1))
print(f"[shot] noiseless m1 estimator: {ng} bond-basis groups, sigma_shot(Ns={NS})={sigma_shot:.4f}")
# fair-test gate: every corruption must clear 5x the shot CI noiselessly
for L in ladder:
    L['fair'] = bool(L['Delta'] > MANIFEST['fair_test_margin'] * Z * sigma_shot)

# ----------------------------- STEP 4-6: exact noisy floor B(eps) via density_matrix -----------------------------
SIM_DM = AerSimulator(method='density_matrix'); SIM_DM.set_options(fusion_enable=False)
SIM_SV = AerSimulator(method='statevector'); SIM_SV.set_options(fusion_enable=False)
PAIRS = [(i, i + 1) for i in range(NQ - 1)]

def _carriers(qc, npairs, coherent_theta=0.0):
    """npairs identity carriers cz;cz (= I). If coherent_theta>0: cz;RZZ(th);cz;RZZ(th) coherent edge."""
    for k in range(npairs):
        a, b = PAIRS[k % len(PAIRS)]
        qc.cz(a, b)
        if coherent_theta: qc.rzz(coherent_theta, a, b)
        qc.cz(a, b)
        if coherent_theta: qc.rzz(coherent_theta, a, b)

def moment_on_rho(rho, Md):
    return float(np.real(np.tensordot(rho, Md.T, axes=([0, 1], [0, 1]))))

def moments_depol(eps, ncz, Mds, seed=0):
    """DENSITY-MATRIX exact: depolarizing 2q noise on cz carriers. Returns [<Md>] for each Md in Mds."""
    qc = QuantumCircuit(NQ); qc.set_statevector(gs)
    _carriers(qc, ncz // 2, coherent_theta=0.0)
    qc.save_density_matrix()
    nm = NoiseModel()
    if eps > 0:
        nm.add_all_qubit_quantum_error(depolarizing_error(eps, 2), 'cz')
    rho = np.asarray(SIM_DM.run(qc, noise_model=nm, seed_simulator=seed).result().data()['density_matrix'])
    return [moment_on_rho(rho, Md) for Md in Mds]

def moment_coherent(eps, ncz, Md):
    """STATEVECTOR (fast): coherent ZZ over-rotation theta=sqrt(eps) per cz, NO decoherence -> pure state."""
    th = np.sqrt(eps)
    qc = QuantumCircuit(NQ); qc.set_statevector(gs)
    _carriers(qc, ncz // 2, coherent_theta=th)
    qc.save_statevector()
    sv = np.asarray(SIM_SV.run(qc).result().data()['statevector'])
    return float(np.real(sv.conj() @ (Md @ sv)))

# ----------------------------- primary sweep at G_headline -----------------------------
t0 = time.time()
sweep = []
for eps in MANIFEST['eps_grid']:
    m1_d, m0_d = moments_depol(eps, GH, [M1d, M0d], seed=MANIFEST['master_seed'])
    m1_c = moment_coherent(eps, GH, M1d)
    peff = 1.0 - (1.0 - eps) ** GH
    B_gd = peff * Bmax                                  # global-depolarizing FOLKLORE overlay
    sweep.append({'eps': eps, 'B_depol': m1_d - m1e, 'B_coherent': m1_c - m1e,
                  'B_folklore': B_gd, 'm0_bias': m0_d - m0e, 'peff': peff,
                  'm1_depol': m1_d, 'm1_coherent': m1_c, 'm0_depol': m0_d})
    print(f"  eps={eps:.1e}  B_depol={m1_d-m1e:+.3f}  B_coher={m1_c-m1e:+.3f}  "
          f"B_folklore={B_gd:+.3f}  m0_bias={m0_d-m0e:+.4f}  [{time.time()-t0:.0f}s]")

# ----------------------------- STEP 11: crossing eps*(Delta) per corruption -----------------------------
def crossing(Bcurve):
    """largest eps with Delta - |B(eps)| > z*sigma_shot, per ladder rung (bias-dominated)."""
    epsg = [s['eps'] for s in sweep]; Bg = [abs(s[Bcurve]) for s in sweep]
    out = []
    for L in ladder:
        D = L['Delta']; thr = Z * sigma_shot
        usable = [e for e, b in zip(epsg, Bg) if (D - b) > thr]
        estar = max(usable) if usable else 0.0
        # interpolate the bias=Delta-band crossing for a smooth number
        cross = None
        for i in range(len(epsg) - 1):
            f0 = D - Bg[i] - thr; f1 = D - Bg[i + 1] - thr
            if f0 > 0 >= f1:
                lo, hi = epsg[i], epsg[i + 1]
                if lo > 0:
                    cross = float(np.exp(np.log(lo) + (np.log(hi) - np.log(lo)) * f0 / (f0 - f1)))
                else:
                    cross = float(hi * f0 / (f0 - f1))
                break
        out.append({'coverage': L['coverage'], 'Delta': D, 'eps_star_grid': estar,
                    'eps_star_interp': cross, 'usable_at_ibmfez': bool(D - abs(
                        next(s[Bcurve] for s in sweep if abs(s['eps'] - MANIFEST['ibm_fez_2q_anchor']) < 1e-9)) > thr)})
    return out

cross_depol = crossing('B_depol'); cross_coher = crossing('B_coherent')
print("\n[crossing] eps*(Delta)  (largest 2q error where the corruption still clears floor+band):")
for cd, cc in zip(cross_depol, cross_coher):
    ed = cd['eps_star_interp']; ec = cc['eps_star_interp']
    print(f"   cov={cd['coverage']:.2f} Delta={cd['Delta']:.3f}  eps*_depol="
          f"{ed if ed else 0:.2e}  eps*_coherent={ec if ec else 0:.2e}  "
          f"usable@ibm_fez(2.5e-3)={cd['usable_at_ibmfez']}")

# ----------------------------- non-folklore deviation test -----------------------------
dev = [{'eps': s['eps'], 'B_depol': s['B_depol'], 'B_folklore': s['B_folklore'],
        'rel_dev': (s['B_depol'] - s['B_folklore']) / max(abs(s['B_folklore']), 1e-9)} for s in sweep if s['eps'] > 0]

# ----------------------------- write results + manifest hash -----------------------------
RES = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '06_results'))
os.makedirs(RES, exist_ok=True)
out = {'_manifest': MANIFEST,
       '_manifest_hash': hashlib.sha256(json.dumps(MANIFEST, sort_keys=True).encode()).hexdigest()[:16],
       'anchors': {'E0': E0, 'm0_exact': m0e, 'm1_exact': m1e, 'm2_exact': m2e,
                   'TrM0_dim': trM0, 'TrM1_dim': trM1, 'ratio_m1': trM1 / m1e,
                   'Bmax': Bmax, 'sigma_shot': sigma_shot, 'n_groups_m0m1': ng},
       'ladder': ladder, 'sweep_G%d' % GH: sweep,
       'crossing_depol': cross_depol, 'crossing_coherent': cross_coher,
       'folklore_deviation': dev,
       'runtime_s': time.time() - t0}
stamp = time.strftime('%Y-%m-%d')
fn = os.path.join(RES, f'{stamp}_offdiag_noise_forecast.json')
json.dump(out, open(fn, 'w'), indent=2)
print(f"\n[write] {fn}   manifest_hash={out['_manifest_hash']}   runtime={out['runtime_s']:.0f}s")
print("[expectation check] honest negative-with-a-number:",
      "off-diagonal UNUSABLE at ibm_fez" if not any(c['usable_at_ibmfez'] for c in cross_depol)
      else "SOME corruptions usable at ibm_fez -- inspect")
