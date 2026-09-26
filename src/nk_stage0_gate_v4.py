# -*- coding: utf-8 -*-
"""GATE v4 — the amendment-validation gate (chair ruling wf_c07c723b-b96, sealed verbatim).

Purpose: validate the ONE permitted amendment of manifest v1 (floors 150/100us -> 100/70us,
backend set B={fez,marrakesh,kingston}, NEW per-qubit readout cap <=0.030, per-edge CZ cap)
by verdict rows at EXACTLY the amended values. If R1^R2^R3 pass -> manifest v2 reseal;
if R1 fails solely via readout -> mechanical ladder 0.030->0.025->0.018; if R1@0.018 fails
or R2 or R3 fail -> amendment DENIED in toto, manifest v1 stands (gate-negative branch).
ONE-AMENDMENT DISCIPLINE: this is the first and last floor amendment on this manifest line.

Verdict rows (ALL must pass; >=10 replica seeds each; criteria i/ii/iii + lambda=0.50
falsifier 10/10 per row; iv-D via the sealed v3-R2 band rule mean+1.645*SE < 0.0266667):
  R1 dominating worst-admissible: T1=100us T2=70us  eps=3.5e-3  ro=0.030 uniform
  R2 mid-bracket (v3 continuity): T1=100us T2=70us  eps=2.85e-3 ro=0.018
  R3 best-case (kingston regime): T1=200us T2=130us eps=2.0e-3  ro=0.010
Non-verdict: S1 stress T1=50 T2=35 eps=3.5e-3 ro=0.030 (recorded, never quotable).
Diagnostics (mirror survival, halfB, CZ-ratio) recorded on seed0 of each row.

SIM-ONLY, $0. Writes data/2026-09-04_nk_v4_row_<tag>.json (never overwriting an existing row
file) and, only when every verdict row runs in one process, data/2026-09-04_nk_stage0_gate_v4.json.
The rows were in fact run as separate processes; the sealed verdict is nk_stage0_gate_v4_combine.py
-> data/2026-09-05_nk_v4_verdict.json (AMENDMENT-DENIED, gate-negative).
"""
import os, sys, json, time, zlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from nk_stage0_gate import (L, FILL, SHOTS, EPS1Q, ed_references, ft_matrix, build_circuits,
                            rotation_matrix_of_circuit, orbital_rotation_circuit,
                            nk_from_statevector)
from nk_stage0_gate_v2 import cz_pauli_table, build_twirled_instance
from nk_stage0_gate_v3 import build_isa
from qiskit.circuit import Delay
from qiskit.quantum_info import Statevector
from qiskit.transpiler import PassManager
from qiskit_aer import AerSimulator
from qiskit_aer.noise import (NoiseModel, depolarizing_error, ReadoutError,
                              thermal_relaxation_error, RelaxationNoisePass)

SEED = 20260904
N_TWIRL, N_REPLICAS = 32, 10
DUR_CZ, DUR_1Q = 84e-9, 32e-9
THR_IV = 0.20*(2.0/3.0)/5.0                 # 0.0266667 (sealed, unchanged)
LAM_ROWS = ['0.20', '0.35', '0.50']         # 0.50 primary; 0.20 recorded, never quotable
ROWS = [
    ('R1', 100e-6,  70e-6, 3.5e-3,  0.030, True),
    ('R2', 100e-6,  70e-6, 2.85e-3, 0.018, True),
    ('R3', 200e-6, 130e-6, 2.0e-3,  0.010, True),
    ('S1',  50e-6,  35e-6, 3.5e-3,  0.030, False),
    # SPECULATIVE PRE-COMPUTATION of the sealed readout-ladder rung 2 (chair ruling: the
    # ladder 0.030 -> 0.025 -> 0.018 is mechanical and pre-committed). Computing it in
    # parallel changes NO decision rule — only the order of execution. Tag matches exactly
    # what the sequential ladder would produce, so the result is identical either way.
    ('R1ro0.025', 100e-6, 70e-6, 3.5e-3, 0.025, True),
    # Ladder rung 3 (final sealed rung), also pre-computed speculatively: rung 2 is
    # statistically decided (4/4 replicas above threshold, sd 0.0012), so the sealed
    # mechanical ladder descends here. No decision rule changes.
    ('R1ro0.018', 100e-6, 70e-6, 3.5e-3, 0.018, True),
]
CIRCS = ('rung0', 'mirrorFT', 'calibFT', 'rungB', 'calibB2')

def make_noise(eps, ro, t1, t2):
    g2 = thermal_relaxation_error(t1, t2, DUR_CZ)
    err2 = depolarizing_error(eps, 2).compose(g2.tensor(g2))
    err1 = depolarizing_error(EPS1Q, 1).compose(thermal_relaxation_error(t1, t2, DUR_1Q))
    nm = NoiseModel()
    nm.add_all_qubit_quantum_error(err2, ['cz'])
    nm.add_all_qubit_quantum_error(err1, ['sx', 'x'])
    nm.add_all_qubit_readout_error(ReadoutError([[1-ro, ro], [ro, 1-ro]]))
    return nm

def pipeline(counts, ro):
    """Post-select -> symmetric readout inversion at the ROW's uniform ro -> spin-avg n_k."""
    tot = sum(counts.values()); kept = 0; occ = np.zeros(12)
    for bstr, c in counts.items():
        b = int(bstr, 2)
        if bin(b & 0x3F).count('1') == 2 and bin((b >> 6) & 0x3F).count('1') == 2:
            kept += c
            for q in range(12):
                if (b >> q) & 1: occ[q] += c
    if kept == 0: return np.full(L, np.nan), 0.0
    marg = np.clip((occ/kept - ro)/(1.0 - 2.0*ro), 0.0, 1.0)
    return np.array([(marg[k] + marg[6+k])/2.0 for k in range(L)]), kept/tot

def run_row(tag, t1, t2, eps, ro, isa, refs, table, relax_qubits, log):
    theta_zz = float(np.sqrt(eps))
    relax_pm = PassManager(RelaxationNoisePass(t1s=[t1]*relax_qubits, t2s=[t2]*relax_qubits,
                                               dt=4e-9, op_types=[Delay]))
    maxdev_B = float(np.max(np.abs(refs['rungB_nk'] - FILL)))
    maxdev_F = float(np.max(np.abs(refs['rung0_nk'] - FILL)))
    reps = []
    for stag in range(N_REPLICAS):
        rng = np.random.default_rng(SEED + 104729*stag + zlib.crc32(tag.encode()) % 9973)
        sim = AerSimulator(method='density_matrix', noise_model=make_noise(eps, ro, t1, t2),
                           seed_simulator=SEED + 1009*stag)
        out = {}
        for nm_ in CIRCS:
            insts = [relax_pm.run(build_twirled_instance(isa[nm_], theta_zz, rng, table))
                     for _ in range(N_TWIRL)]
            res = sim.run(insts, shots=SHOTS//N_TWIRL).result()
            pool = {}
            for i in range(N_TWIRL):
                for k_, v_ in res.get_counts(i).items(): pool[k_] = pool.get(k_, 0) + v_
            nk, keep = pipeline(pool, ro)
            out[nm_] = (nk, keep)
        def p_fit(nm_, refkey):
            ref = refs[refkey]; d = ref - FILL
            return float(np.clip(np.sum((ref - out[nm_][0])*d)/np.sum(d*d), 0, 1))
        p2, pfB = p_fit('calibB2', 'calibB2_nk'), p_fit('rungB', 'rungB_nk')
        pFT = p_fit('calibFT', 'rung0_nk')
        ntB = (out['rungB'][0] - p2*FILL)/(1-p2)
        B_B = float(np.max(np.abs(ntB - refs['rungB_nk'])))
        ntF = (out['mirrorFT'][0] - pFT*FILL)/(1-pFT)
        B_F = float(np.max(np.abs(ntF - refs['rung0_nk'])))
        keepB = out['rungB'][1]
        sig_shot = 0.5/np.sqrt(SHOTS*max(keepB, 1e-3))/np.sqrt(2)
        sig_eff = float(np.sqrt((sig_shot/(1-p2))**2 + (maxdev_B/(1-p2)**2*0.02)**2))
        lam = {ls: {'crit_i': bool(B_B < float(ls)*maxdev_B/5.0),
                    'crit_ii': bool(float(ls)*maxdev_B - B_B > 1.96*sig_eff)} for ls in LAM_ROWS}
        reps.append({'seed': stag, 'B_rungB': B_B, 'B_FT_ivD': B_F, 'p_calibB2': p2,
                     'p_truth_rungB': pfB, 'crit_iii': bool(abs(p2-pfB) < 0.02),
                     'keep_rungB': keepB, 'keep_ok': bool(keepB >= 0.25), 'lambdas': lam})
        log(f"[{tag} seed{stag}] B={B_B:.4f} ivD={B_F:.4f} p2={p2:.3f} truth={pfB:.3f} keep={keepB:.2f}")
    # row verdict
    BF = np.array([r['B_FT_ivD'] for r in reps])
    ucb = float(BF.mean() + 1.645*BF.std(ddof=1)/np.sqrt(len(BF)))
    iv_pass = bool(ucb < THR_IV)
    all_iii = all(r['crit_iii'] for r in reps)
    all_keep = all(r['keep_ok'] for r in reps)
    lam50 = all(r['lambdas']['0.50']['crit_i'] and r['lambdas']['0.50']['crit_ii'] for r in reps)
    row_pass = bool(iv_pass and all_iii and all_keep and lam50)
    fail_channels = [c for c, ok in (('iv_readout_or_model', iv_pass), ('iii', all_iii),
                                     ('keep', all_keep), ('lambda50', lam50)) if not ok]
    return {'replicas': reps, 'ivD_mean': float(BF.mean()), 'ivD_ucb': ucb,
            'ivD_thr': THR_IV, 'row_pass': row_pass, 'fail_channels': fail_channels}

def main():
    t0 = time.time()
    lines = []
    def log(s): print(s, flush=True); lines.append(s)
    refs, _ = ed_references()
    W = ft_matrix().conj().T
    assert np.max(np.abs(rotation_matrix_of_circuit(orbital_rotation_circuit, W) - W)) < 1e-9
    circs = build_circuits(W)
    val_map = {'rung0': 'rung0_nk', 'mirrorFT': 'rung0_nk', 'rungB': 'rungB_nk',
               'calibB2': 'calibB2_nk', 'calibFT': 'rung0_nk'}
    isa, czc, meta, dt = build_isa(circs, refs, val_map)
    log('[validate] logical+ISA < 1e-8 (build_isa asserts)')
    table = cz_pauli_table()
    nq = max(v.num_qubits for v in isa.values())

    row_filter = [r for r in os.environ.get('ROWS', '').split(',') if r]
    results = {}
    for (tag, t1, t2, eps, ro, verdict) in ROWS:
        if row_filter and tag not in row_filter: continue
        log(f"=== ROW {tag}: T1={t1*1e6:.0f}us T2={t2*1e6:.0f}us eps={eps:g} ro={ro} "
            f"({'VERDICT' if verdict else 'non-verdict'}) ===")
        results[tag] = run_row(tag, t1, t2, eps, ro, isa, refs, table, nq, log)
        results[tag]['params'] = {'T1_us': t1*1e6, 'T2_us': t2*1e6, 'eps': eps, 'ro': ro,
                                  'verdict_row': verdict}
        log(f"=== {tag}: {'PASS' if results[tag]['row_pass'] else 'FAIL'} "
            f"(ivD ucb={results[tag]['ivD_ucb']:.5f} vs {THR_IV:.5f}; fails={results[tag]['fail_channels']}) ===")

    res_dir = os.path.normpath(os.path.join(HERE, '..', 'data'))
    for tag_, r_ in results.items():
        _rp = os.path.join(res_dir, f'2026-09-04_nk_v4_row_{tag_}.json')
        if os.path.exists(_rp):   # committed row files are inputs of the sealed verdict
            log(f'NOT overwriting existing {os.path.basename(_rp)} (record of the sealed gate)')
            continue
        with open(_rp, 'w') as f:
            json.dump(r_, f, indent=1)
    # ---- amendment verdict + mechanical readout ladder (only in the R1 process; the
    #      cross-row combination is done by the combiner once all row files exist) ----
    if 'R1' not in results:
        log('row subset complete (no amendment verdict in this process)')
        return
    ladder_note = None
    r1 = results['R1']
    if not r1['row_pass'] and r1['fail_channels'] == ['iv_readout_or_model']:
        for ro2 in (0.025, 0.018):
            log(f"=== readout fallback ladder: re-running R1 at ro={ro2} ===")
            rr = run_row(f'R1ro{ro2}', 100e-6, 70e-6, 3.5e-3, ro2, isa, refs, table, nq, log)
            rr['params'] = {'T1_us': 100, 'T2_us': 70, 'eps': 3.5e-3, 'ro': ro2, 'verdict_row': True}
            results[f'R1ro{ro2}'] = rr
            if rr['row_pass']:
                ladder_note = f'R1 passed at tightened readout cap {ro2}; seal {ro2} into rule text'
                r1 = rr
                break
    if 'R2' not in results or 'R3' not in results:
        with open(os.path.join(res_dir, '2026-09-04_nk_v4_row_R1final.json'), 'w') as f:
            json.dump({'r1_final': r1, 'ladder_note': ladder_note}, f, indent=1)
        log(f'R1 branch complete (ladder: {ladder_note}); combiner pending')
        return
    amendment_pass = bool(r1['row_pass'] and results['R2']['row_pass'] and results['R3']['row_pass'])
    verdict = ('AMENDMENT-VALIDATED' if amendment_pass else 'AMENDMENT-DENIED (manifest v1 stands)')
    # re-measured iv-D headroom line at the amended floors (mandatory for manifest v2)
    floor_rows = [results[k] for k in results if results[k]['params']['T1_us'] == 100]
    frac = [r['ivD_mean']/THR_IV for r in floor_rows]
    headroom_line = (f'At the amended floors the iv-D canary runs at '
                     f'{min(frac)*100:.0f}-{max(frac)*100:.0f}% of its threshold '
                     f'(measured in gate v4 across the floor rows); sub-shot-noise trips remain '
                     f'indeterminate by the sealed band rule.')
    log(f"\n=== GATE v4 VERDICT: {verdict}" + (f' | {ladder_note}' if ladder_note else '') + " ===")
    log(headroom_line)

    out = {'_provenance': {'script': 'nk_stage0_gate_v4.py', 'seed': SEED, 'sim_only': True,
                           'qpu_spent': 0, 'sealed_by': 'chair ruling wf_c07c723b-b96 (adopted in full)',
                           'one_amendment_discipline': 'first and LAST floor amendment; a second is '
                           'forbidden threshold-shopping; fallback after v4 = gate-negative under v1'},
           'rows': results, 'verdict': verdict, 'ladder_note': ladder_note,
           'headroom_line_v2': headroom_line, 'runtime_s': round(time.time()-t0, 1)}
    path = os.path.join(os.path.normpath(os.path.join(HERE, '..', 'data')),
                        '2026-09-04_nk_stage0_gate_v4.json')
    with open(path, 'w') as f: json.dump(out, f, indent=1)
    log('wrote ' + path)

if __name__ == '__main__':
    main()
