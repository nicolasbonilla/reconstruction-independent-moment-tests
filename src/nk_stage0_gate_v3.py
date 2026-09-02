# -*- coding: utf-8 -*-
"""STAGE-0 GATE v3 — final adjudicated build (chair ruling wf_90a0b575-9ab).

Everything below was FIXED BY THE ADJUDICATION BEFORE ANY v3 NUMBER EXISTED (anti-forking-paths):
  ISA        : heavy-hex FakeFez, sabre layout/routing (seed-scanned), compacted to the 12 active
               qubits; scheduled ASAP with explicit Delays. The v2 line ISA (281 CZ, 12.46 us) is
               superseded — it wasted ~27% wall-clock the hardware would never spend.
  NOISE      : depolarizing (cz eps, 1q 2.5e-4) + asymmetric readout as before; T1/T2 via
               SCHEDULED, DELAY-AWARE relaxation (RelaxationNoisePass on Delays + gate-duration
               thermal relaxation composed into gate errors) at MANIFEST FLOORS T1=150us, T2=100us
               (verdict rows) and a half-floor stress row (75/50, non-verdict). The uniform-
               multiplier abstraction of v2 was STRUCK (per-qubit exposure support 1.85-18.2).
  CRITERION iv (adjudicated semantics D primary, B fallback): calibFT = an independent,
               pulse-identical acquisition of the mirrorFT ISA (own shots, own twirl dressings);
               the deployed chain runs VERBATIM on the known-truth rung: p measured on calibFT,
               applied to mirrorFT, post-mitigation residual vs exact ED, thresholds unchanged
               (crit-i B < Delta_n(0.20)/5 = 0.0267; crit-ii). The v2 CZ-ratio transfer is STRUCK
               (condemned depth-scaling; systematic bias confirmed in all 13 v2 configs) and is
               recorded as a non-verdict diagnostic only. Self-fit (B) recorded as fallback.
  LAMBDA     : sealed claim rule frozen NOW: lambda*=0.50 iff all verdict rows pass at 0.50;
               STRENGTHEN to 0.35 only iff every verdict row passes 0.35 with B < 0.8*threshold
               (pre-registered margin condition). No other upgrade path exists.
  DAY-OF RULE (manifest, verbatim): accept the chain iff min over chain qubits of day-of
               calibrated T1_q >= 150 us AND T2_q >= 100 us; otherwise re-select or forfeit.

SIM-ONLY, $0. Writes 06_results/2026-09-02_nk_stage0_gate_v3.json.
"""
import os, sys, json, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from nk_stage0_gate import (L, U, FN, THETA, PHI, FILL, KS, EPS_K, SHOTS, EPS1Q, RO_01, RO_10,
                            ed_references, ft_matrix, build_circuits,
                            rotation_matrix_of_circuit, orbital_rotation_circuit,
                            nk_from_statevector, counts_pipeline)
from nk_stage0_gate_v2 import cz_pauli_table, build_twirled_instance
from qiskit import QuantumCircuit, transpile
from qiskit.circuit import Delay
from qiskit.quantum_info import Statevector
from qiskit.transpiler import PassManager
from qiskit_aer import AerSimulator
from qiskit_aer.noise import (NoiseModel, depolarizing_error, ReadoutError,
                              thermal_relaxation_error, RelaxationNoisePass)
from qiskit_ibm_runtime.fake_provider import FakeFez

SEED = 20260902
N_TWIRL = 32
T1_FLOOR, T2_FLOOR = 150e-6, 100e-6          # manifest floors (verdict rows)
T1_STRESS, T2_STRESS = 75e-6, 50e-6          # half-floor stress (non-verdict)
LAMBDAS = ['0.35', '0.50']                   # sealed ladder (0.20 never quotable — audit 1)
CIRCS = ('rung0', 'mirrorFT', 'calibFT', 'rungB', 'calibB2')
# (eps, mode, seed, relax): relax in (None, 'floor', 'stress')
CONFIGS = [
    (1.5e-3,  None,   0, None), (2.85e-3, None, 0, None), (2.85e-3, None, 1, None),
    (3.5e-3,  None,   0, None), (3.5e-3,  None, 1, None),
    (3.5e-3,  'tw32', 0, None),
    (2.85e-3, 'tw32', 0, 'floor'),
    (3.5e-3,  'tw32', 0, 'floor'), (3.5e-3, 'tw32', 1, 'floor'),   # BINDING verdict stack
    (3.5e-3,  'tw32', 0, 'stress'),                                 # non-verdict
    (2.85e-3, 'raw',  0, None),                                     # record-only
]
NONVERDICT = {'stress', 'raw'}

def compact(tq):
    """Reduce a 156q transpiled circuit to its active qubits (ASAP scheduling delays idle
    qubits too — active = qubits touched by any NON-Delay op; idle-qubit Delays dropped)."""
    act = sorted({tq.find_bit(q).index for inst in tq.data
                  if inst.operation.name not in ('delay', 'barrier')
                  for q in inst.qubits})
    actset = set(act)
    cmap = {p: i for i, p in enumerate(act)}
    out = QuantumCircuit(len(act), tq.num_clbits)
    for inst in tq.data:
        idxs = [tq.find_bit(q).index for q in inst.qubits]
        if inst.operation.name in ('delay', 'barrier') and not all(i in actset for i in idxs):
            continue
        qs = [out.qubits[cmap[i]] for i in idxs]
        cs = [out.clbits[tq.find_bit(c).index] for c in inst.clbits]
        out.append(inst.operation, qs, cs)
    return out, cmap

def build_isa(circs, refs, val_map):
    fez = FakeFez(); dt = fez.target.dt
    isa, cz_counts, meta = {}, {}, {}
    for nm_, qc in circs.items():
        if nm_ in ('rung0p', 'halfB', 'calibB', 'mirrorB'): continue   # not in v3 circuit set
        qm = qc.copy(); qm.measure_all()
        best = None
        for seed in range(5):
            tq = transpile(qm, backend=fez, optimization_level=3, seed_transpiler=seed,
                           scheduling_method='asap')
            n2 = sum(1 for i in tq.data if i.operation.num_qubits == 2)
            if best is None or n2 < best[1]: best = (tq, n2, seed)
        tq, n2, sd = best
        cq, cmap = compact(tq)
        # logical->compact permutation for validation
        fil = tq.layout.final_index_layout()
        perm = [cmap[fil[i]] for i in range(12)]
        # ISA statevector validation (strip measures; Delays are identity)
        nom = cq.remove_final_measurements(inplace=False)
        sv = Statevector.from_instruction(nom)
        probs = np.abs(sv.data)**2
        nq = cq.num_qubits
        occ = np.zeros(nq)
        for idx, p in enumerate(probs):
            if p < 1e-16: continue
            for q in range(nq):
                if (idx >> q) & 1: occ[q] += p
        occ_log = np.array([occ[perm[i]] for i in range(12)])
        nk_isa = np.array([(occ_log[k] + occ_log[6+k])/2.0 for k in range(L)])
        e = float(np.max(np.abs(nk_isa - refs[val_map[nm_]])))
        assert e < 1e-8, f'ISA validation failed {nm_}: {e}'
        # wall-clock from schedule
        wall_dt = tq.duration if tq.duration else 0
        isa[nm_], cz_counts[nm_] = cq, n2
        meta[nm_] = {'cz': n2, 'seed': sd, 'active_qubits': cq.num_qubits,
                     'wall_us': round(wall_dt*dt*1e6, 3) if wall_dt else None,
                     'isa_val_err': e}
    # calibFT = independent pulse-identical acquisition of the mirrorFT ISA (criterion iv-D)
    isa['calibFT'] = isa['mirrorFT'].copy()
    cz_counts['calibFT'] = cz_counts['mirrorFT']
    meta['calibFT'] = dict(meta['mirrorFT']); meta['calibFT']['twin_of'] = 'mirrorFT'
    return isa, cz_counts, meta, dt

def make_noise(eps, relax, dur_cz, dur_1q):
    err2 = depolarizing_error(eps, 2)
    err1 = depolarizing_error(EPS1Q, 1)
    if relax:
        t1, t2 = (T1_FLOOR, T2_FLOOR) if relax == 'floor' else (T1_STRESS, T2_STRESS)
        g2 = thermal_relaxation_error(t1, t2, dur_cz)
        err2 = err2.compose(g2.tensor(g2))
        err1 = err1.compose(thermal_relaxation_error(t1, t2, dur_1q))
    nm = NoiseModel()
    nm.add_all_qubit_quantum_error(err2, ['cz'])
    nm.add_all_qubit_quantum_error(err1, ['sx', 'x'])
    nm.add_all_qubit_readout_error(ReadoutError([[1-RO_01, RO_01], [RO_10, 1-RO_10]]))
    return nm

def run_gate_v3():
    t0 = time.time()
    refs, _ = ed_references()
    W = ft_matrix().conj().T
    assert np.max(np.abs(rotation_matrix_of_circuit(orbital_rotation_circuit, W) - W)) < 1e-9
    circs = build_circuits(W)
    val_map = {'rung0': 'rung0_nk', 'mirrorFT': 'rung0_nk', 'rungB': 'rungB_nk',
               'calibB2': 'calibB2_nk', 'calibFT': 'rung0_nk'}
    for nm_ in ('rung0', 'mirrorFT', 'rungB', 'calibB2'):
        e = np.max(np.abs(nk_from_statevector(Statevector.from_instruction(circs[nm_]))
                          - refs[val_map[nm_]]))
        assert e < 1e-8
    print('[validate/logical] PASS')
    isa, czc, meta, dt = build_isa(circs, refs, val_map)
    print('[validate/ISA] PASS on FakeFez heavy-hex (all < 1e-8)')
    print('[isa]', {k: (v['cz'], v['wall_us']) for k, v in meta.items()})
    dur_cz = 84e-9; dur_1q = 32e-9        # fez-typical (cz 21 dt, 1q 8 dt at dt=4ns)

    table = cz_pauli_table()
    results = {}
    for (eps, mode, stag, relax) in CONFIGS:
        rng = np.random.default_rng(SEED + 7919*stag + (101 if relax else 0))
        sim = AerSimulator(method='density_matrix', noise_model=make_noise(eps, relax, dur_cz, dur_1q),
                           seed_simulator=SEED + 1000*stag + (7 if relax == 'floor' else 11 if relax == 'stress' else 0))
        theta_zz = float(np.sqrt(eps)) if mode in ('raw', 'tw32') else 0.0
        relax_pm = None
        if relax:
            t1, t2 = (T1_FLOOR, T2_FLOOR) if relax == 'floor' else (T1_STRESS, T2_STRESS)
            n_active = max(v.num_qubits for v in isa.values())
            relax_pm = PassManager(RelaxationNoisePass(
                t1s=[t1]*n_active, t2s=[t2]*n_active, dt=dt, op_types=[Delay]))
        out = {}
        for nm_ in CIRCS:
            base = isa[nm_]
            def prep_instance(qc_):
                return relax_pm.run(qc_) if relax_pm else qc_
            if mode == 'tw32':
                pool = {}
                for _ in range(N_TWIRL):
                    ti = build_twirled_instance(base, theta_zz, rng, table)
                    ti = prep_instance(ti)
                    for k_, v_ in sim.run(ti, shots=SHOTS//N_TWIRL).result().get_counts().items():
                        pool[k_] = pool.get(k_, 0) + v_
                counts = pool
            elif mode == 'raw':
                ti = QuantumCircuit(*base.qregs, *base.cregs)
                for inst in base.data:
                    ti.append(inst.operation, inst.qubits, inst.clbits)
                    if inst.operation.num_qubits == 2 and inst.operation.name == 'cz':
                        ti.rzz(theta_zz, inst.qubits[0], inst.qubits[1])
                counts = sim.run(prep_instance(ti), shots=SHOTS).result().get_counts()
            else:
                counts = sim.run(prep_instance(base.copy()), shots=SHOTS).result().get_counts()
            nk, keep = counts_pipeline(counts, (RO_01, RO_10))
            out[nm_] = {'nk_hat': nk.tolist(), 'keep_frac': keep}

        def p_fit(nm_, refkey):
            nk = np.array(out[nm_]['nk_hat']); ref = refs[refkey]; d = ref - FILL
            return float(np.clip(np.sum((ref - nk)*d)/np.sum(d*d), 0, 1))
        pf_B, pf_cal2 = p_fit('rungB', 'rungB_nk'), p_fit('calibB2', 'calibB2_nk')
        pf_FT_twin = p_fit('calibFT', 'rung0_nk')     # iv-D: measured on the pulse-identical twin
        pf_FT_self = p_fit('mirrorFT', 'rung0_nk')    # iv-B fallback (recorded)
        keepB = max(out['rungB']['keep_frac'], 1e-3)
        sigma_shot = 0.5/np.sqrt(SHOTS*keepB)/np.sqrt(2)
        sigma_p = 0.02

        def resid(nm_, refkey, p_use):
            nk = np.array(out[nm_]['nk_hat']); ref = refs[refkey]
            nk_t = (nk - p_use*FILL)/(1.0 - p_use)
            return nk_t, float(np.max(np.abs(nk_t - ref)))
        crit = {}
        # headline: rungB mitigated by the calibB2 twin
        nkB, BB = resid('rungB', 'rungB_nk', pf_cal2)
        maxdevB = float(np.max(np.abs(refs['rungB_nk'] - FILL)))
        sigB = float(np.sqrt((sigma_shot/(1-pf_cal2))**2 + (maxdevB/(1-pf_cal2)**2*sigma_p)**2))
        rows = {}
        for lam_s in LAMBDAS:
            lam = float(lam_s); dn = lam*maxdevB
            rows[lam_s] = {'Delta_n': dn, 'B_resid_max': BB, 'sigma_eff': sigB,
                           'crit_i': bool(BB < dn/5.0), 'crit_ii': bool(dn - BB > 1.96*sigB),
                           'strengthen_margin_ok': bool(BB < 0.8*dn/5.0)}
        crit['rungB'] = {'p_used': pf_cal2, 'B_resid_max': BB, 'lambdas': rows}
        # iv-D: mirrorFT mitigated by calibFT-measured p (deployed chain verbatim)
        nkF, BF = resid('mirrorFT', 'rung0_nk', pf_FT_twin)
        maxdevF = float(np.max(np.abs(refs['rung0_nk'] - FILL)))
        dnF = 0.20*maxdevF
        sigF = float(np.sqrt((sigma_shot/(1-pf_FT_twin))**2 + (maxdevF/(1-pf_FT_twin)**2*sigma_p)**2))
        crit['iv_D'] = {'p_calibFT': pf_FT_twin, 'B_resid_max': BF, 'threshold': dnF/5.0,
                        'crit_i': bool(BF < dnF/5.0), 'crit_ii': bool(dnF - BF > 1.96*sigF)}
        _, BFs = resid('mirrorFT', 'rung0_nk', pf_FT_self)
        crit['iv_B_fallback'] = {'p_selffit': pf_FT_self, 'B_resid_max': BFs,
                                 'crit_i': bool(BFs < dnF/5.0)}
        crit['iv_A_diagnostic_nonverdict'] = {'p_czratio': pf_cal2*czc['mirrorFT']/czc['rungB']}
        crit['crit_iii'] = {'p_hat_calibB2': pf_cal2, 'p_fit_rungB_truth': pf_B,
                            'pass': bool(abs(pf_cal2 - pf_B) < 0.02)}
        crit['sigma_shot_from_keep'] = sigma_shot; crit['keep_frac_rungB'] = keepB
        tag = (f'eps={eps:g}' + (f'+{mode}' if mode else '') + (f'+seed{stag}' if stag else '')
               + (f'+{relax}' if relax else ''))
        results[tag] = {'criteria': crit, 'raw_nk': {k: v['nk_hat'] for k, v in out.items()},
                        'keeps': {k: v['keep_frac'] for k, v in out.items()},
                        'nonverdict': bool(relax == 'stress' or mode == 'raw')}
        iv_ok = crit['iv_D']['crit_i'] and crit['iv_D']['crit_ii']
        print(f"[{tag}] p2={pf_cal2:.3f} truth={pf_B:.3f} B(B)={BB:.4f} "
              f"ivD:B={BF:.4f}({'ok' if iv_ok else 'FAIL'}) iii={crit['crit_iii']['pass']} keep={keepB:.2f}")

    # ---- sealed verdict (rule frozen in the adjudication, before any v3 number) ----
    verdict_rows = [t for t, r in results.items() if not r['nonverdict']]
    def row_state(t):
        c = results[t]['criteria']
        if not c['crit_iii']['pass']: return 'STOP'
        if not (c['iv_D']['crit_i'] and c['iv_D']['crit_ii']): return 'STOP'
        return 'OK'
    stops = [t for t in verdict_rows if row_state(t) == 'STOP']
    def all_pass(lam_s, margin=False):
        key = 'strengthen_margin_ok' if margin else 'crit_i'
        return all(results[t]['criteria']['rungB']['lambdas'][lam_s][key]
                   and results[t]['criteria']['rungB']['lambdas'][lam_s]['crit_ii']
                   for t in verdict_rows)
    if stops:
        final = ('STOP-MODEL-CHAIN', None, stops)
    elif all_pass('0.50'):
        lam_seal = 0.35 if all_pass('0.35', margin=True) else 0.50
        final = ('PROCEED', lam_seal, [])
    else:
        final = ('GATE-NEGATIVE', None, [])
    print(f"\n=== STAGE-0 v3 VERDICT: {final[0]}"
          + (f" — SEAL lambda*={final[1]}" if final[1] else "")
          + (f"  (stops: {final[2]})" if final[2] else "") + " ===")

    out_json = {
        '_provenance': {
            'script': 'nk_stage0_gate_v3.py', 'seed': SEED, 'sim_only': True, 'qpu_spent': 0,
            'adjudication': 'wf_90a0b575-9ab chair ruling: ISA=FakeFez heavy-hex sabre (line ISA '
                            'superseded); T1 multipliers STRUCK -> scheduled delay-aware relaxation '
                            'at manifest floors T1>=150us T2>=100us (+half-floor stress, non-verdict); '
                            'criterion iv = D (calibFT pulse-identical twin, deployed chain verbatim) '
                            'with B self-fit fallback, CZ-ratio transfer STRUCK to diagnostic; '
                            'sealed lambda ladder frozen pre-run: 0.50 primary, strengthen to 0.35 '
                            'only if all verdict rows pass with B < 0.8*threshold.',
            'day_of_chain_rule': 'accept iff min over chain qubits of day-of T1_q >= 150us AND '
                                 'T2_q >= 100us; else re-select chain or forfeit the window.',
            'v2_record': 'v1 PROCEED@0.35 fell (T1 omission); v2 STOP stands on the line ISA — '
                         'attributed to genuine T1-model breakdown at realistic exposure (adjudicated).',
        },
        'frozen_refs': {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in refs.items()},
        'isa_meta': meta, 'shots': SHOTS, 'n_twirl': N_TWIRL,
        'noise_results': results,
        'verdict': {'rows': {t: row_state(t) for t in verdict_rows}, 'final': final},
        'runtime_s': round(time.time()-t0, 1),
    }
    path = os.path.join(os.path.normpath(os.path.join(HERE, '..', '06_results')),
                        '2026-09-02_nk_stage0_gate_v3.json')
    with open(path, 'w') as f: json.dump(out_json, f, indent=1)
    print('wrote', path, f'({out_json["runtime_s"]}s)')

if __name__ == '__main__':
    run_gate_v3()
