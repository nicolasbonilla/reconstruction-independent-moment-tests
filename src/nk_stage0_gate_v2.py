# -*- coding: utf-8 -*-
"""STAGE-0 GATE v2 — chair-mandated rerun after the adversarial audit (wf_348330ab-00a).

v1 verdict (PROCEED at lambda*=0.35) FELL: the audit proved (a) omitted amplitude damping is
NOT conservative — scheduled-exposure thermal relaxation brings B_resid to 0.035-0.041 vs the
0.0382 threshold, and criterion (iii) is structurally blind to it (the twin shares the damping);
(b) the analytic infinite-twirl channel understates the finite-32-twirl coherent residue (~1/sqrt32);
(c) v1 appended measure_all AFTER transpile — routing can permute measured bits (D9).
v2 blocking fixes (chair R1-R6):
  R1  thermal-relaxation verdict rows: T1=180us/T2=120us scaled by exposure multipliers m in {3,6}
      (audit: realistic wall-clock exposure sits between x3 and x6; sealing requires passing x6).
  R2  FINITE 32-instance Pauli twirl of the coherent-ZZ error (explicit sampled dressings,
      rzz(sqrt(eps)) after each CZ, twirl Paulis as noiseless 'unitary' insertions = hardware
      virtual merging), pooled counts 32 x SHOTS/32.
  R3  measurements added BEFORE transpile (routing-tracked clbits) + ISA-level statevector
      validation via final_index_layout permutation, < 1e-8 vs ED, as a hard gate.
  R4  sigma_shot recomputed from ACTUAL kept counts per config (0.5/sqrt(N_kept)/sqrt(2)).
  R5  provenance deviations_v2 (complete, incl. estimator-selection history and the
      never-quote-lambda-0.20 note).
  R6  seed replicas on the binding twirl configs.
Estimator frozen per R8: calibB2 twin = SOLE criterion-(iii) verdict estimator (not independent —
a maximally-matched noise replica, pulse-identical up to virtual-Z frames); calibB / halfB-scaled /
mirror-survival recorded as non-verdict diagnostics whatever they say. Discarded candidates and
their failure modes are disclosed in provenance (both FAIL the binding config).

SIM-ONLY, $0. Writes 06_results/2026-09-02_nk_stage0_gate_v2.json.
"""
import os, sys, json, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from nk_stage0_gate import (L, U, FN, THETA, PHI, FILL, KS, EPS_K, SHOTS, EPS1Q, RO_01, RO_10,
                            ed_references, ft_matrix, build_circuits, rotation_matrix_of_circuit,
                            orbital_rotation_circuit, nk_from_statevector, counts_pipeline)
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Statevector, Operator
from qiskit.transpiler import CouplingMap
from qiskit_aer import AerSimulator
from qiskit_aer.noise import (NoiseModel, depolarizing_error, ReadoutError, pauli_error,
                              thermal_relaxation_error)

SEED = 20260902
T1_US, T2_US = 180.0, 120.0            # fez-typical; day-of values go into the sealed manifest
DUR_CZ, DUR_1Q = 60e-9, 30e-9
N_TWIRL = 32
LAMBDAS = [0.20, 0.35, 0.50]
CIRCS_LIGHT = ('rung0', 'rung0p', 'mirrorFT', 'rungB', 'calibB', 'calibB2', 'halfB', 'mirrorB')
CIRCS_HEAVY = ('rung0', 'mirrorFT', 'rungB', 'calibB2')     # verdict-relevant only (runtime)

# (eps, mode, seed, t1_mult): mode None=depol | 'raw'=untwirled coherent (record-only)
#                             | 'tw32'=finite-32-twirl coherent + depol
CONFIGS = [
    (1.5e-3,  None,   0, 0),
    (2.85e-3, None,   0, 0), (2.85e-3, None, 1, 0), (2.85e-3, None, 2, 0),
    (3.5e-3,  None,   0, 0), (3.5e-3,  None, 1, 0), (3.5e-3,  None, 2, 0),
    (2.85e-3, 'tw32', 0, 0),
    (3.5e-3,  'tw32', 0, 0), (3.5e-3, 'tw32', 1, 0),
    (3.5e-3,  'tw32', 0, 3), (3.5e-3, 'tw32', 0, 6),     # R1 binding stack: twirl + T1
    (2.85e-3, 'raw',  0, 0),                              # record-only (twirling mandatory)
]

def cz_pauli_table():
    """CZ (P0 x P1) CZ = phase * (P0' x P1'): numeric table over the 16 2q Paulis."""
    from qiskit.quantum_info import Pauli
    CZ = Operator(np.diag([1, 1, 1, -1])).data
    labels = [a+b for a in 'IXYZ' for b in 'IXYZ']
    table = {}
    for lab in labels:
        M = CZ @ Operator(Pauli(lab)).data @ CZ
        for lab2 in labels:
            P2 = Operator(Pauli(lab2)).data
            r = np.trace(P2.conj().T @ M) / 4.0
            if abs(abs(r) - 1.0) < 1e-9:
                table[lab] = lab2      # phase irrelevant for twirling statistics
                break
    return table

def apply_pauli_ops(qc, lab, qpair):
    """Insert Pauli lab = 'AB' on (q0,q1) as NOISELESS unitary insertions (hardware merges
    twirl Paulis into adjacent 1q frames; labeling avoids the noise model's x/sx hooks)."""
    P = {'I': np.eye(2), 'X': np.array([[0,1],[1,0]]),
         'Y': np.array([[0,-1j],[1j,0]]), 'Z': np.diag([1,-1])}
    for ch, q in zip(lab[::-1], qpair):    # label 'AB': A on q1, B on q0 (qiskit order)
        if ch != 'I':
            qc.unitary(P[ch], [q], label='twirlP')

def build_twirled_instance(isa_qc, theta_zz, rng, table):
    """One Pauli-twirl dressing of every 2q gate, with the coherent rzz error kept physical
    (inserted between the gate and the compensating Pauli, as on hardware)."""
    out = QuantumCircuit(*isa_qc.qregs, *isa_qc.cregs)
    for inst in isa_qc.data:
        op, qs, cs = inst.operation, inst.qubits, inst.clbits
        if op.num_qubits == 2 and op.name == 'cz':
            lab = 'IXYZ'[rng.integers(4)] + 'IXYZ'[rng.integers(4)]
            apply_pauli_ops(out, lab, qs)
            out.append(op, qs, cs)
            if theta_zz:
                out.rzz(theta_zz, qs[0], qs[1])
            apply_pauli_ops(out, table[lab], qs)
        else:
            out.append(op, qs, cs)
    return out

def make_noise(eps, t1_mult):
    err2 = depolarizing_error(eps, 2)
    err1 = depolarizing_error(EPS1Q, 1)
    if t1_mult > 0:
        t1, t2 = T1_US*1e-6, T2_US*1e-6
        tre1_cz = thermal_relaxation_error(t1, t2, DUR_CZ*t1_mult)
        tre1_1q = thermal_relaxation_error(t1, t2, DUR_1Q*t1_mult)
        err2 = err2.compose(tre1_cz.tensor(tre1_cz))
        err1 = err1.compose(tre1_1q)
    nm = NoiseModel()
    nm.add_all_qubit_quantum_error(err2, ['cz'])
    nm.add_all_qubit_quantum_error(err1, ['sx', 'x'])
    nm.add_all_qubit_readout_error(ReadoutError([[1-RO_01, RO_01], [RO_10, 1-RO_10]]))
    return nm

def run_gate_v2():
    t0 = time.time()
    refs, _ = ed_references()
    F = ft_matrix()
    # convention + logical validation (identical to v1's validated path)
    W = F.conj().T
    M = rotation_matrix_of_circuit(orbital_rotation_circuit, W)
    assert np.max(np.abs(M - W)) < 1e-9
    circs = build_circuits(W)
    val_map = {'rung0': 'rung0_nk', 'rung0p': 'rung0p_nk', 'mirrorFT': 'rung0_nk',
               'rungB': 'rungB_nk', 'calibB': 'calibB_nk', 'calibB2': 'calibB2_nk',
               'halfB': 'rung0p_nk'}
    for nm_, key in val_map.items():
        err = np.max(np.abs(nk_from_statevector(Statevector.from_instruction(circs[nm_]))
                            - refs[key]))
        assert err < 1e-8, f'logical validation failed {nm_}: {err}'
    print('[validate/logical] all rungs < 1e-8 vs ED  PASS')

    # ---- R3: measure BEFORE transpile; ISA validation via final_index_layout ----
    basis = ['cz', 'rz', 'sx', 'x']
    cm = CouplingMap([[i, i+1] for i in range(11)])
    layouts = [list(range(12)), [0,2,4,6,8,10,1,3,5,7,9,11]]
    isa, cz_counts, isa_val = {}, {}, {}
    for nm_, qc in circs.items():
        qm = qc.copy(); qm.measure_all()
        best = None
        for lay in layouts:
            tq = transpile(qm, coupling_map=cm, basis_gates=basis, optimization_level=3,
                           initial_layout=lay, seed_transpiler=SEED)
            n2 = sum(1 for inst in tq.data if inst.operation.num_qubits == 2)
            if best is None or n2 < best[1]: best = (tq, n2)
        tq, n2 = best
        isa[nm_], cz_counts[nm_] = tq, n2
        # ISA-level statevector validation (strip measures, permute physical->logical)
        tq_nom = tq.remove_final_measurements(inplace=False)
        sv = Statevector.from_instruction(tq_nom)
        perm = tq.layout.final_index_layout()          # logical i -> physical perm[i]
        probs = np.abs(sv.data)**2
        occ_phys = np.zeros(12)
        for idx, p in enumerate(probs):
            if p < 1e-16: continue
            for q in range(12):
                if (idx >> q) & 1: occ_phys[q] += p
        occ_log = np.array([occ_phys[perm[i]] for i in range(12)])
        nk_isa = np.array([(occ_log[k] + occ_log[6+k])/2.0 for k in range(L)])
        if nm_ == 'mirrorB':
            isa_val[nm_] = {'survival': float(probs[0])}
            assert probs[0] > 1 - 1e-8
        else:
            e = float(np.max(np.abs(nk_isa - refs[val_map[nm_]])))
            isa_val[nm_] = {'max_err_vs_ED': e}
            assert e < 1e-8, f'ISA validation failed {nm_}: {e}'
    print('[validate/ISA] all rungs < 1e-8 vs ED after routing-permutation  PASS')
    print('[isa] CZ:', cz_counts)

    table = cz_pauli_table()
    results = {}
    for (eps, mode, stag, t1m) in CONFIGS:
        rng = np.random.default_rng(SEED + 7919*stag + 13*int(t1m))
        sim = AerSimulator(method='density_matrix', noise_model=make_noise(eps, t1m),
                           seed_simulator=SEED + 1000*stag + 17*int(t1m))
        theta_zz = float(np.sqrt(eps)) if mode in ('raw', 'tw32') else 0.0
        which = CIRCS_HEAVY if (mode == 'tw32' or t1m > 0) else CIRCS_LIGHT
        out = {}
        for nm_ in which:
            if mode == 'tw32':
                pool = {}
                for _ in range(N_TWIRL):
                    ti = build_twirled_instance(isa[nm_], theta_zz, rng, table)
                    cts = sim.run(ti, shots=SHOTS//N_TWIRL).result().get_counts()
                    for k_, v_ in cts.items(): pool[k_] = pool.get(k_, 0) + v_
                counts = pool
            elif mode == 'raw':
                # raw = coherent rzz after each CZ, NO twirl dressing (record-only config)
                ti = QuantumCircuit(*isa[nm_].qregs, *isa[nm_].cregs)
                for inst in isa[nm_].data:
                    ti.append(inst.operation, inst.qubits, inst.clbits)
                    if inst.operation.num_qubits == 2 and inst.operation.name == 'cz':
                        ti.rzz(theta_zz, inst.qubits[0], inst.qubits[1])
                counts = sim.run(ti, shots=SHOTS).result().get_counts()
            else:
                counts = sim.run(isa[nm_], shots=SHOTS).result().get_counts()
            if nm_ == 'mirrorB':
                zero = '0'*12
                s_raw = counts.get(zero, 0)/SHOTS
                out[nm_] = {'survival_corr': min(1.0, s_raw/((1-RO_01)**12))}
            else:
                nk, keep = counts_pipeline(counts, (RO_01, RO_10))
                out[nm_] = {'nk_hat': nk.tolist(), 'keep_frac': keep}

        def p_fit(nm_, ref):
            nk = np.array(out[nm_]['nk_hat']); d = ref - FILL
            return float(np.clip(np.sum((ref - nk)*d)/np.sum(d*d), 0, 1))
        pf_B    = p_fit('rungB', refs['rungB_nk'])
        pf_cal2 = p_fit('calibB2', refs['calibB2_nk'])
        pf_FT   = p_fit('mirrorFT', refs['rung0_nk'])
        diags = {}
        if 'calibB' in out: diags['p_fit_calibB'] = p_fit('calibB', refs['calibB_nk'])
        if 'halfB' in out:  diags['p_fit_halfB'] = p_fit('halfB', refs['rung0p_nk'])
        if 'mirrorB' in out: diags['p_mirror_total'] = float(1-np.sqrt(max(out['mirrorB']['survival_corr'],1e-12)))

        # R4: sigma_shot from ACTUAL kept counts of the test circuit
        keepB = max(out['rungB']['keep_frac'], 1e-3)
        sigma_shot = 0.5/np.sqrt(SHOTS*keepB)/np.sqrt(2)
        sigma_p = 0.02
        crit = {}
        # mirrorFT is mitigated with the twin p scaled by the CZ ratio: this is the
        # transfer-model chain test (criterion iv) — an exactly-known rung mitigated by a
        # parameter measured elsewhere.
        p_FT_transfer = pf_cal2*cz_counts['mirrorFT']/max(cz_counts['rungB'], 1)
        for nm_, ref, p_use, lam_list in (('rungB', refs['rungB_nk'], pf_cal2, LAMBDAS),
                                          ('mirrorFT', refs['rung0_nk'], p_FT_transfer, [0.20])):
            nk = np.array(out[nm_]['nk_hat'])
            nk_t = (nk - p_use*FILL)/(1.0 - p_use)
            B = float(np.max(np.abs(nk_t - ref)))
            maxdev = float(np.max(np.abs(ref - FILL)))
            sig_eff = float(np.sqrt((sigma_shot/(1-p_use))**2 + (maxdev/(1-p_use)**2*sigma_p)**2))
            rows = {}
            for lam in lam_list:
                dn = lam*maxdev
                rows[f'{lam:.2f}'] = {'Delta_n': dn, 'B_resid_max': B, 'sigma_eff': sig_eff,
                                      'crit_i': bool(B < dn/5.0),
                                      'crit_ii': bool(dn - B > 1.96*sig_eff)}
            crit[nm_] = {'p_used': float(p_use), 'B_resid_max': B, 'lambdas': rows,
                         'nk_mitigated': nk_t.tolist()}
        crit['crit_iii'] = {'p_hat_calibB2': pf_cal2, 'p_fit_rungB_truth': pf_B,
                            'pass': bool(abs(pf_cal2 - pf_B) < 0.02)}
        crit['crit_iv_mirrorFT_lam020'] = bool(crit['mirrorFT']['lambdas']['0.20']['crit_i']
                                               and crit['mirrorFT']['lambdas']['0.20']['crit_ii'])
        crit['sigma_shot_from_keep'] = float(sigma_shot); crit['keep_frac_rungB'] = float(keepB)
        crit['diagnostics_nonverdict'] = diags
        tag = (f'eps={eps:g}' + (f'+{mode}' if mode else '') + (f'+seed{stag}' if stag else '')
               + (f'+T1x{t1m}' if t1m else ''))
        results[tag] = {'criteria': crit,
                        'raw_nk': {k: v.get('nk_hat') for k, v in out.items() if 'nk_hat' in v},
                        'keeps': {k: v.get('keep_frac') for k, v in out.items() if 'keep_frac' in v}}
        print(f"[{tag}] p2={pf_cal2:.3f} truth={pf_B:.3f} B={crit['rungB']['B_resid_max']:.4f} "
              f"keep={keepB:.2f} iii={crit['crit_iii']['pass']} iv={crit['crit_iv_mirrorFT_lam020']}")

    # ---- verdict: every non-raw config in each family must pass; ladder rises to worst ----
    def verdict_at(tag):
        c = results[tag]['criteria']
        if not c['crit_iv_mirrorFT_lam020']: return ('STOP', None)
        if not c['crit_iii']['pass']: return ('STOP', None)
        for lam in ('0.20', '0.35', '0.50'):
            row = c['rungB']['lambdas'].get(lam)
            if row and row['crit_i'] and row['crit_ii']: return ('PROCEED', float(lam))
        return ('GATE-NEGATIVE', None)
    fam_central = [t for t in results if t.startswith('eps=0.00285') and 'raw' not in t]
    fam_pess    = [t for t in results if t.startswith('eps=0.0035') and 'raw' not in t]
    def worst(tags):
        vs = [verdict_at(t) for t in tags]
        if any(v[0] == 'STOP' for v in vs): return ('STOP', None)
        if any(v[0] == 'GATE-NEGATIVE' for v in vs): return ('GATE-NEGATIVE', None)
        return ('PROCEED', max(v[1] for v in vs))
    v_c, v_p = worst(fam_central), worst(fam_pess)
    if v_c[0] == 'PROCEED' and v_p[0] == 'PROCEED':
        final = ('PROCEED', max(v_c[1], v_p[1]))
    elif 'STOP' in (v_c[0], v_p[0]):
        final = ('STOP-MODEL-CHAIN', None)
    else:
        final = ('GATE-NEGATIVE', None)
    print(f"\n=== STAGE-0 v2 VERDICT: {final[0]}" + (f" at lambda*={final[1]}" if final[1] else "") + " ===")

    out_json = {
        '_provenance': {
            'script': 'nk_stage0_gate_v2.py', 'seed': SEED, 'sim_only': True, 'qpu_spent': 0,
            'supersedes': '2026-09-02_nk_stage0_gate.json (v1: PROCEED@0.35 RULED a model artifact '
                          'by the adversarial audit wf_348330ab-00a — T1 omission non-conservative, '
                          'finite-twirl residue, D9 measure-after-transpile bug)',
            'deviations_v2': [
                'estimator history DISCLOSED: mirror-survival and calibB (Slater-input) were tried and '
                'rejected pre-seal for stated physical defects (over/under-estimation); BOTH would fail '
                'the binding config; calibB2 twin is the sole frozen verdict estimator (R8). The Stage-0 '
                'outcome is a DESIGN result, not a pre-registered test pass; the pre-registered test is '
                'the device run alone.',
                'calibB2 is NOT independent: a maximally-matched noise replica, pulse-identical to rungB '
                'up to virtual-Z frames; criterion (iii) = residual-coherence detector; the falsifier is '
                'insensitive to corruption components proportional to (n_ref - fill) by design.',
                'lambda=0.20 passes are inside shot noise and are NEVER quotable (audit 1).',
                'line-topology ISA (conservative CZ); uniform per-edge eps (day-of chain-acceptance rule '
                'in the manifest closes the tail-edge gap).',
                'symmetric readout inversion as TREX proxy (sim harder than device in asymmetry).',
                'leakage (~3%/circuit) accepted unmodeled; named in the manifest.',
                'T1 exposure modeled as gate-duration x multiplier {3,6} (audit method: realistic '
                'wall-clock exposure between x3 and x6); day-of per-qubit T1/T2 go into the manifest.',
                'v1 JSON rung0p diagnostics may be routing-scrambled (D9); superseded by this file.',
            ],
        },
        'frozen_refs': {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in refs.items()},
        'theta_phi': [THETA, PHI], 'shots': SHOTS, 'n_twirl': N_TWIRL,
        'isa_cz_counts': cz_counts, 'isa_validation': isa_val,
        'noise_results': results,
        'verdict': {'central_family': {t: verdict_at(t) for t in fam_central},
                    'pessimistic_family': {t: verdict_at(t) for t in fam_pess},
                    'raw_coherent_record': {t: verdict_at(t) for t in results if 'raw' in t},
                    'final': final,
                    'note': 'raw coherent-ZZ excluded from verdict families: gate-level Pauli '
                            'twirling (N>=32, independently compiled) is a sealed job requirement '
                            '(R7); if the job cannot guarantee it, the sealed verdict is STOP.'},
        'runtime_s': round(time.time()-t0, 1),
    }
    path = os.path.join(os.path.normpath(os.path.join(HERE, '..', '06_results')),
                        '2026-09-02_nk_stage0_gate_v2.json')
    with open(path, 'w') as f: json.dump(out_json, f, indent=1)
    print('wrote', path, f'({out_json["runtime_s"]}s)')

if __name__ == '__main__':
    run_gate_v2()
