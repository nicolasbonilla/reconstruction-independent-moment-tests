# -*- coding: utf-8 -*-
"""v3-R2 — the SEALED resolution protocol for the v3 indeterminate (chair ruling wf_e79137d4-be0).

Sealed BEFORE running (verbatim from the ruling):
  - K = 8 new replicas, seeds FIXED and immutable: 2,3,4,5,6,7,8,9. Identical config to the
    binding row (eps=3.5e-3, tw32, manifest floors T1>=150us/T2>=100us, same heavy-hex FakeFez
    ISA build, same calibFT pulse-identical twin, independent shots+dressings, same chain).
  - Replica set = ALL of {seed0, seed1 (from the v3 JSON, exact recorded values), seed2..seed9}.
    Dropping either existing seed is a forbidden fork.
  - Statistic: mean(B_FT) over the 10; SE = sd/sqrt(10). Terminal and symmetric:
      QUOTABLE PASS iff mean + 1.645*SE < 0.0266667  -> STOP lifted, v3 = 10/10 PASS,
                                                        gate PROCEED at lambda*=0.50.
      QUOTABLE FAIL iff mean - 1.645*SE > 0.0266667  -> STOP CONFIRMED (10-replica weight).
      OTHERWISE (straddle)                            -> STOP STANDS by pre-commitment
                                                        (gate-negative side; priced as likely).
  - Mandatory report line (every outcome): at manifest floors the iv-D canary runs at 74-100%
    of its threshold (0.0198/0.0263/0.0267/stress 0.0259) — headroom <=26%, within ~1 sigma of
    tripping; on hardware this canary will trip stochastically unless a gate v4 aligns the
    canary threshold to the claim level (Delta_n(0.50)/5) ex ante.

SIM-ONLY, $0. Writes 06_results/2026-09-02_nk_stage0_v3R2.json.
"""
import os, sys, json, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from nk_stage0_gate import (L, FILL, SHOTS, RO_01, RO_10, ed_references, ft_matrix,
                            build_circuits, counts_pipeline)
from nk_stage0_gate_v2 import cz_pauli_table, build_twirled_instance
from nk_stage0_gate_v3 import build_isa, make_noise, T1_FLOOR, T2_FLOOR
from qiskit.circuit import Delay
from qiskit.transpiler import PassManager
from qiskit_aer import AerSimulator
from qiskit_aer.noise import RelaxationNoisePass

SEED = 20260902
EPS = 3.5e-3
N_TWIRL = 32
THR = 0.20*(2.0/3.0)/5.0                     # 0.0266667 (Delta_n(0.20)/5, mirrorFT maxdev 2/3)
NEW_SEEDS = [2, 3, 4, 5, 6, 7, 8, 9]         # sealed, immutable

def run():
    t0 = time.time()
    refs, _ = ed_references()
    W = ft_matrix().conj().T
    circs = build_circuits(W)
    val_map = {'rung0': 'rung0_nk', 'mirrorFT': 'rung0_nk', 'rungB': 'rungB_nk',
               'calibB2': 'calibB2_nk', 'calibFT': 'rung0_nk'}
    isa, czc, meta, dt = build_isa(circs, refs, val_map)   # identical deterministic ISA build
    print('[isa] rebuilt identically:', {k: v['cz'] for k, v in meta.items()})
    table = cz_pauli_table()
    dur_cz, dur_1q = 84e-9, 32e-9

    # exact recorded seed0/seed1 values from the v3 JSON (mandatory inclusion)
    res_dir = os.path.normpath(os.path.join(HERE, '..', '06_results'))
    v3 = json.load(open(os.path.join(res_dir, '2026-09-02_nk_stage0_gate_v3.json')))
    B_list, provenance = [], {}
    for tag, want in (('eps=0.0035+tw32+floor', 'seed0'), ('eps=0.0035+tw32+seed1+floor', 'seed1')):
        b = v3['noise_results'][tag]['criteria']['iv_D']['B_resid_max']
        B_list.append(b); provenance[want] = b
    print(f'[existing] seed0 B={B_list[0]:.9f}  seed1 B={B_list[1]:.9f}')

    n_active = max(v.num_qubits for v in isa.values())
    relax_pm = PassManager(RelaxationNoisePass(t1s=[T1_FLOOR]*n_active, t2s=[T2_FLOOR]*n_active,
                                               dt=dt, op_types=[Delay]))
    theta_zz = float(np.sqrt(EPS))
    for stag in NEW_SEEDS:
        rng = np.random.default_rng(SEED + 7919*stag + 101)          # v3 formula, floor branch
        sim = AerSimulator(method='density_matrix', noise_model=make_noise(EPS, 'floor', dur_cz, dur_1q),
                           seed_simulator=SEED + 1000*stag + 7)
        out = {}
        for nm_ in ('mirrorFT', 'calibFT'):    # rung0 consumes no draws (0 CZ) — stream identical
            pool = {}
            for _ in range(N_TWIRL):
                ti = relax_pm.run(build_twirled_instance(isa[nm_], theta_zz, rng, table))
                for k_, v_ in sim.run(ti, shots=SHOTS//N_TWIRL).result().get_counts().items():
                    pool[k_] = pool.get(k_, 0) + v_
            nk, keep = counts_pipeline(pool, (RO_01, RO_10))
            out[nm_] = (np.array(nk), keep)
        ref = refs['rung0_nk']; d = ref - FILL
        p_cal = float(np.clip(np.sum((ref - out['calibFT'][0])*d)/np.sum(d*d), 0, 1))
        nk_t = (out['mirrorFT'][0] - p_cal*FILL)/(1.0 - p_cal)
        B = float(np.max(np.abs(nk_t - ref)))
        B_list.append(B); provenance[f'seed{stag}'] = B
        print(f'[seed{stag}] p_calibFT={p_cal:.4f} B_FT={B:.6f} keep={out["mirrorFT"][1]:.2f}')

    B = np.array(B_list)
    mean, sd = float(B.mean()), float(B.std(ddof=1))
    se = sd/np.sqrt(len(B))
    ucb, lcb = mean + 1.645*se, mean - 1.645*se
    if ucb < THR:   verdict = ('QUOTABLE-PASS', 'STOP lifted: v3 = 10/10 PASS; gate PROCEED at lambda*=0.50')
    elif lcb > THR: verdict = ('QUOTABLE-FAIL', 'STOP CONFIRMED with 10-replica weight; gate-negative')
    else:           verdict = ('STRADDLE', 'STOP STANDS by pre-commitment (gate-negative side)')
    print(f'\n=== v3-R2: mean={mean:.6f} sd={sd:.6f} SE={se:.6f} '
          f'[LCB,UCB]=[{lcb:.6f},{ucb:.6f}] vs thr={THR:.7f} -> {verdict[0]} ===')
    print(verdict[1])
    report_line = ('MANDATORY: at manifest floors the iv-D canary runs at 74-100% of its threshold '
                   '(headroom <=26%, within ~1 sigma); on hardware it will trip stochastically '
                   'unless gate v4 aligns the canary threshold to the claim level ex ante.')
    print(report_line)
    out_json = {
        '_provenance': {'script': 'nk_stage0_v3R2.py', 'sim_only': True, 'qpu_spent': 0,
                        'sealed_by': 'chair ruling wf_e79137d4-be0 (rule fixed before any replica ran)',
                        'rule': 'mean+1.645*SE<thr PASS | mean-1.645*SE>thr FAIL | else STOP stands',
                        'threshold': THR, 'seeds': ['seed0(v3)', 'seed1(v3)'] + [f'seed{s}' for s in NEW_SEEDS]},
        'B_FT_replicas': provenance, 'mean': mean, 'sd': sd, 'se': se,
        'lcb': lcb, 'ucb': ucb, 'verdict': verdict, 'mandatory_report_line': report_line,
        'runtime_s': round(time.time()-t0, 1),
    }
    path = os.path.join(res_dir, '2026-09-02_nk_stage0_v3R2.json')
    with open(path, 'w') as f: json.dump(out_json, f, indent=1)
    print('wrote', path)

if __name__ == '__main__':
    run()
