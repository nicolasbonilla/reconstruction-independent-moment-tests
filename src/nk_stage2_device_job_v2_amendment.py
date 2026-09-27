# -*- coding: utf-8 -*-
"""AMENDMENT VERSION (day_of_rule_v2, 2026-09-04). The amendment was DENIED on 2026-09-05
(data/2026-09-05_nk_v4_verdict.json); manifest v2 was never sealed; this job was never run.
The job that matches the sealed manifest v1 is nk_stage2_device_job.py. One editing artifact was
fixed on publication (2026-09-26): a literal "\n" in the manifest-v2 assert, which made the file fail
to parse, is now a line continuation; nothing else in the logic was changed. On 2026-09-27 the docstring,
the printed messages and the assert message were rescoped to state that manifest v2 was never sealed;
the logic is unchanged.

STAGE 2 — the on-device n_k discriminating-falsifier job, amendment version (never run).

Reads the circuits, thresholds and analysis from the sealed manifest v1 (data/manifest_nk_device_v1.json,
SHA-256 recorded), but submits shots only if a manifest v2 (data/manifest_nk_device_v2.json) exists. Manifest
v2 was NEVER written or sealed, because the amendment it would have recorded was denied by its gate
(nk_stage0_gate_v4_combine.py -> data/2026-09-05_nk_v4_verdict.json, AMENDMENT-DENIED), so this job can
only dry-run. It applies day_of_rule_v2 (pick_chain_v2 below), the day-of rule PROPOSED by the amendment
(fixed before the gate-v4 verdict rows ran), not the sealed v1 rule: backends ibm_fez / ibm_marrakesh /
ibm_kingston with CZ median <= 3.5e-3; a 12-qubit path with T1 >= 100 us, T2 >= 70 us and readout <= ro_cap
on every qubit, and edge CZ <= min(3.5e-3, 2x median). Without a manifest v2, ro_cap = 0.030, the cap the
amendment proposed and the first rung of its readout ladder; the gate's ladder went 0.030 -> 0.025 -> 0.018
and the 0.018 row still failed, so no readout cap was ever adopted.
Modes:
  default        : DRY RUN — loads credentials, applies the day-of chain-acceptance rule, transpiles
                   onto the accepted chain, prints the plan + estimated QPU time. SUBMITS NOTHING.
  RUN=1          : would submit the Batch (~3-4 QPU min of the 10-min window), but ABORTS at the manifest-v2
                   assert, since manifest v2 does not exist. Nothing was ever submitted.
  ANALYZE=<file> : runs the frozen analysis on a retained-counts JSON (device or dry-run replay).

Credentials: QiskitRuntimeService.save_account(channel='ibm_quantum_platform', token=<YOUR NEW KEY>,
instance=<YOUR CRN>) once, locally — NEVER commit a token; revoke any previously exposed token first.

Job spec (from the sealed manifest v1): 5 circuits (rung0, mirrorFT, calibFT, rungB, calibB2), 50k shots each,
GATE-LEVEL Pauli twirling num_randomizations=32 (hard requirement) + measurement twirling (the manifest's
wording 'TREX/measure twirling' is a misnomer: SamplerV2 returns twirled raw counts and no readout-error
mitigation is applied) + DD XpXm, one Batch, randomized interleaving. (The v1 day-of chain rule, min T1 >= 150us AND min T2 >= 100us
AND no chain CZ error > 2x device median, is applied by nk_stage2_device_job.py, not by this file.)
"""
import os, sys, json, time, hashlib, io
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
RES = os.path.normpath(os.path.join(HERE, '..', 'data'))
from nk_stage0_gate import (L, FILL, SHOTS, RO_01, RO_10, ed_references, ft_matrix,
                            build_circuits, counts_pipeline)

MANIFEST = json.load(open(os.path.join(RES, 'manifest_nk_device_v1.json')))
MHASH = open(os.path.join(RES, 'manifest_nk_device_v1.sha256')).read().strip()
CIRC_NAMES = MANIFEST['job_spec']['circuits']
TAU = MANIFEST['criteria']['falsifier_decision']['tau']
THR_I = 0.20*(2.0/3.0)/5.0   # iv threshold (band rule per manifest)

def frozen_logical_circuits():
    refs, _ = ed_references()
    W = ft_matrix().conj().T
    circs = build_circuits(W)
    circs['calibFT'] = circs['mirrorFT'].copy()
    out = {}
    for nm in CIRC_NAMES:
        qc = circs[nm].copy(); qc.measure_all()
        out[nm] = qc
    return refs, out

def pick_chain_v2(svc, ro_cap=0.030):
    """day_of_rule_v2 (proposed by the denied amendment; fixed before the gate-v4 verdict rows ran; never
    sealed into a manifest): backends B={fez,marrakesh,kingston};
    qualification CZmed<=3.5e-3; chain = 12-path with T1>=100us, T2>=70us, RO<=ro_cap per
    qubit, edge CZ<=min(3.5e-3, 2*med); deterministic selection: min sum CZ, ties by
    (1) min sum RO, (2) min sum 1/T2, (3) lexicographic. Returns (backend_name, chain,
    score, snapshot) or (None,...) = FORFEIT."""
    import collections
    B = ('ibm_fez', 'ibm_marrakesh', 'ibm_kingston')
    def safe(fn, *a):
        try: return fn(*a)
        except Exception: return None
    best = None; snapshots = {}
    for bname in B:
        bk = svc.backend(bname); props = bk.properties()
        cz = {}
        for (x, y) in bk.coupling_map:
            e = safe(props.gate_error, 'cz', [x, y])
            if e is not None: cz[(x, y)] = e
        med = float(np.median(list(cz.values())))
        snapshots[bname] = {'cal_time': str(props.last_update_date), 'cz_median': med}
        if med > 3.5e-3: continue                      # backend disqualified
        cap = min(3.5e-3, 2*med)
        T1 = {q: safe(props.t1, q) for q in range(bk.num_qubits)}
        T2 = {q: safe(props.t2, q) for q in range(bk.num_qubits)}
        RO = {q: safe(props.readout_error, q) for q in range(bk.num_qubits)}
        ok = {q for q in range(bk.num_qubits)
              if (T1[q] or 0) >= 100e-6 and (T2[q] or 0) >= 70e-6 and (RO[q] or 1) <= ro_cap}
        adj = collections.defaultdict(set)
        def ee(a, b): return cz.get((a, b), cz.get((b, a), np.nan))
        for (x, y), e in cz.items():
            if x in ok and y in ok and e <= cap: adj[x].add(y); adj[y].add(x)
        found = []
        def dfs(path, vis):
            if len(path) == 12:
                found.append(list(path)); return len(found) >= 4000
            for n in sorted(adj[path[-1]]):
                if n not in vis:
                    vis.add(n); path.append(n)
                    if dfs(path, vis): return True
                    path.pop(); vis.remove(n)
            return False
        for st in sorted(ok):
            if dfs([st], {st}): break
        for ch in found:
            sc = (sum(ee(ch[i], ch[i+1]) for i in range(11)),
                  sum(RO[q] for q in ch), sum(1.0/T2[q] for q in ch), tuple(ch))
            if best is None or sc < best[2]:
                best = (bname, ch, sc)
    if best is None: return None, None, None, snapshots
    return best[0], best[1], best[2], snapshots

def analyze(counts_by_circuit, refs, out_path):
    res = {}
    for nm in CIRC_NAMES:
        nk, keep = counts_pipeline(counts_by_circuit[nm], (RO_01, RO_10))
        res[nm] = {'nk_hat': list(map(float, nk)), 'keep': keep}
    ref_B = refs['rungB_nk']; ref_FT = refs['rung0_nk']
    def p_fit(nk, ref):
        d = ref - FILL; return float(np.clip(np.sum((ref - np.array(nk))*d)/np.sum(d*d), 0, 1))
    p2 = p_fit(res['calibB2']['nk_hat'], refs['calibB2_nk'])
    ntilde = (np.array(res['rungB']['nk_hat']) - p2*FILL)/(1 - p2)
    D_honest = float(np.max(np.abs(ntilde - ref_B)))
    lam = MANIFEST['experiment']['claim_level_lambda']
    nbar_corr = (1 - lam)*ref_B + lam*FILL
    D_corr = float(np.max(np.abs(ntilde - nbar_corr)))
    pFT = p_fit(res['calibFT']['nk_hat'], ref_FT)
    ivB = float(np.max(np.abs((np.array(res['mirrorFT']['nk_hat']) - pFT*FILL)/(1 - pFT) - ref_FT)))
    Ttilde = 2.0*float(np.sum((-2*np.cos(2*np.pi*np.arange(L)/L))*ntilde))
    verdict = {
        'honest_accepted': bool(D_honest < TAU), 'corruption_rejected': bool(D_corr > TAU),
        'D_honest': D_honest, 'D_corrupted_lam': D_corr, 'tau': TAU,
        'p_calibB2': p2, 'canary_ivD_B': ivB, 'canary_threshold': THR_I,
        'canary_state': ('pass' if ivB < THR_I - 1.96*2e-3 else
                         'fail-above-band' if ivB > THR_I + 1.96*2e-3 else 'indeterminate-in-band'),
        'kinetic_anchor_Ttilde': Ttilde, 'Ttilde_ED': refs['rungB_Ttilde'],
        'keep_rungB': res['rungB']['keep'],
        'kill_keep_frac': bool(res['rungB']['keep'] < 0.25),
        'FIRES': bool(D_honest < TAU and D_corr > TAU and res['rungB']['keep'] >= 0.25),
    }
    out = {'manifest_sha256': MHASH, 'results': res, 'verdict': verdict,
           'utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
    json.dump(out, open(out_path, 'w'), indent=1)
    print(json.dumps(verdict, indent=1))
    print('wrote', out_path)
    return out

def main():
    print(f'[manifest] {MHASH}')
    refs, logical = frozen_logical_circuits()
    if os.environ.get('ANALYZE'):
        blob = json.load(open(os.environ['ANALYZE']))
        analyze(blob['counts_by_circuit'], refs,
                os.path.join(RES, time.strftime('%Y-%m-%d') + '_nk_device_analysis.json'))
        return
    from qiskit_ibm_runtime import QiskitRuntimeService, Batch, SamplerV2
    svc = QiskitRuntimeService()
    v2p = os.path.join(RES, 'manifest_nk_device_v2.json')
    ro_cap = 0.030
    if os.path.exists(v2p):
        ro_cap = json.load(open(v2p)).get('noise_floors', {}).get('readout_cap', 0.030)
    else:
        print('WARNING: manifest v2 was never sealed (the amendment was denied) — shots are FORBIDDEN (dry run only).')
    bname, chain, score, snaps = pick_chain_v2(svc, ro_cap)
    if bname is None:
        print('day_of_rule_v2: NO candidate chain on any qualifying backend -> FORFEIT (per the proposed rule).')
        print('snapshots:', json.dumps(snaps)); return
    backend = svc.backend(bname)
    snap_path = os.path.join(RES, time.strftime('%Y-%m-%d_%H%M') + '_chain_snapshot.json')
    json.dump({'backend': bname, 'chain': chain, 'score': score[:3], 'snapshots': snaps,
               'ro_cap': ro_cap}, open(snap_path, 'w'), indent=1)
    print(f'[chain] {bname} {chain} sumCZ={score[0]:.4f}; snapshot -> {snap_path}')
    from qiskit import transpile
    isa_jobs = {}
    for nm, qc in logical.items():
        best_t = None
        for seed in range(5):
            tq = transpile(qc, backend=backend, initial_layout=None, optimization_level=3,
                           seed_transpiler=seed)
            n2 = sum(1 for i in tq.data if i.operation.num_qubits == 2)
            if best_t is None or n2 < best_t[1]: best_t = (tq, n2)
        isa_jobs[nm] = best_t[0]
        print(f'  {nm}: {best_t[1]} 2q gates')
    if os.environ.get('RUN') != '1':
        print('\nDRY RUN complete (nothing submitted). RUN=1 would abort: manifest v2 was never sealed '
              f'(planned Batch, ~3-4 QPU min: 5 circuits x {SHOTS} shots, gate twirling N=32, measurement twirling, DD).')
        return
    assert os.path.exists(os.path.join(RES, 'manifest_nk_device_v2.json')), \
        'ABORT: shots forbidden without a sealed manifest v2 (the amendment was denied; manifest v1 stands)'
    # ---- HARD $0 GUARD (standing user constraint: strictly within the free tier) ----
    insts = svc.instances()
    assert any(i.get('pricing_type') == 'free' or i.get('plan') == 'open' for i in insts), \
        'ABORT: no free/open instance found — this job may only run on the $0 plan'
    u = svc.usage()
    rem = u['usage_remaining_seconds']
    EST_S = 110   # measured anchor: 0.283 s/kshot x 250 kshots ~ 71 s, + batch overhead margin
    assert rem >= EST_S + 20, \
        f'ABORT: only {rem}s free-tier QPU remaining < estimate {EST_S}s + margin — never exceed $0'
    print(f'[cost-guard] free plan OK; remaining {rem}s >= {EST_S+20}s needed — $0 assured')
    # ---- SUBMISSION (unreachable: manifest v2 was never sealed) ----
    order = list(CIRC_NAMES); np.random.default_rng(20260902).shuffle(order)
    with Batch(backend=backend) as batch:
        sampler = SamplerV2(mode=batch)
        sampler.options.twirling.enable_gates = True          # HARD requirement (manifest)
        sampler.options.twirling.enable_measure = True        # measurement twirling
        sampler.options.twirling.num_randomizations = 32
        sampler.options.twirling.shots_per_randomization = SHOTS // 32
        sampler.options.dynamical_decoupling.enable = True
        sampler.options.dynamical_decoupling.sequence_type = 'XpXm'
        jobs = {nm: sampler.run([isa_jobs[nm]], shots=SHOTS) for nm in order}
        print('[submitted]', {nm: j.job_id() for nm, j in jobs.items()})
        counts = {}
        for nm, j in jobs.items():
            r = j.result()[0]
            counts[nm] = {k: int(v) for k, v in r.data.meas.get_counts().items()}
    raw_path = os.path.join(RES, time.strftime('%Y-%m-%d') + '_nk_device_counts.json')
    json.dump({'manifest_sha256': MHASH, 'chain': chain,
               'job_ids': {nm: j.job_id() for nm, j in jobs.items()},
               'counts_by_circuit': counts}, open(raw_path, 'w'), indent=1)
    print('retained raw counts ->', raw_path)
    analyze(counts, refs, os.path.join(RES, time.strftime('%Y-%m-%d') + '_nk_device_analysis.json'))

if __name__ == '__main__':
    main()
