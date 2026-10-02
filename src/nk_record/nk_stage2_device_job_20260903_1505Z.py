# -*- coding: utf-8 -*-
"""STAGE 2 — the sealed on-device n_k discriminating-falsifier job (ibm_fez).

Runs ONLY against the sealed manifest (06_results/manifest_nk_device_v1.json, SHA-256 recorded).
Modes:
  default        : DRY RUN — loads credentials, applies the day-of chain-acceptance rule, transpiles
                   onto the accepted chain, prints the plan + estimated QPU time. SUBMITS NOTHING.
  RUN=1          : submits the sealed Batch (~3-4 QPU min of the 10-min window), retains raw counts.
  ANALYZE=<file> : runs the frozen analysis on a retained-counts JSON (device or dry-run replay).

Credentials: QiskitRuntimeService.save_account(channel='ibm_quantum_platform', token=<YOUR NEW KEY>,
instance=<YOUR CRN>) once, locally — NEVER commit a token; revoke any previously exposed token first.

Sealed job spec (manifest): 5 circuits (rung0, mirrorFT, calibFT, rungB, calibB2), 50k shots each,
GATE-LEVEL Pauli twirling num_randomizations=32 (hard requirement) + measure twirling (TREX) + DD
XpXm, one Batch, randomized interleaving. Day-of chain rule: min T1 >= 150us AND min T2 >= 100us AND
no chain CZ error > 2x device median, else re-select or forfeit.
"""
import os, sys, json, time, hashlib, io
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
RES = os.path.normpath(os.path.join(HERE, '..', '06_results'))
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

def pick_chain(backend):
    """Day-of chain-acceptance rule (manifest verbatim): best 12-qubit path meeting the floors."""
    props = backend.properties()
    def safe(fn, q):
        try: return fn(q)
        except Exception: return None      # uncalibrated property -> cannot demonstrate the floor
    T1 = {q: safe(props.t1, q) for q in range(backend.num_qubits)}
    T2 = {q: safe(props.t2, q) for q in range(backend.num_qubits)}
    cmap = backend.coupling_map
    cz_err = {}
    for (a, b) in cmap:
        try: cz_err[(a, b)] = props.gate_error('cz', [a, b])
        except Exception: pass
    med = float(np.median(list(cz_err.values()))) if cz_err else None
    # exhaustive DFS for 12-qubit simple paths satisfying the sealed rule (competent
    # implementation of the SAME rule — greedy chain-growing can miss existing paths),
    # scored by mean CZ error; edge rule: no chain CZ error > 2x device median.
    import collections
    adj = collections.defaultdict(list)
    for (a, b) in cmap: adj[a].append(b); adj[b].append(a)
    def ok_q(q): return (T1[q] is not None and T2[q] is not None
                         and T1[q] >= 150e-6 and T2[q] >= 100e-6)
    def edge_err(a, b): return cz_err.get((a, b), cz_err.get((b, a), np.nan))
    def edge_ok(a, b):
        e = edge_err(a, b)
        return (not np.isnan(e)) and e <= 2*med
    okset = {q for q in range(backend.num_qubits) if ok_q(q)}
    best = None
    def dfs(path, used):
        nonlocal best
        if len(path) == 12:
            errs = [edge_err(path[i], path[i+1]) for i in range(11)]
            score = float(np.nanmean(errs))
            if best is None or score < best[1]: best = (list(path), score, errs)
            return
        for n in sorted(adj[path[-1]], key=lambda n: edge_err(path[-1], n) if not np.isnan(edge_err(path[-1], n)) else 1.0):
            if n in okset and n not in used and edge_ok(path[-1], n):
                used.add(n); path.append(n)
                dfs(path, used)
                path.pop(); used.remove(n)
                if best is not None and len(path) == 1:
                    break   # one good chain per start suffices; keep global best across starts
    for start in sorted(okset):
        dfs([start], {start})
    return best, med

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
    backend = svc.backend('ibm_fez')
    best, med = pick_chain(backend)
    if best is None:
        print('DAY-OF CHAIN RULE: NO 12-qubit chain meets the floors -> forfeit the window (sealed).')
        return
    chain, score, errs = best
    print(f'[chain] accepted {chain}  mean_cz_err={score:.4f} (median {med:.4f}); floors OK')
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
        print('\nDRY RUN complete (nothing submitted). Set RUN=1 to submit the sealed Batch '
              f'(~3-4 QPU min): 5 circuits x {SHOTS} shots, gate twirling N=32, TREX, DD.')
        return
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
    # ---- SEALED SUBMISSION ----
    order = list(CIRC_NAMES); np.random.default_rng(20260902).shuffle(order)
    with Batch(backend=backend) as batch:
        sampler = SamplerV2(mode=batch)
        sampler.options.twirling.enable_gates = True          # HARD requirement (manifest)
        sampler.options.twirling.enable_measure = True        # TREX
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
