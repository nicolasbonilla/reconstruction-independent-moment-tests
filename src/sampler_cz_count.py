# -*- coding: utf-8 -*-
r"""CZ COUNT OF THE SAMPLER CIRCUITS ON ibm_fez  ->  data/sampler_cz_count.json

Two independent sources for the two-qubit (CZ) gate count of the Krylov-Trotter sampler circuits
(K = 7 circuits, k = 0..6; k = 0 has dt = 0 and therefore no entangling gate):

(A) REPRODUCTION.  `build_circuit(k)` is imported from hardware_matched_job.py (L = 6, 12 qubits,
    the companion's circuit family) and hardware_matched_job_L8.py (L = 8, 16 qubits, our device
    runs) and transpiled to the ibm_fez target (FakeFez of the installed qiskit-ibm-runtime) with
    the exact `transpile` call of each job script (optimization_level=3; the L = 8 script also passes
    translation_method='translator'), seed_transpiler = 0..4.

(B) DEPLOYED ISA CIRCUITS.  The ISA circuits actually submitted to ibm_fez, decoded from the job
    inputs (`params.pubs`) of
      - the four L = 8 Sampler jobs (2026-08-29 UTC) whose counts are in
        data/heron_counts_matched_L8_<job_id>.json, and
      - the two 12-qubit Sampler jobs of 2026-08-21 (da44t8u1vhnc73fkie30, da459juaa69c739jcb80).
    Each deployed circuit is identified by simulation: its measured-bit output distribution is
    compared (total-variation distance) with that of build_circuit(k) at the job's parameters.

    Source of the job inputs, in order of preference:
      1. env P3_JOB_JSON_DIR: a directory of REST responses <job_id>_job.json (GET /v1/jobs/<id>)
         and <job_id>_metrics.json (GET /v1/jobs/<id>/metrics);
      2. env P3_FETCH_LIVE=1: QiskitRuntimeService(name=IBM_ACCOUNT).job(id).inputs (read-only);
      3. the encoded pubs already deposited in data/sampler_cz_count.json (offline re-count).
    Only the encoded circuits (`pubs`), the backend name, the creation time and the client versions
    are stored; no account data.

Environment used for the deposited output: qiskit 2.5.2, qiskit-ibm-runtime 0.50.0, numpy 2.4.6,
scipy 1.17.1 (no qiskit-aer needed: identification uses qiskit.quantum_info.Statevector).
Run:  python sampler_cz_count.py        (about 2-4 min, dominated by the L = 8 ED at import)
"""
import os, sys, json, glob, time, hashlib, importlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.normpath(os.path.join(HERE, '..', 'data'))
OUT = os.path.join(RES, 'sampler_cz_count.json')
SEEDS = [0, 1, 2, 3, 4]
ACCOUNT = os.environ.get('IBM_ACCOUNT', 'open-instance')

L8_JOBS = ['da983qse74ec73ajfr2g', 'da98cksjbipc73fej1tg', 'da98icsjbipc73fej78g', 'da98k16rbfbs73chckvg']
L6_JOBS = ['da44t8u1vhnc73fkie30', 'da459juaa69c739jcb80']

import qiskit
import qiskit_ibm_runtime
from qiskit import transpile, QuantumCircuit
from qiskit.quantum_info import Statevector
from qiskit_ibm_runtime import RuntimeDecoder, RuntimeEncoder
from qiskit_ibm_runtime.fake_provider import FakeFez

t0 = time.time()
log = lambda *a: print(f"[{time.time()-t0:6.1f}s]", *a, flush=True)


def two_q(qc):
    ops = qc.count_ops()
    return {'cz': int(ops.get('cz', 0)),
            'n2q': int(sum(1 for ci in qc.data if ci.operation.num_qubits == 2)),
            'depth': int(qc.depth()),
            'depth_2q': int(qc.depth(lambda ci: ci.operation.num_qubits == 2))}


def stats(v):
    v = np.asarray(v, float)
    return {'n': int(v.size), 'min': float(v.min()), 'max': float(v.max()),
            'median': float(np.median(v)), 'mean': float(v.mean())}


# ---------------------------------------------------------------- circuit families
log('importing build_circuit (L=6) ...')
hm6 = importlib.import_module('hardware_matched_job')
log('importing build_circuit (L=8; runs the L=8 ED at import) ...')
hm8 = importlib.import_module('hardware_matched_job_L8')
# build_circuit only reads the scalar model parameters; release the ED arrays built at import (memory).
for _m in (hm6, hm8):
    for _a in ('C', 'Cd', 'H', 'v', 'w', 'Vn', 'En', 'coef', 'psi0', 'phi', 'Nocc', 'Szocc', 'A_exact', 'pos'):
        if hasattr(_m, _a):
            delattr(_m, _a)
import gc; gc.collect()


def logical(mod, k, U=None):
    if U is None:
        return mod.build_circuit(k)
    U0 = mod.U
    try:
        mod.U = U
        return mod.build_circuit(k)
    finally:
        mod.U = U0


# ---------------------------------------------------------------- (A) reproduction on FakeFez
fez = FakeFez()
props_path = os.path.join(os.path.dirname(qiskit_ibm_runtime.fake_provider.__file__), 'backends', 'fez', 'props_fez.json')
fez_snapshot = json.load(open(props_path)).get('last_update_date')
CALLS = {6: dict(optimization_level=3), 8: dict(optimization_level=3, translation_method='translator')}
repro = {}
for Lval, mod in ((6, hm6), (8, hm8)):
    circs = [logical(mod, k) for k in range(mod.K)]
    rows = []
    for s in SEEDS:
        isa = transpile(circs, backend=fez, seed_transpiler=s, **CALLS[Lval])
        for k, qc in enumerate(isa):
            rows.append({'seed': s, 'k': k, **two_q(qc)})
    cz_k = [r['cz'] for r in rows if r['k'] >= 1]
    repro[f'L{Lval}'] = {
        'n_qubits': 2 * Lval, 'transpile_kwargs': CALLS[Lval], 'seeds': SEEDS,
        'logical_2q_gates_per_circuit_k_ge_1': int(sum(1 for ci in circs[1].data if ci.operation.num_qubits == 2)),
        'per_circuit': rows,
        'cz_k_ge_1': stats(cz_k),
        'cz_k_ge_1_median_per_seed': {str(s): float(np.median([r['cz'] for r in rows if r['seed'] == s and r['k'] >= 1]))
                                       for s in SEEDS},
    }
    log(f"FakeFez L={Lval}: CZ (k>=1, 5 seeds) {repro[f'L{Lval}']['cz_k_ge_1']}")


# ---------------------------------------------------------------- (B) deployed ISA circuits
def load_prev():
    try:
        return json.load(open(OUT)).get('deployed', {}).get('jobs', {})
    except Exception:
        return {}


def fetch_job(jid, prev):
    d = os.environ.get('P3_JOB_JSON_DIR')
    if d and os.path.exists(os.path.join(d, f'{jid}_job.json')):
        j = json.load(open(os.path.join(d, f'{jid}_job.json')))
        mpath = os.path.join(d, f'{jid}_metrics.json')
        m = json.load(open(mpath)) if os.path.exists(mpath) else {}
        return {'source': 'REST job JSON (GET /v1/jobs/<id>)', 'backend': j['backend'], 'created': j['created'],
                'client': m.get('qiskit_version'), 'pubs_encoded': j['params']['pubs']}
    if os.environ.get('P3_FETCH_LIVE') == '1':
        from qiskit_ibm_runtime import QiskitRuntimeService
        job = QiskitRuntimeService(name=ACCOUNT).job(jid)
        pubs = json.loads(json.dumps(job.inputs['pubs'], cls=RuntimeEncoder))
        return {'source': 'QiskitRuntimeService.job(id).inputs', 'backend': job.backend().name,
                'created': str(job.creation_date), 'client': None, 'pubs_encoded': pubs}
    if jid in prev:
        tag = ' [re-read from deposited payload]'
        e = dict(prev[jid]); e['source'] = e.get('source', '').replace(tag, '') + tag
        return e
    return None


def measured_distribution(qc):
    """Output distribution over the classical bits (clbit 0 = least significant)."""
    used = sorted({qc.find_bit(q).index for ci in qc.data if ci.operation.name not in ('barrier',)
                   for q in ci.qubits})
    pos = {q: i for i, q in enumerate(used)}
    red = QuantumCircuit(len(used))
    meas = {}
    for ci in qc.data:
        name = ci.operation.name
        if name == 'barrier':
            continue
        qi = [pos[qc.find_bit(q).index] for q in ci.qubits]
        if name == 'measure':
            meas[qc.find_bit(ci.clbits[0]).index] = qi[0]
            continue
        red.append(ci.operation, qi)
    qargs = [meas[c] for c in sorted(meas)]
    return Statevector(red).probabilities(qargs), len(used)


def tvd(p, q):
    return float(0.5 * np.abs(np.asarray(p) - np.asarray(q)).sum())


_plog = {}


def logical_distribution(mod, k, U):
    key = (mod.__name__, k, U)
    if key not in _plog:
        _plog[key] = measured_distribution(logical(mod, k, U))[0]
    return _plog[key]


prev = load_prev()
deployed = {}
JOBSETS = [(8, hm8, L8_JOBS, [None]), (6, hm6, L6_JOBS, [4.0, 3.5])]
for Lval, mod, jids, Ucands in JOBSETS:
    for jid in jids:
        e = fetch_job(jid, prev)
        if e is None:
            log(f'job {jid}: no source available, skipped')
            continue
        pubs = json.loads(json.dumps(e['pubs_encoded']), cls=RuntimeDecoder)
        rows = []
        for k, pub in enumerate(pubs):
            qc = pub[0] if isinstance(pub, (list, tuple)) else pub
            shots = pub[2] if isinstance(pub, (list, tuple)) and len(pub) > 2 else None
            p_isa, nact = measured_distribution(qc)
            gc.collect()
            tv = {float(Uc if Uc is not None else mod.U): tvd(p_isa, logical_distribution(mod, k, Uc))
                  for Uc in Ucands}
            Ubest = min(tv, key=tv.get)
            # k=0 (dt=0) does not depend on U: candidates tie and U is not identified
            ident = len(tv) == 1 or sorted(tv.values())[1] - tv[Ubest] > 1e-6
            rows.append({'k': k, **two_q(qc), 'active_qubits': nact, 'pub_shots': shots,
                         'tvd_vs_build_circuit_by_U': {str(u): d for u, d in tv.items()},
                         'matches_build_circuit_U': Ubest if ident else None,
                         'tvd_vs_build_circuit': tv[Ubest]})
        e_out = {kk: vv for kk, vv in e.items() if kk != 'pubs_encoded'}
        e_out.update({'L': Lval, 'per_circuit': rows,
                      'cz_k_ge_1': stats([r['cz'] for r in rows if r['k'] >= 1]),
                      'pubs_sha256': hashlib.sha256(json.dumps(e['pubs_encoded'], sort_keys=True).encode()).hexdigest(),
                      'pubs_encoded': e['pubs_encoded']})
        e_out['source'] = e['source']
        deployed[jid] = e_out
        log(f"job {jid} (L={Lval}, {e['created']}): CZ k=0..6 {[r['cz'] for r in rows]}  "
            f"max TVD vs build_circuit {max(r['tvd_vs_build_circuit'] for r in rows):.2e}  "
            f"U(k>=1)={sorted({r['matches_build_circuit_U'] for r in rows if r['k'] >= 1})}")

pooled = {}
for Lval in (8, 6):
    cz = [r['cz'] for jid, e in deployed.items() if e['L'] == Lval for r in e['per_circuit'] if r['k'] >= 1]
    if cz:
        vals, cnts = np.unique(cz, return_counts=True)
        pooled[f'L{Lval}'] = {**stats(cz), 'jobs': [j for j, e in deployed.items() if e['L'] == Lval],
                              'values_sorted': sorted(int(x) for x in cz),
                              'mode': int(vals[np.argmax(cnts)])}
        log(f"deployed L={Lval}: CZ (k>=1) {pooled[f'L{Lval}']}")

out = {
    'description': ('Transpiled CZ count of the K=7 Krylov-Trotter sampler circuits (k=0..6) on ibm_fez. '
                    'k=0 has no entangling gate and is excluded from every summary (k>=1). '
                    "'reproduction': FakeFez transpilation of build_circuit(k) with each job script's transpile call, "
                    "seeds 0-4. 'deployed': the ISA circuits actually submitted, decoded from the job inputs and "
                    'identified against build_circuit(k) by the total-variation distance of their output distributions.'),
    'generator': 'src/sampler_cz_count.py',
    'run_date': time.strftime('%Y-%m-%d'),
    'versions': {'python': sys.version.split()[0], 'qiskit': qiskit.__version__,
                 'qiskit_ibm_runtime': qiskit_ibm_runtime.__version__, 'numpy': np.__version__},
    'fakefez_props_last_update': fez_snapshot,
    'reproduction': repro,
    'deployed': {'summary_k_ge_1': pooled, 'jobs': deployed},
    'headline': {
        'G_L8_deployed_median': pooled.get('L8', {}).get('median'),
        'G_L6_deployed_median': pooled.get('L6', {}).get('median'),
        'note': ('G_L8_deployed_median is the reference depth used by offdiag_gsurface.py and '
                 'offdiag_noise_forecast.py (median over the 24 deployed L=8 circuits with k>=1).'),
    },
}
os.makedirs(RES, exist_ok=True)
with open(OUT, 'w') as fh:
    json.dump(out, fh, indent=1)
log(f'wrote {OUT}')
