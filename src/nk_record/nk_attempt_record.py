# -*- coding: utf-8 -*-
"""Build data/2026-09-05_nk_device_attempt_record.json: the record of the sealed n_k device attempt.

What it records
  1. The sealed pre-registration (data/manifest_nk_device_v1.json, SHA-256 recomputed here) and its
     day-of chain rule, quoted verbatim from the manifest.
  2. The deposited scripts in src/nk_record/ with their SHA-256 and provenance.
  3. The logged outputs of 2026-09-03 .. 2026-09-05 (src/nk_record/nk_logged_outputs.py), verbatim,
     with the transcript UTC timestamps: two RUN=1 launches on 2026-09-03 that submitted nothing, the
     ibm_kingston window of 2026-09-04 under the proposed (never sealed) amended rule, the gate-v4
     denial of 2026-09-05 and the v1-rule checks of 2026-09-05.
  4. A re-query, today, of the historical calibrations (backend.properties(datetime=...)) at the
     logged instants, with the sealed v1 day-of chain rule applied to each snapshot: qualifying-qubit
     counts, largest connected qualifying set, longest qualifying chain, verdict. The per-qubit and
     per-edge calibration tables used are stored in the record, so `--offline` recomputes every
     number without network access.
  5. The jobs created on the account after the seal (expected: none).

Modes
  python src/nk_record/nk_attempt_record.py            live: read-only IBM queries (0 QPU), then write
  python src/nk_record/nk_attempt_record.py --offline  recompute from the tables already in the record

Credentials: QiskitRuntimeService() reads the locally saved account; no token is read, printed or
written by this script. Needs qiskit-ibm-runtime for the live mode only.
"""
import os, sys, json, hashlib, datetime as dt, statistics, collections, importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(ROOT, 'data')
OUT = os.path.join(DATA, '2026-09-05_nk_device_attempt_record.json')
MANIFEST = os.path.join(DATA, 'manifest_nk_device_v1.json')
VERDICT_V4 = os.path.join(DATA, '2026-09-05_nk_v4_verdict.json')

NQ = 12                         # chain length: L = 6 sites x 2 spins
T1_MIN, T2_MIN = 150e-6, 100e-6  # manifest v1 noise_floors
# proposed amended rule (day_of_rule_v2; denied by gate v4, never sealed) - used ONLY to check that the
# 2026-09-04 ibm_kingston snapshot is the one the watcher saw
V2_T1, V2_T2, V2_RO, V2_CZ_ABS = 100e-6, 70e-6, 0.030, 3.5e-3

SNAPSHOTS = [  # (key, backend, query instant UTC, why)
    ('fez_2026-09-03_launch1', 'ibm_fez', '2026-09-03T14:27:21Z', 'instant of RUN=1 launch 1 (logged L01)'),
    ('fez_2026-09-03_launch2', 'ibm_fez', '2026-09-03T14:28:29Z', 'instant of RUN=1 launch 2 (logged L03)'),
    ('kingston_2026-09-04_window', 'ibm_kingston', '2026-09-04T17:13:53Z',
     'watcher v3 exit on the ibm_kingston window under the proposed amended rule (logged L18)'),
    ('fez_2026-09-05', 'ibm_fez', '2026-09-05T13:08:27Z', 'instant of v1_chain_check.py output (logged L23)'),
    ('marrakesh_2026-09-05', 'ibm_marrakesh', '2026-09-05T13:08:27Z', 'instant of v1_chain_check.py output (logged L23)'),
    ('kingston_2026-09-05', 'ibm_kingston', '2026-09-05T13:08:27Z', 'instant of v1_chain_check.py output (logged L23)'),
]
KINGSTON_WINDOW_CHAIN = [13, 14, 15, 19, 35, 34, 33, 32, 31, 30, 29, 28]   # logged L18

SCRIPTS = {
    'src/nk_record/nk_stage2_device_job_launch1.py':
        'verbatim: written 2026-09-02T20:09:07Z, committed in 9bdb1d5 (pre-registration commit, '
        '2026-09-02T20:09:52Z) and run unmodified as launch 1 (2026-09-03T14:27:21Z); its traceback '
        'line numbers (48, 123, 168) match this file',
    'src/nk_record/nk_stage2_device_job_launch2.py':
        'reconstructed: launch1 plus the two Edit calls of 2026-09-03T14:28:08Z and 14:28:18Z '
        '(tool ids toolu_01HiVEcDLyNpDeibGT4MfXHq, toolu_016CwCzfZ2CSAWNKnTHgMaGC), applied verbatim; '
        'the state run as launch 2 (2026-09-03T14:28:29Z). Only change vs launch1: a qubit with an '
        'uncalibrated T1 or T2 fails the floors instead of raising',
    'src/nk_record/nk_stage2_device_job_20260903_1505Z.py':
        'reconstructed: launch2 plus the Edit calls of 2026-09-03T14:49:37Z (exhaustive DFS for the '
        'same v1 rule) and 15:05:40Z (free-tier guard); the module whose pick_chain the watchers v1/v2 '
        'imported. Never run with RUN=1. Applying the 2026-09-04T13:42:42Z patch to this file reproduces '
        'the archived work copy except one line, whose literal "\\n" (a parse error) the archived file '
        'carries; that later work copy is published, with that line repaired, as '
        'src/nk_stage2_device_job_v2_amendment.py',
    'src/nk_record/green_light_watch.py': 'verbatim (archived work copy; equals the heredoc of 2026-09-03T14:50:01Z)',
    'src/nk_record/green_light_watch_v2.py': 'verbatim (archived work copy; equals the heredoc of 2026-09-03T14:57:23Z)',
    'src/nk_record/green_light_watch_v3.py':
        'verbatim up to line endings (archived work copy; equals the heredoc of 2026-09-04T11:25:52Z '
        'after the patch of 2026-09-04T17:16:11Z)',
    'src/nk_record/v1_chain_check.py': 'verbatim (archived work copy; equals the Write of 2026-09-05T13:08:07Z)',
    'src/nk_record/v1_chain_diagnose.py': 'verbatim (archived work copy; equals the Write of 2026-09-05T13:09:00Z)',
    'src/nk_record/nk_logged_outputs.py': 'logged outputs, transcribed verbatim (see its docstring)',
    'src/nk_record/nk_attempt_record.py': 'this generator',
}


def sha256_file(p, normalize_eol=False):
    b = open(p, 'rb').read()
    if normalize_eol:
        b = b.replace(b'\r\n', b'\n')
    return hashlib.sha256(b).hexdigest()


def utc(s):
    return dt.datetime.fromisoformat(s.replace('Z', '+00:00'))


# ----------------------------------------------------------------------------- live query
def fetch_table(svc, backend, when):
    """Reduce backend.properties(datetime=when) to the quantities the rules use."""
    be = svc.backend(backend)
    props = be.properties(datetime=utc(when))
    n = be.num_qubits

    def get(fn, *a):
        try:
            v = fn(*a)
            return float(v) if v is not None and v > 0 else None
        except Exception:
            return None

    edges = sorted({tuple(sorted(e)) for e in be.coupling_map.get_edges()})
    cz = []
    for a, b in edges:
        e = get(props.gate_error, 'cz', [a, b])
        if e is None:
            e = get(props.gate_error, 'cz', [b, a])
        cz.append([a, b, e])
    directed = [list(e) for e in be.coupling_map.get_edges()]
    cz_directed = [get(props.gate_error, 'cz', e) for e in directed]
    return {
        'backend': backend, 'query_utc': when,
        'calibration_last_update': str(props.last_update_date),
        'num_qubits': n,
        't1_s': [get(props.t1, q) for q in range(n)],
        't2_s': [get(props.t2, q) for q in range(n)],
        'readout_error': [get(props.readout_error, q) for q in range(n)],
        'cz_error_undirected': cz,
        'cz_median_directed_list': (statistics.median([x for x in cz_directed if x is not None])
                                    if any(x is not None for x in cz_directed) else None),
    }


def jobs_after_seal(svc, seal_utc):
    out = []
    for j in svc.jobs(limit=200, created_after=utc(seal_utc)):
        try:
            created = str(j.creation_date)
        except Exception:
            created = None
        out.append({'job_id': j.job_id(), 'created': created, 'backend': j.backend().name})
    latest = None
    js = svc.jobs(limit=1, descending=True)
    if js:
        latest = {'job_id': js[0].job_id(), 'created': str(js[0].creation_date)}
    return out, latest


# ----------------------------------------------------------------------------- rule application
def components(nodes, adj):
    seen, sizes = set(), []
    for s in sorted(nodes):
        if s in seen:
            continue
        comp, st = {s}, [s]
        while st:
            x = st.pop()
            for y in adj.get(x, ()):
                if y not in comp:
                    comp.add(y); st.append(y)
        seen |= comp
        sizes.append(len(comp))
    return sorted(sizes, reverse=True)


def longest_path(adj, cap=NQ):
    """Exhaustive longest simple path, stopped as soon as one of length `cap` is found."""
    best = [0, None]

    def walk(path, seen):
        if len(path) > best[0]:
            best[0], best[1] = len(path), list(path)
        if best[0] >= cap:
            return True
        for v in sorted(adj.get(path[-1], ())):
            if v not in seen:
                seen.add(v); path.append(v)
                if walk(path, seen):
                    return True
                path.pop(); seen.discard(v)
        return False

    for s in sorted(adj):
        if walk([s], {s}):
            break
    return best[0], best[1]


def graph(edges, qok, eok):
    adj = collections.defaultdict(set)
    for a, b, e in edges:
        if qok(a) and qok(b) and eok(e):
            adj[a].add(b); adj[b].add(a)
    return adj


def apply_v1(tab):
    n = tab['num_qubits']
    t1, t2 = tab['t1_s'], tab['t2_s']
    edges = tab['cz_error_undirected']
    errs = [e for _, _, e in edges if e is not None]
    med = statistics.median(errs)
    cap = 2.0 * med
    qual = {q for q in range(n) if t1[q] is not None and t2[q] is not None
            and t1[q] >= T1_MIN and t2[q] >= T2_MIN}
    e_ok = lambda e: e is not None and e <= cap
    g_floor = graph(edges, lambda q: q in qual, lambda e: True)        # T1/T2 floors only
    g_cz = graph(edges, lambda q: True, e_ok)                          # CZ homogeneity only
    g_v1 = graph(edges, lambda q: q in qual, e_ok)                     # the full v1 rule
    L_floor, _ = longest_path(g_floor)
    L_cz, _ = longest_path(g_cz)
    L_v1, path_v1 = longest_path(g_v1)
    t1v = [x for x in t1 if x is not None]; t2v = [x for x in t2 if x is not None]
    return {
        'calibration_last_update': tab['calibration_last_update'],
        'n_qubits': n,
        'n_T1_calibrated': len(t1v), 'n_T2_calibrated': len(t2v),
        'T1_median_us': round(statistics.median(t1v) * 1e6, 1),
        'T2_median_us': round(statistics.median(t2v) * 1e6, 1),
        'n_T1_ge_150us': sum(1 for x in t1v if x >= T1_MIN),
        'n_T2_ge_100us': sum(1 for x in t2v if x >= T2_MIN),
        'n_qualifying_qubits': len(qual),
        'n_edges': len(edges), 'n_edges_with_cz_error': len(errs),
        'cz_median': med, 'cz_cap_2x_median': cap,
        'cz_median_directed_list': tab.get('cz_median_directed_list'),
        'n_edges_within_cap': sum(1 for e in errs if e <= cap),
        'largest_connected_qualifying_set_coupling_graph': components(qual, g_floor)[:1][0] if qual else 0,
        'component_sizes_coupling_graph_top6': components(qual, g_floor)[:6],
        'largest_connected_qualifying_set_v1_graph': (components(qual, g_v1)[:1] or [0])[0],
        'component_sizes_v1_graph_top6': components(qual, g_v1)[:6],
        'longest_chain_floors_only': L_floor,
        'longest_chain_cz_cap_only': L_cz,
        'longest_chain_v1_rule': L_v1,
        'longest_chain_v1_rule_qubits': path_v1,
        'v1_chain_of_12_exists': L_v1 >= NQ,
        'verdict': ('a 12-qubit chain satisfies the sealed v1 rule' if L_v1 >= NQ else
                    'FORFEIT: no 12-qubit chain satisfies the sealed v1 rule'),
        'binding_constraint': (None if L_v1 >= NQ else
                               'T1/T2 floors' if L_floor < NQ <= L_cz else
                               'CZ homogeneity' if L_cz < NQ <= L_floor else
                               'both' if max(L_floor, L_cz) < NQ else 'their combination'),
    }


def check_chain(tab, chain):
    """Per-qubit and per-edge values of one chain under a snapshot (used for the kingston window)."""
    t1, t2, ro = tab['t1_s'], tab['t2_s'], tab['readout_error']
    ed = {(min(a, b), max(a, b)): e for a, b, e in tab['cz_error_undirected']}
    errs = [e for e in ed.values() if e is not None]
    med = statistics.median(errs)
    ce = [ed.get((min(chain[i], chain[i + 1]), max(chain[i], chain[i + 1]))) for i in range(len(chain) - 1)]
    is_path = all(e is not None for e in ce)
    r = {
        'chain': chain, 'is_coupled_path': is_path,
        'T1_min_us': round(min(t1[q] for q in chain) * 1e6, 1),
        'T2_min_us': round(min(t2[q] for q in chain) * 1e6, 1),
        'readout_max': max(ro[q] for q in chain),
        'cz_max_on_chain': max(ce) if is_path else None,
        'n_chain_qubits_failing_v1_floors': sum(1 for q in chain if not (t1[q] >= T1_MIN and t2[q] >= T2_MIN)),
    }
    cap2 = min(V2_CZ_ABS, 2 * med)
    r['satisfies_proposed_amended_rule'] = bool(
        is_path and all(t1[q] >= V2_T1 and t2[q] >= V2_T2 and ro[q] <= V2_RO for q in chain)
        and all(e <= cap2 for e in ce) and med <= V2_CZ_ABS)
    r['satisfies_sealed_v1_rule'] = bool(
        is_path and r['n_chain_qubits_failing_v1_floors'] == 0 and all(e <= 2 * med for e in ce))
    return r


# ----------------------------------------------------------------------------- logged comparison
LOGGED_EXPECT = {  # numbers as printed in the logged outputs (key -> (logged id, value))
    'fez_2026-09-03_launch2': {
        'calibration_last_update': ('L07', '2026-09-03 15:43:34+02:00'),
        'n_qualifying_qubits': ('L05', 31),
        'component_sizes_coupling_graph_top6': ('L05', [9, 4, 2, 1, 1, 1]),
        'T1_median_us_rounded': ('L05', 129), 'T2_median_us_rounded': ('L05', 96),
        'v1_chain_of_12_exists': ('L04', False),
    },
    'fez_2026-09-05': {'cz_median_e3': ('L23', 2.779), 'n_qualifying_qubits': ('L23', 18),
                       'n_edges_within_cap': ('L24', 142), 'longest_chain_floors_only': ('L24', 3),
                       'longest_chain_cz_cap_only': ('L24', 12), 'longest_chain_v1_rule': ('L24', 2)},
    'marrakesh_2026-09-05': {'cz_median_e3': ('L23', 3.409), 'n_qualifying_qubits': ('L23', 39),
                             'n_edges_within_cap': ('L24', 146), 'longest_chain_floors_only': ('L24', 4),
                             'longest_chain_cz_cap_only': ('L24', 12), 'longest_chain_v1_rule': ('L24', 4)},
    'kingston_2026-09-05': {'cz_median_e3': ('L23', 2.040), 'n_qualifying_qubits': ('L23', 87),
                            'n_edges_within_cap': ('L24', 145), 'longest_chain_floors_only': ('L24', 9),
                            'longest_chain_cz_cap_only': ('L24', 12), 'longest_chain_v1_rule': ('L24', 7)},
}


def compare(key, res):
    exp = LOGGED_EXPECT.get(key)
    if not exp:
        return None
    view = dict(res)
    view['cz_median_e3'] = round(res['cz_median'] * 1e3, 3)
    view['T1_median_us_rounded'] = round(res['T1_median_us'])
    view['T2_median_us_rounded'] = round(res['T2_median_us'])
    out = {}
    for k, (lid, v) in exp.items():
        out[k] = {'logged': v, 'logged_in': lid, 'requeried': view[k], 'match': view[k] == v}
    out['all_match'] = all(x['match'] for x in out.values() if isinstance(x, dict))
    return out


# ----------------------------------------------------------------------------- replay of the launched code
class _Props:
    """Stand-in for BackendProperties built from a stored table: raises where IBM's object raises."""
    def __init__(self, tab):
        self.t = tab
        self.cz = {}
        for a, b, e in tab['cz_error_undirected']:
            if e is not None:
                self.cz[(a, b)] = self.cz[(b, a)] = e

    def _q(self, key, q):
        v = self.t[key][q]
        if v is None:
            raise KeyError(f"Couldn't find the property for qubit {q}.")
        return v

    def t1(self, q): return self._q('t1_s', q)
    def t2(self, q): return self._q('t2_s', q)
    def readout_error(self, q): return self._q('readout_error', q)

    def gate_error(self, gate, qubits):
        return self.cz[tuple(qubits)]


class _Backend:
    def __init__(self, tab):
        self._p = _Props(tab)
        self.num_qubits = tab['num_qubits']
        self.coupling_map = [[a, b] for a, b, _ in tab['cz_error_undirected']] + \
                            [[b, a] for a, b, _ in tab['cz_error_undirected']]

    def properties(self, **kw): return self._p


def replay_pick_chain(script, tab):
    """Run the pick_chain function of a deposited launch script, unmodified, on a stored table."""
    import ast, numpy as np
    src = open(os.path.join(HERE, script), encoding='utf-8').read()
    fn = next(n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef) and n.name == 'pick_chain')
    ns = {'np': np}
    exec(compile(ast.Module(body=[fn], type_ignores=[]), script, 'exec'), ns)
    try:
        best, med = ns['pick_chain'](_Backend(tab))
        return {'raised': None, 'chain_found': best is not None, 'cz_median': med,
                'printed_by_main': ('DAY-OF CHAIN RULE: NO 12-qubit chain meets the floors -> forfeit the window (sealed).'
                                    if best is None else None)}
    except Exception as e:
        return {'raised': f'{type(e).__name__}: {e}', 'chain_found': None}


# ----------------------------------------------------------------------------- main
def main():
    offline = '--offline' in sys.argv
    spec = importlib.util.spec_from_file_location('nk_logged_outputs', os.path.join(HERE, 'nk_logged_outputs.py'))
    lo = importlib.util.module_from_spec(spec); spec.loader.exec_module(lo)
    man = json.load(open(MANIFEST, encoding='utf-8'))
    man_sha = sha256_file(MANIFEST)
    man_sealed = open(MANIFEST.replace('.json', '.sha256')).read().split()[0]
    v4 = json.load(open(VERDICT_V4, encoding='utf-8'))

    if offline:
        old = json.load(open(OUT, encoding='utf-8'))
        tables = old['calibration_tables']
        jobs, latest = old['jobs_after_seal']['jobs'], old['jobs_after_seal']['most_recent_job_on_account']
        queried_utc = old['requery']['queried_utc']
    else:
        import warnings; warnings.filterwarnings('ignore')
        from qiskit_ibm_runtime import QiskitRuntimeService
        svc = QiskitRuntimeService()
        queried_utc = dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
        tables = {}
        for key, bname, when, _ in SNAPSHOTS:
            tables[key] = fetch_table(svc, bname, when)
            print(f'[query] {key}: {bname} @ {when} -> calibration {tables[key]["calibration_last_update"]}')
        jobs, latest = jobs_after_seal(svc, man['sealed_utc'])

    snaps = {}
    for key, bname, when, why in SNAPSHOTS:
        res = apply_v1(tables[key])
        res.update({'backend': bname, 'query_utc': when, 'why': why})
        cmp_ = compare(key, res)
        if cmp_ is not None:
            res['comparison_with_logged'] = cmp_
        snaps[key] = res
        print(f'[v1] {key:28s} cal {res["calibration_last_update"]}  qualifying {res["n_qualifying_qubits"]:3d}/{res["n_qubits"]}'
              f'  largest set {res["largest_connected_qualifying_set_v1_graph"]:3d}  longest v1 chain {res["longest_chain_v1_rule"]:2d}'
              f'  -> {res["verdict"]}' + ('' if cmp_ is None else f'  [logged match: {cmp_["all_match"]}]'))
    kw = check_chain(tables['kingston_2026-09-04_window'], KINGSTON_WINDOW_CHAIN)
    print('[kingston window chain]', json.dumps(kw))
    t1_missing = [q for q, v in enumerate(tables['fez_2026-09-03_launch1']['t1_s']) if v is None]
    t2_missing = [q for q, v in enumerate(tables['fez_2026-09-03_launch1']['t2_s']) if v is None]
    replay = {
        'note': ('pick_chain of each deposited launch script, executed unmodified on the re-queried table of '
                 'its launch instant (a stand-in properties object raises on a missing value, as IBM\'s does)'),
        'fez_qubits_without_T1': t1_missing, 'fez_qubits_without_T2': t2_missing,
        'launch1': replay_pick_chain('nk_stage2_device_job_launch1.py', tables['fez_2026-09-03_launch1']),
        'launch2': replay_pick_chain('nk_stage2_device_job_launch2.py', tables['fez_2026-09-03_launch2']),
        'watcher_module_20260903_1505Z': replay_pick_chain('nk_stage2_device_job_20260903_1505Z.py',
                                                           tables['fez_2026-09-03_launch2']),
    }
    print('[replay]', json.dumps(replay))

    logged = {r['id']: r for r in lo.LOGGED}
    rec = {
        'title': 'n_k on-device discriminating falsifier: record of the sealed device attempt (2026-09-03 to 2026-09-05)',
        'license': 'CC BY 4.0',
        'generator': 'src/nk_record/nk_attempt_record.py (live re-query; --offline recomputes from calibration_tables)',
        'summary': {
            'jobs_submitted': 0, 'shots': 0, 'qpu_seconds': 0,
            'statement': ('The sealed job was launched twice on 2026-09-03 against live ibm_fez. Launch 1 '
                          '(14:27:21Z) stopped in the chain-selection step on an uncalibrated T2 (qubit 72); '
                          'launch 2 (14:28:29Z) applied the sealed day-of chain rule and forfeited the window. '
                          'Neither reached the submission code. On 2026-09-04 a window opened on ibm_kingston '
                          'only under a proposed amended rule (looser T1/T2 floors plus a readout cap) that '
                          'required a manifest v2; it could not fire. On 2026-09-05 gate v4 denied the '
                          'amendment, so manifest v1 stands, and no backend offered a v1-compliant chain. No '
                          'device job was created after the seal.'),
            'requery_verdict': ('every re-queried snapshot (ibm_fez 2026-09-03 at both launches; ibm_fez, '
                                'ibm_marrakesh, ibm_kingston 2026-09-05; ibm_kingston 2026-09-04 window) '
                                'forfeits under the sealed v1 rule' if not any(s['v1_chain_of_12_exists'] for s in snaps.values())
                                else 'at least one re-queried snapshot admits a v1-compliant chain: see requery'),
        },
        'manifest_v1': {
            'path': 'data/manifest_nk_device_v1.json', 'sha256_recomputed': man_sha,
            'sha256_sealed': man_sealed, 'sha256_match': man_sha == man_sealed,
            'sealed_utc': man['sealed_utc'],
            'preregistration_commit': {'hash': '9bdb1d5', 'committed': '2026-09-02T20:09:52Z'},
            'day_of_chain_rule_verbatim': man['noise_floors']['day_of_chain_rule'],
            'kill_rule_verbatim': man['criteria']['kill_rules'][1],
            'chain_length_qubits': NQ,
            'implementation_note': ('the rule is applied as: every chain qubit has calibrated T1 >= 150 us and '
                                    'T2 >= 100 us (a qubit without a calibrated T1 or T2 fails); every one of the 11 '
                                    'chain edges has CZ error <= 2 x the median CZ error over the coupled pairs of the '
                                    'same snapshot; the chain is a simple path in the coupling map'),
        },
        'events': [
            {'utc': '2026-09-02T20:08:16Z', 'event': 'manifest v1 sealed (SHA-256 496203db...)', 'source': 'manifest sealed_utc'},
            {'utc': logged['L01']['utc'], 'event': 'launch 1, RUN=1, ibm_fez, script nk_stage2_device_job_launch1.py',
             'outcome': 'BackendPropertyError: T2 not calibrated for qubit 72, raised in pick_chain before any submission; 0 jobs',
             'logged': ['L01', 'L02']},
            {'utc': logged['L03']['utc'], 'event': 'launch 2, RUN=1, ibm_fez, script nk_stage2_device_job_launch2.py',
             'outcome': 'chain rule printed "DAY-OF CHAIN RULE: NO 12-qubit chain meets the floors -> forfeit the window (sealed)."; 0 jobs',
             'logged': ['L03', 'L04']},
            {'utc': logged['L05']['utc'], 'event': 'exhaustive check of the same calibration',
             'outcome': '31/156 qubits pass both floors; largest connected set 9; no 12-path', 'logged': ['L05', 'L07', 'L08']},
            {'utc': logged['L06']['utc'], 'event': 'QPU usage at the time',
             'outcome': 'open (free) plan, 447 s used, 153 s remaining of 600 s', 'logged': ['L06']},
            {'utc': logged['L09']['utc'] + ' .. ' + logged['L14']['utc'],
             'event': 'watchers v1 and v2 (v1 rule, ibm_fez only; observe, never submit)',
             'outcome': 'no GREEN line; watcher v2 flushed no output before it was stopped', 'logged': ['L09', 'L10', 'L11', 'L14']},
            {'utc': logged['L12']['utc'], 'event': 'v1-rule checks the next morning',
             'outcome': 'ibm_fez 31 qualifying, chain none; at 150/100 us no backend has a 12-chain (fez 31, marrakesh 40, kingston 87 qualifying)',
             'logged': ['L12', 'L13']},
            {'utc': logged['L18']['utc'], 'event': '2026-09-04 ibm_kingston window',
             'outcome': ('watcher v3, testing the PROPOSED amended rule (T1 >= 100 us, T2 >= 70 us, readout <= 0.030, '
                         'CZ <= min(3.5e-3, 2x median); never sealed), found chain [13,14,15,19,35,34,33,32,31,30,29,28] '
                         'at 19:13 CEST; it could not fire because manifest v2 did not exist. That chain has T1 min 131 us '
                         'and T2 min 71 us, so it fails the sealed v1 floors. Further kingston windows of the same kind '
                         'were logged overnight (chain [0,1,2,3,16,23,22,21,36,41,42,43])'),
             'logged': ['L16', 'L17', 'L18', 'L19', 'L20', 'L21']},
            {'utc': logged['L22']['utc'], 'event': 'gate v4 denies the amendment',
             'outcome': ('AMENDMENT-DENIED (manifest v1 stands; gate-negative): readout ladder 0.030 -> 0.025 -> 0.018 '
                         'exhausted on R1 (0.018 row: mean 0.02804, ucb 0.02867 > threshold 0.0266667)'),
             'record': 'data/2026-09-05_nk_v4_verdict.json', 'record_verdict': v4['verdict'],
             'record_sha256': sha256_file(VERDICT_V4), 'logged': ['L22']},
            {'utc': logged['L23']['utc'], 'event': 'v1_chain_check.py on the three backends',
             'outcome': 'no v1-compliant 12-chain (qualifying qubits: fez 18, marrakesh 39, kingston 87)', 'logged': ['L23']},
            {'utc': logged['L24']['utc'], 'event': 'v1_chain_diagnose.py',
             'outcome': 'longest v1-compliant chain: kingston 7, fez 2, marrakesh 4; binding constraint T1/T2 on all three',
             'logged': ['L24']},
        ],
        'requery': {
            'queried_utc': queried_utc,
            'method': ('QiskitRuntimeService().backend(name).properties(datetime=query_utc): the calibration '
                       'in force at the logged instant; read-only, 0 QPU'),
            'snapshots': snaps,
            'kingston_window_chain_under_window_snapshot': kw,
            'replay_of_launched_chain_selection': replay,
        },
        'jobs_after_seal': {
            'created_after': man['sealed_utc'], 'count': len(jobs), 'jobs': jobs,
            'most_recent_job_on_account': latest,
        },
        'scripts': {p: {'sha256_lf': sha256_file(os.path.join(ROOT, p), normalize_eol=True), 'provenance': s}
                    for p, s in SCRIPTS.items() if p != 'src/nk_record/nk_attempt_record.py'},
        'scripts_not_redeposited': {
            'src/nk_stage2_device_job_v2_amendment.py': 'the 2026-09-04 work copy (proposed amended rule), already in src/',
            'archived_work_copy_sha256': {
                'nk_stage2_device_job.py (2026-09-04T13:42Z, does not parse: literal "\\n" at line 170)':
                    '2794d092ecb9012062b64f89e5ddc2b51f130f1e847de1c8c520e87512c38bd3'},
        },
        'logged_outputs': {
            'source': 'src/nk_record/nk_logged_outputs.py (verbatim from session transcript ' + lo.SESSION_ID + ')',
            'entries': lo.LOGGED,
        },
        'calibration_tables': tables,
    }
    json.dump(rec, open(OUT, 'w', encoding='utf-8', newline='\n'), indent=1, ensure_ascii=False)
    print('wrote', os.path.relpath(OUT, ROOT), '| jobs after seal:', len(jobs))


if __name__ == '__main__':
    main()
