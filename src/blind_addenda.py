# -*- coding: utf-8 -*-
"""Addenda to the sealed blinded test (P3 revision, 2026-09-27; plan items R2, B1, M4, M8, m7).

Reads ONLY the sealed record and the frozen pipeline; it never writes to any of them:
    data/prereg.json (+ prereg.sha256)          the seal: battery, thresholds, catalogue, seed
    data/blind_instances_public.json            m_hat, m_bar, var_loc (what the classifier saw)
    data/blind_labels_sealed.json               the labels (unsealed after scoring)
    data/blind_verdicts.json                    the frozen verdicts
    data/2026-08-24_blind_harness_score.json    the sealed score
    src/blind_classify.py                       the frozen decision rule (run VERBATIM, see below)
    src/blind_generate.py, src/spectral_lanczos.py   the sealed generator's physics (part c only)

Output: data/2026-09-27_blind_addenda.json (a NEW file; no pre-existing file is touched).

(0) Frozen rule. The decision rule lives inline in blind_classify.main(). To use it verbatim,
    frozen_classify() copies the sealed prereg into a temporary directory, writes the instances
    to be scored there as a public file, points blind_classify.OUT at that directory and runs
    blind_classify.main() unchanged. On the 300 sealed instances this reproduces the sealed
    verdicts exactly. The Monte-Carlo in part (b) uses a vectorised copy of the rule, which
    is checked against the verbatim classifier on the first redraws of each scenario.
(a) Primary endpoint (TPR/FPR, Wilson 95%), per-class and per-subclass rates with Wilson
    intervals, and the two sealed secondary predictions that failed.
(b) FPR sensitivity. The simulated estimator of the sealed test is unbiased, while the frozen
    threshold budgets a 2% bias. m_hat is redrawn for the 112 clean instances, with no bias
    and with a realised +2% bias on each moment (a -2% bias is reported as a supplement),
    under the frozen threshold.
(c) Shared-state rescoring of the 56 truncations. The sealed truncation instances are first
    REBUILT EXACTLY. The truncation keeps the top-d configurations by |psi0|^2, and |psi0|^2
    has 24-, 12- and 6-fold degenerate shells (lattice symmetry), so which members of a
    boundary shell are kept is decided by floating-point noise in psi0. The sealed run
    computed psi0 with ARPACK's default random start (the fixed v0 in spectral_lanczos._sub_gs
    came later); ARPACK draws that start from an internal generator that is seeded once per
    process. Calling eigsh without v0 for U = 3, 4, 5, 6 (the prereg U_grid order) as the
    FIRST such calls of a fresh process reproduces the sealed tie order. Every rebuilt instance
    is checked against the sealed m_bar and var_loc, and against the sealed m_hat by
    replaying the generator's RNG stream. The recovered tie order is stored in the output,
    so it does not depend on ARPACK in future. The truncations are then rescored with the
    shared-state estimator: full operators evaluated on the SAME truncated ground state g
    that produced the reconstruction, m_hat_k = <J g|(H - e0)^k|J g> + (the sealed shot-noise
    draw), e0 = <g|H|g>, scored by the frozen rule. On a device these would be ground-state
    expectation values of the prepared state; here they are computed classically from g, and
    no computational-basis samples are involved (this is not a same-sample estimate).
(d) Under-converged Krylov: m_bar vs m_exact per moment, and residual/threshold ratios. nl=1
    is exact in m0 and m1 by construction and misses m2.
(e) Calibrated scope of the sealed test: the retained ground-state support and weight in the
    truncations, the spurious-atom weights and frequencies, the U/t values and the probe.

Run from src/:  python blind_addenda.py      (about 10-20 s on one core)
"""
import contextlib
import hashlib
import io
import json
import math
import os
import platform
import shutil
import sys
import tempfile
import time
import warnings
from collections import defaultdict

import numpy as np
import scipy
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
DATA = os.path.normpath(os.path.join(HERE, '..', 'data'))
OUT_JSON = os.path.join(DATA, '2026-09-27_blind_addenda.json')

SEED = 20260927                 # all new randomness in this script derives from this seed
R_FPR = 20000                   # redraws of the 112 clean instances per scenario, part (b)
R_FPR_VERBATIM = 40             # redraws per scenario also run through the verbatim classifier
K_SHARED = 4000                 # fresh-noise redraws per truncation instance, part (c)

SEALED_INPUTS = ['prereg.json', 'prereg.sha256', 'blind_instances_public.json',
                 'blind_labels_sealed.json', 'blind_verdicts.json',
                 '2026-08-24_blind_harness_score.json']
CODE_INPUTS = ['blind_classify.py', 'blind_generate.py', 'spectral_lanczos.py',
               'blind_addenda.py']


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------
def sha256(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


def wilson(k, n, z=1.96):
    """Wilson score interval; identical to blind_score.wilson."""
    if n == 0:
        return [float('nan'), float('nan')]
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [max(0.0, c - h), min(1.0, c + h)]


def rate(k, n):
    return {'k': int(k), 'n': int(n), 'rate': (k / n) if n else float('nan'),
            'wilson95': wilson(k, n)}


def load_sealed():
    blob = open(os.path.join(DATA, 'prereg.json'), 'rb').read()
    want = open(os.path.join(DATA, 'prereg.sha256')).read().strip()
    if hashlib.sha256(blob).hexdigest() != want:
        raise SystemExit('SEAL BROKEN: prereg.json does not match prereg.sha256')
    pr = json.loads(blob)
    pub = json.load(open(os.path.join(DATA, 'blind_instances_public.json')))
    ver = json.load(open(os.path.join(DATA, 'blind_verdicts.json')))
    lab = json.load(open(os.path.join(DATA, 'blind_labels_sealed.json')))['labels']
    score = json.load(open(os.path.join(DATA, '2026-08-24_blind_harness_score.json')))
    for name, sha in (('public', pub['_prereg_sha256']), ('verdicts', ver['_prereg_sha256']),
                      ('score', score['_provenance']['prereg_sha256'])):
        if sha != want:
            raise SystemExit(f'{name} file was produced under a different prereg')
    P = {x['id']: x for x in pub['instances']}
    V = {x['id']: x for x in ver['verdicts']}
    Lb = {x['id']: x for x in lab}
    assert sorted(P) == sorted(V) == sorted(Lb) == list(range(pr['n_instances']))
    return pr, want, P, V, Lb, score


def subclass(lab):
    m = lab['mode']
    if m == 'krylov':
        return f"krylov nl={lab['params']['nl']}"
    if m == 'ac':
        return 'ac preserve_m0' if lab['params'].get('preserve_m0') else 'ac raw'
    return m


# ---------------------------------------------------------------------------------------------
# (0) the frozen rule, verbatim and vectorised
# ---------------------------------------------------------------------------------------------
def frozen_classify(instances, prereg_sha):
    """Run blind_classify.main() UNCHANGED on `instances` (list of dicts with id, m_hat, m_bar,
    var_loc) inside a temporary copy of the sealed prereg; return {id: verdict dict}."""
    import blind_classify as bc
    old_out = bc.OUT
    with tempfile.TemporaryDirectory(prefix='p3_frozen_') as tmp:
        for f in ('prereg.json', 'prereg.sha256'):
            shutil.copyfile(os.path.join(DATA, f), os.path.join(tmp, f))
        json.dump({'_prereg_sha256': prereg_sha, 'instances': instances},
                  open(os.path.join(tmp, 'blind_instances_public.json'), 'w'))
        bc.OUT = tmp
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                bc.main()
            res = json.load(open(os.path.join(tmp, 'blind_verdicts.json')))['verdicts']
        finally:
            bc.OUT = old_out
    return {v['id']: v for v in res}


def rule_vectorised(m_hat, m_bar, var, battery):
    """Vectorised copy of blind_classify.py:23-30 (validated against frozen_classify).
    m_hat, m_bar, var: arrays (..., 3). Returns reject, moment_fire, hankel_fire, fires(...,3)."""
    Ns, z, bias, tol = battery['Ns'], battery['z'], battery['assumed_bias_frac'], battery['hankel_tol']
    delta = z * np.sqrt(var / Ns + (bias * np.abs(m_hat)) ** 2)
    fires = np.abs(m_hat - m_bar) > delta
    moment_fire = fires.any(axis=-1)
    a, b, c = m_bar[..., 0], m_bar[..., 1], m_bar[..., 2]
    lam_min = 0.5 * (a + c) - np.sqrt(0.25 * (a - c) ** 2 + b ** 2)    # min eig of [[a,b],[b,c]]
    hankel_fire = (lam_min < -tol) | (b < -tol)
    return moment_fire | hankel_fire, moment_fire, hankel_fire, fires


def part0_rule_check(pr, sha, P, V):
    insts = [{'id': i, 'm_hat': P[i]['m_hat'], 'm_bar': P[i]['m_bar'], 'var_loc': P[i]['var_loc']}
             for i in sorted(P)]
    Vv = frozen_classify(insts, sha)
    same = sum(Vv[i]['reject'] == V[i]['reject'] and Vv[i]['moment_fire'] == V[i]['moment_fire']
               and Vv[i]['hankel_fire'] == V[i]['hankel_fire'] for i in V)
    dg = max(abs(a - b) for i in V for a, b in zip(Vv[i]['g'], V[i]['g']))
    dd = max(abs(a - b) for i in V for a, b in zip(Vv[i]['delta'], V[i]['delta']))
    arr = lambda key: np.array([P[i][key] for i in sorted(P)])
    rej, _, hk, _ = rule_vectorised(arr('m_hat'), arr('m_bar'), arr('var_loc'), pr['battery'])
    vec_same = int(sum(bool(rej[i]) == V[i]['reject'] for i in sorted(V)))
    # Hankel branch: min eigenvalue of H1(m_bar) over the catalogue (M2: can it ever fire?)
    mb = arr('m_bar')
    lam = 0.5 * (mb[:, 0] + mb[:, 2]) - np.sqrt(0.25 * (mb[:, 0] - mb[:, 2]) ** 2 + mb[:, 1] ** 2)
    return {
        'verbatim_classifier_reproduces_sealed_verdicts': f'{same}/{len(V)}',
        'max_abs_diff_g': dg, 'max_abs_diff_delta': dd,
        'vectorised_rule_reproduces_sealed_verdicts': f'{vec_same}/{len(V)}',
        'hankel_branch': {
            'fired_in_sealed_record': int(sum(v['hankel_fire'] for v in V.values())),
            'min_eig_H1_of_m_bar_over_300': float(lam.min()),
            'min_m1_bar_over_300': float(mb[:, 1].min()),
            'note': 'the rule evaluates H1 at the reconstruction moments m_bar only '
                    '(blind_classify.py:28-29); every catalogue reconstruction is a positive '
                    'measure, so H1(m_bar) is PSD and the branch cannot fire by construction.'},
    }


# ---------------------------------------------------------------------------------------------
# (a) primary endpoint and per-class rates
# ---------------------------------------------------------------------------------------------
def part_a(pr, V, Lb, score):
    tp = fn = fp = tn = 0
    cls = defaultdict(lambda: [0, 0])
    sub = defaultdict(lambda: [0, 0])
    fire_by = defaultdict(lambda: np.zeros(3, int))
    for i, lab in Lb.items():
        rej = V[i]['reject']
        c = lab['corrupted']
        tp += c and rej
        fn += c and not rej
        fp += (not c) and rej
        tn += (not c) and not rej
        cls[lab['mode']][0] += rej
        cls[lab['mode']][1] += 1
        sub[subclass(lab)][0] += rej
        sub[subclass(lab)][1] += 1
        fire_by[lab['mode']] += (np.array(V[i]['g']) > np.array(V[i]['delta'])).astype(int)
    ncorr, nclean = tp + fn, fp + tn
    tpr, fpr = rate(tp, ncorr), rate(fp, nclean)
    agree = (tp, fn, fp, tn) == tuple(score['confusion'][k] for k in ('TP', 'FN', 'FP', 'TN'))
    nk = [i for i, l in Lb.items() if l['mode'] != 'krylov' and l['corrupted']]
    tp_nk = sum(V[i]['reject'] for i in nk)
    return {
        'source': 'recomputed from blind_verdicts.json + blind_labels_sealed.json; must equal '
                  '2026-08-24_blind_harness_score.json',
        'agrees_with_sealed_score_file': bool(agree),
        'prereg_decision_text': pr['decision_prereg'],
        'confusion': {'TP': tp, 'FN': fn, 'FP': fp, 'TN': tn},
        'primary_endpoint': {'TPR': tpr, 'FPR': fpr,
                             'nominal_joint_FPR_target_1_minus_0.95^3': 1 - 0.95 ** 3},
        'per_class': {m: rate(*cls[m]) for m in sorted(cls)},
        'per_subclass': {s: rate(*sub[s]) for s in sorted(sub)},
        'per_class_moment_fire_counts_[m0,m1,m2]': {m: fire_by[m].tolist() for m in sorted(fire_by)},
        'post_hoc_TPR_excluding_krylov': dict(rate(tp_nk, len(nk)),
                                              label='post hoc decomposition, not a sealed endpoint'),
        'sealed_secondary_predictions_that_failed': [
            {'prediction': "per-mode catch rate 'high for ... nl=1 krylov'",
             'outcome': rate(*sub['krylov nl=1']), 'verdict': 'FAILED (caught 0/15)'},
            {'prediction': "per-mode catch rate 'LOW (by design) for ... m0-preserving ac'",
             'outcome': rate(*sub['ac preserve_m0']),
             'verdict': 'FAILED (caught 20/34 = 59%; an outcome the pre-registration did not '
                        'anticipate, not a designed capability)'},
        ],
        'sealed_secondary_predictions_that_held': {
            'trunc (expected high)': rate(*sub['trunc']),
            'ac raw (expected high)': rate(*sub['ac raw']),
            'krylov nl>=2 (expected low)': rate(sum(sub[f'krylov nl={n}'][0] for n in (2, 3, 4)),
                                                sum(sub[f'krylov nl={n}'][1] for n in (2, 3, 4))),
        },
    }


# ---------------------------------------------------------------------------------------------
# (b) FPR sensitivity of the 112 clean instances
# ---------------------------------------------------------------------------------------------
def part_b(pr, sha, P, V, Lb):
    b = pr['battery']
    ids = sorted(i for i, l in Lb.items() if l['mode'] == 'clean')
    m_ex = np.array([Lb[i]['m_exact'] for i in ids])
    m_bar = np.array([P[i]['m_bar'] for i in ids])
    var = np.array([P[i]['var_loc'] for i in ids])
    Us = np.array([Lb[i]['params']['U'] for i in ids])
    se = np.sqrt(var / b['Ns'])
    rng = np.random.default_rng([SEED, 2])
    Z = rng.standard_normal((R_FPR, len(ids), 3))            # common random numbers
    scen = {'unbiased (as in the sealed test)': 0.0, 'realized +2% bias on m0,m1,m2': +0.02,
            'supplementary: realized -2% bias on m0,m1,m2': -0.02}
    out = {'n_clean': len(ids), 'redraws_per_scenario': R_FPR, 'rng': f'default_rng([{SEED}, 2])',
           'm_bar_equals_m_exact_max_abs_diff': float(np.abs(m_bar - m_ex).max()),
           'sealed_realized': rate(sum(V[i]['reject'] for i in ids), len(ids)),
           'estimator_model': 'm_hat_k = (1 + b) m_exact_k + se_k Z, se_k = sqrt(var_loc_k/Ns); '
                              'frozen threshold delta_k = z sqrt(var_loc_k/Ns + (0.02 |m_hat_k|)^2)',
           'scenarios': {}}
    verb_checked = verb_same = 0
    for name, bias in scen.items():
        m_hat = (1.0 + bias) * m_ex[None] + se[None] * Z
        rej, mf, hk, fires = rule_vectorised(m_hat, m_bar[None], var[None], b)
        nfp = rej.sum(axis=1)                                   # FP count per redraw (of 112)
        fpr = rej.mean()
        # verbatim frozen classifier on the first R_FPR_VERBATIM redraws
        insts = []
        for r in range(R_FPR_VERBATIM):
            for j, i in enumerate(ids):
                insts.append({'id': r * 1000 + i, 'm_hat': m_hat[r, j].tolist(),
                              'm_bar': m_bar[j].tolist(), 'var_loc': var[j].tolist()})
        Vv = frozen_classify(insts, sha)
        for r in range(R_FPR_VERBATIM):
            for j, i in enumerate(ids):
                verb_checked += 1
                verb_same += Vv[r * 1000 + i]['reject'] == bool(rej[r, j])
        per_U = {}
        for U in sorted(set(Us.tolist())):
            sel = Us == U
            per_U[f'U={U:g}'] = {'n_instances': int(sel.sum()),
                                 'per_instance_FPR': float(rej[:, sel].mean()),
                                 'per_moment_fire_rate': fires[:, sel].mean(axis=(0, 1)).tolist()}
        out['scenarios'][name] = {
            'bias_b': bias,
            'FPR_mean': float(fpr),
            'FPR_mc_se': float(nfp.std(ddof=1) / len(ids) / math.sqrt(R_FPR)),
            'FP_count_of_112': {'mean': float(nfp.mean()), 'sd': float(nfp.std(ddof=1)),
                                'p2.5': float(np.percentile(nfp, 2.5)),
                                'median': float(np.median(nfp)),
                                'p97.5': float(np.percentile(nfp, 97.5))},
            'P(FP_count <= 1)': float((nfp <= 1).mean()),
            'P(FP_count >= 1)': float((nfp >= 1).mean()),
            'per_moment_fire_rate_[m0,m1,m2]': fires.mean(axis=(0, 1)).tolist(),
            'hankel_fire_rate': float(hk.mean()),
            'per_U': per_U,
        }
    out['vectorised_rule_vs_verbatim_classifier'] = f'{verb_same}/{verb_checked} identical'
    return out


# ---------------------------------------------------------------------------------------------
# (c) exact rebuild of the sealed truncations + shared-state rescoring
# ---------------------------------------------------------------------------------------------
def sealed_era_ground_states(pr):
    """psi0 per U with ARPACK's DEFAULT random start, called for U in prereg U_grid order.
    MUST be the first v0-less eigsh calls of the process (ARPACK seeds its start generator once
    per process); main() calls this before anything else."""
    import blind_generate as bg
    out = {}
    for U in pr['system']['U_grid']:
        Hs, Js, Hc, psi_fixed = bg.build_U(pr['system']['L'], U)   # build_U itself uses a fixed v0
        e, v = eigsh(Hs, k=1, which='SA')                          # sealed-era call: no v0
        out[U] = dict(Hs=Hs, Js=Js, Hc=Hc, psi0=v[:, 0], psi_fixed=psi_fixed, E0=float(e[0]),
                      nsup=int((np.abs(v[:, 0]) ** 2 > 1e-14).sum()))
    return out


def trunc_instance(C, d, psi0):
    """the generator's truncation step (blind_generate.py:106-114), for a given psi0 tie order."""
    import spectral_lanczos as sl
    import blind_generate as bg
    Hs, Js, Hc = C['Hs'], C['Js'], C['Hc']
    prob = np.abs(psi0) ** 2
    idx = np.sort(np.argsort(prob)[::-1][:d])
    Hsub = Hs[idx][:, idx]
    Jsub = Js[idx][:, idx]
    e0, g = sl._sub_gs(Hsub)
    Hcs = (Hsub - e0 * sp.identity(d)).tocsr()
    Jg = Jsub @ g
    m_bar = np.array([np.vdot(Jg, Jg).real, np.vdot(Jg, Hcs @ Jg).real,
                      np.vdot(Jg, Hcs @ (Hcs @ Jg)).real])
    gfull = np.zeros_like(psi0)
    gfull[idx] = g
    with warnings.catch_warnings():            # local_var casts O_k psi (real-valued, complex
        warnings.simplefilter('ignore')        # dtype: J = i*A) to real; the imaginary part is 0
        var = bg.local_var(gfull, Js, Hc)
    Hc_sh = (Hs - e0 * sp.identity(Hs.shape[0])).tocsr()
    m_sh = bg.moments_of_state(gfull, Js, Hc_sh)       # full operators on the SAME state g
    return dict(idx=idx, e0=float(e0), m_bar=m_bar, var=var, m_shared=m_sh)


def replay_rng(pr, i, nsup):
    """replay the generator's RNG stream for instance i up to (and returning) the noise draw
    function; returns (U, mode, d, rng) with rng positioned just before rng.normal(0, se)."""
    cat = pr['error_catalogue']
    modes = list(cat)
    probs = np.array([cat[m]['prob'] for m in modes])
    probs /= probs.sum()
    rng = np.random.default_rng(pr['master_seed'] + 1000 + i)
    U = float(rng.choice(pr['system']['U_grid']))
    mode = str(rng.choice(modes, p=probs))
    lo, hi = cat['trunc']['d_frac_range']
    d = int(np.clip(rng.uniform(lo, hi) * nsup, 4, nsup))
    return U, mode, d, rng


def shells_of(prob, rel=1e-9):
    """group the supported configurations into shells of equal |psi0|^2."""
    order = np.argsort(prob)[::-1]
    shells, cur = [], [order[0]]
    for a in order[1:]:
        if prob[a] < 1e-14:
            break
        if abs(prob[a] - prob[cur[-1]]) <= rel * prob[cur[-1]]:
            cur.append(a)
        else:
            shells.append(cur)
            cur = [a]
    shells.append(cur)
    return shells


def part_c(pr, sha, P, V, Lb, GS):
    b = pr['battery']
    Ns = b['Ns']
    ids = sorted(i for i, l in Lb.items() if l['mode'] == 'trunc')
    rows, insts_sh, insts_naive = [], [], []
    n_mbar = n_var = n_hat = n_naive_diff = 0
    worst = {'m_bar': 0.0, 'var_loc': 0.0, 'm_hat': 0.0}
    tie_order = {}
    for U, C in GS.items():
        prob = np.abs(C['psi0']) ** 2
        tie_order[f'U={U:g}'] = [int(a) for a in np.argsort(prob)[::-1][:C['nsup']]]
        C['shells'] = shells_of(prob)
        C['cum'] = np.cumsum(np.sort(prob)[::-1]) / prob.sum()
    for i in ids:
        lab = Lb[i]
        U, d = lab['params']['U'], lab['params']['d']
        C = GS[U]
        U2, mode2, d2, rng = replay_rng(pr, i, C['nsup'])
        assert (U2, mode2, d2) == (U, 'trunc', d), f'RNG replay disagrees for id {i}'
        T = trunc_instance(C, d, C['psi0'])
        m_ex = np.array(lab['m_exact'])
        noise = rng.normal(0.0, np.sqrt(T['var'] / Ns))              # the generator's draw
        m_hat_replayed = m_ex + noise
        sb, sv, sh = (np.array(P[i][k]) for k in ('m_bar', 'var_loc', 'm_hat'))
        e_mbar = float(np.max(np.abs(T['m_bar'] - sb) / np.maximum(np.abs(sb), 1e-12)))
        e_var = float(np.max(np.abs(T['var'] - sv) / np.maximum(np.abs(sv), 1e-12)))
        e_hat = float(np.max(np.abs(m_hat_replayed - sh) / np.abs(sh)))
        worst = {k: max(worst[k], v) for k, v in zip(worst, (e_mbar, e_var, e_hat))}
        ok_m, ok_v, ok_h = e_mbar < 1e-8, e_var < 1e-8, e_hat < 1e-12
        n_mbar += ok_m
        n_var += ok_v
        n_hat += ok_h
        # naive regeneration (today's generator, fixed-v0 psi0): which instances differ?
        Tn = trunc_instance(C, d, C['psi_fixed'])
        naive_same = bool(np.allclose(Tn['m_bar'], sb, rtol=1e-7, atol=1e-9))
        n_naive_diff += not naive_same
        # shared-state rescoring on the REBUILT sealed instance, with the SEALED noise draw
        sealed_noise = sh - m_ex
        m_hat_sh = T['m_shared'] + sealed_noise
        insts_sh.append({'id': i, 'm_hat': m_hat_sh.tolist(), 'm_bar': sb.tolist(),
                         'var_loc': sv.tolist()})
        # same rescoring on the naive regeneration (cross-check of the 2026-09-26 scratch 56/56)
        _, _, _, rng2 = replay_rng(pr, i, C['nsup'])
        noise_n = rng2.normal(0.0, np.sqrt(Tn['var'] / Ns))
        insts_naive.append({'id': i, 'm_hat': (Tn['m_shared'] + noise_n).tolist(),
                            'm_bar': Tn['m_bar'].tolist(), 'var_loc': Tn['var'].tolist()})
        # boundary shell of the cut
        pos = 0
        for sh_ in C['shells']:
            if pos + len(sh_) >= d:
                cut = {'shell_size': len(sh_), 'kept_from_shell': d - pos,
                       'ambiguous': 0 < d - pos < len(sh_)}
                break
            pos += len(sh_)
        rows.append({'id': i, 'U': U, 'd': d, 'd_over_nsup': d / C['nsup'],
                     'retained_gs_weight': float(C['cum'][d - 1]), 'boundary_shell': cut,
                     'rebuild_rel_err': {'m_bar': e_mbar, 'var_loc': e_var, 'm_hat': e_hat},
                     'naive_regeneration_matches_sealed': naive_same,
                     'e0_subspace': T['e0'], 'm_exact': m_ex.tolist(), 'm_bar': sb.tolist(),
                     'm_shared_mean': T['m_shared'].tolist(), 'm_hat_shared': m_hat_sh.tolist()})
    Vsh = frozen_classify(insts_sh, sha)
    Vnv = frozen_classify(insts_naive, sha)
    # fresh-noise redraws of the shared-state estimator (vectorised rule, validated in part b)
    det_prob = []
    for r_, row in zip(insts_sh, rows):
        i = row['id']
        sv = np.array(P[i]['var_loc'])
        rng = np.random.default_rng([SEED, 3, i])
        mh = np.array(row['m_shared_mean'])[None] + np.sqrt(sv / Ns)[None] * rng.standard_normal((K_SHARED, 3))
        rej, _, _, fires = rule_vectorised(mh, np.array(P[i]['m_bar'])[None], sv[None], b)
        det_prob.append(float(rej.mean()))
        row['shared_fresh_noise_detection_prob'] = float(rej.mean())
        v = Vsh[i]
        row['shared_sealed_noise_verdict'] = {'reject': v['reject'],
                                              'r_k': [g / dl for g, dl in zip(v['g'], v['delta'])]}
    fires_sh = np.array([[g > dl for g, dl in zip(Vsh[i]['g'], Vsh[i]['delta'])] for i in ids])
    Vsd = V
    ambiguous = sum(r['boundary_shell']['ambiguous'] for r in rows)
    return {
        'rebuild': {
            'method': "psi0 from scipy eigsh(k=1, which='SA') with ARPACK's default start, first "
                      "v0-less calls of a fresh process, U in prereg U_grid order; truncation = "
                      "top-d of argsort(|psi0|^2) (numpy default sort), exactly blind_generate.py",
            'instances': len(ids),
            'cut_inside_a_degenerate_shell': ambiguous,
            'match_sealed_m_bar_(rel<1e-8)': f'{n_mbar}/{len(ids)}',
            'match_sealed_var_loc_(rel<1e-8)': f'{n_var}/{len(ids)}',
            'match_sealed_m_hat_by_RNG_replay_(rel<1e-12)': f'{n_hat}/{len(ids)}',
            'worst_rel_err': worst,
            'exact_rebuild_possible': bool(n_mbar == n_var == n_hat == len(ids)),
            'naive_regeneration_(fixed-v0 psi0)_differs_from_sealed_on': f'{n_naive_diff}/{len(ids)}',
            'fragility_note': 'the sealed tie order is set by ARPACK start-vector noise; calling '
                              'the four eigsh in another order reproduces only a minority of the '
                              'sealed instances (checked 2026-09-27: 19/56 with order 6,5,4,3). The '
                              'recovered order is stored below so later rescoring need not rely '
                              'on ARPACK state.',
            'recovered_tie_order_top_nsup_basis_indices': tie_order,
            'basis': 'index = iu*Dd + id over sorted up/down bit-strings (spectral_lanczos.strings, L=6, n=2)',
        },
        'shared_state_rescoring': {
            'estimator': 'm_hat_k = <J g|(H - e0)^k|J g> (full operators, g = truncated-subspace ground '
                         'state, e0 = <g|H|g>) + the sealed shot-noise draw (m_hat_sealed - m_exact); '
                         'var_loc as sealed (already the local variance on g); frozen rule, verbatim',
            'rejected_sealed_noise': f"{sum(Vsh[i]['reject'] for i in ids)}/{len(ids)}",
            'rejected_idealized_sealed_estimator': '56/56 (sealed record)',
            'fires_per_moment_[m0,m1,m2]': fires_sh.sum(axis=0).tolist(),
            'idealized_sealed_fires_per_moment_[m0,m1,m2]': np.array(
                [[P_g > P_d for P_g, P_d in zip(Vsd[i]['g'], Vsd[i]['delta'])] for i in ids]).sum(axis=0).tolist(),
            'rejected_by_m0': int(fires_sh[:, 0].sum()),
            'min_r0_over_56': float(min(Vsh[i]['g'][0] / Vsh[i]['delta'][0] for i in ids)),
            'idealized_sealed_min_r0_over_56': float(min(V[i]['g'][0] / V[i]['delta'][0] for i in ids)),
            'reading': 'with the shared state, m_hat_0 - m_bar_0 = <g|J Q_S J|g>, the weight J g leaks '
                       'out of the kept subspace S, so every truncation fires on m0 (r0 >= min_r0_over_56)',
            'rejected_only_off_diagonally': int((~fires_sh[:, 0] & fires_sh[:, 1:].any(axis=1)).sum()),
            'fresh_noise': {'redraws_per_instance': K_SHARED, 'rng': f'default_rng([{SEED}, 3, id])',
                            'min_detection_prob': float(min(det_prob)),
                            'mean_detection_prob': float(np.mean(det_prob)),
                            'expected_caught_of_56': float(np.sum(det_prob))},
            'cross_check_naive_regeneration_shared_state':
                f"{sum(Vnv[i]['reject'] for i in ids)}/{len(ids)} rejected (not the sealed instances)",
        },
        'per_instance': rows,
    }


# ---------------------------------------------------------------------------------------------
# (d) Krylov
# ---------------------------------------------------------------------------------------------
def part_d(pr, P, V, Lb):
    b = pr['battery']
    out = {}
    for nl in sorted({l['params']['nl'] for l in Lb.values() if l['mode'] == 'krylov'}):
        ids = sorted(i for i, l in Lb.items() if l['mode'] == 'krylov' and l['params']['nl'] == nl)
        rel = np.array([(np.array(P[i]['m_bar']) - Lb[i]['m_exact']) / np.array(Lb[i]['m_exact'])
                        for i in ids])
        r = np.array([np.array(V[i]['g']) / np.array(V[i]['delta']) for i in ids])
        shot = np.array([np.array(V[i]['g']) / (b['z'] * np.sqrt(np.array(P[i]['var_loc']) / b['Ns']))
                         for i in ids])
        out[f'nl={nl}'] = {
            'n': len(ids), 'rejected': int(sum(V[i]['reject'] for i in ids)),
            'rel_err_m_bar_vs_m_exact_min_[m0,m1,m2]': rel.min(axis=0).tolist(),
            'rel_err_m_bar_vs_m_exact_max_[m0,m1,m2]': rel.max(axis=0).tolist(),
            'max_residual_over_frozen_threshold_[r0,r1,r2]': r.max(axis=0).tolist(),
            'max_residual_over_shot_only_threshold_[m0,m1,m2]': shot.max(axis=0).tolist(),
            'would_fire_with_shot_only_threshold_(no_2pct_bias_budget)': int((shot > 1).any(axis=1).sum()),
        }
    return {'per_nl': out,
            'note': 'an nl-node Lanczos representation is exact through m_{2nl-1}: nl=1 matches m0,m1 '
                    'and misses m2; nl>=2 matches m0..m2 to rounding. Shot-only threshold = '
                    'z sqrt(var_loc/Ns) (the frozen threshold without its 2% bias term).'}


# ---------------------------------------------------------------------------------------------
# (e) calibrated scope
# ---------------------------------------------------------------------------------------------
def part_e(pr, Lb, GS, rows_c):
    cat = pr['error_catalogue']
    tr = rows_c
    ac = [l for l in Lb.values() if l['mode'] == 'ac']
    wfrac = np.array([l['params']['w_s'] / l['m_exact'][0] for l in ac])
    oms = np.array([l['params']['om_s'] for l in ac])
    band = {f'U={U:g}': None for U in GS}
    for U, C in GS.items():
        Ev = np.linalg.eigvalsh(C['Hs'].toarray())
        om = Ev - C['E0']
        band[f'U={U:g}'] = {'om_lo': float(om[om > 1e-6].min()), 'om_hi': float(om.max())}
    Ucount = defaultdict(lambda: defaultdict(int))
    for l in Lb.values():
        Ucount[l['mode']][f"U={l['params']['U']:g}"] += 1
    nl_count = defaultdict(int)
    for l in Lb.values():
        if l['mode'] == 'krylov':
            nl_count[f"nl={l['params']['nl']}"] += 1
    n_sector = GS[pr['system']['U_grid'][0]]['Hs'].shape[0]
    return {
        'system': {'L': pr['system']['L'], 'boundary': 'periodic ring (spectral_lanczos.hop, PBC)',
                   'filling': pr['system']['filling'] + ' -> N_up = N_dn = 2 (4 electrons on 6 sites)',
                   'sector_dim': int(n_sector),
                   'gs_support_nsup': int(GS[pr['system']['U_grid'][0]]['nsup']),
                   'probe': pr['system']['probe'],
                   'U_over_t_values': pr['system']['U_grid'],
                   'instances_per_mode_and_U': {m: dict(v) for m, v in Ucount.items()},
                   'battery': {k: pr['battery'][k] for k in ('Ns', 'z', 'assumed_bias_frac', 'hankel_tol')}},
        'truncations': {
            'prereg_d_frac_range': cat['trunc']['d_frac_range'],
            'realized_d_range': [min(r['d'] for r in tr), max(r['d'] for r in tr)],
            'realized_fraction_of_gs_support_kept': [min(r['d_over_nsup'] for r in tr),
                                                     max(r['d_over_nsup'] for r in tr)],
            'realized_fraction_of_gs_weight_kept': [min(r['retained_gs_weight'] for r in tr),
                                                    max(r['retained_gs_weight'] for r in tr)],
            'note': 'the kept |psi0|^2 weight depends only on d (tied shells carry equal weight)'},
        'spurious_features': {
            'prereg_spurious_weight_range_(fraction_of_m0)': cat['ac']['spurious_weight_range'],
            'realized_w_s_over_m0': [float(wfrac.min()), float(wfrac.max())],
            'realized_om_s_range': [float(oms.min()), float(oms.max())],
            'om_s_drawn_uniform_in_[om_lo, om_hi]_of_each_U': band,
            'realized_atom_share_of_reconstructed_m0': {
                'weight_changing_(w_s/(m0+w_s))': [float(x) for x in (lambda a: [a.min(), a.max()])(
                    np.array([l['params']['w_s'] / (l['m_exact'][0] + l['params']['w_s'])
                              for l in ac if not l['params'].get('preserve_m0')]))],
                'weight_preserving_(w_s/(m0+w_s), after rescaling to m0)': [float(x) for x in (lambda a: [a.min(), a.max()])(
                    np.array([l['params']['w_s'] / (l['m_exact'][0] + l['params']['w_s'])
                              for l in ac if l['params'].get('preserve_m0')]))]},
            'n_preserve_m0': int(sum(1 for l in ac if l['params'].get('preserve_m0'))),
            'n_weight_changing': int(sum(1 for l in ac if not l['params'].get('preserve_m0'))),
            'note': 'for the m0-preserving subclass the whole measure is rescaled by '
                    'm0/(m0 + w_s) after adding the atom, so the atom\'s share of m0 is '
                    'w_s/(m0 + w_s)'},
        'krylov': {'nl_counts': dict(nl_count), 'prereg_nl_choices': cat['krylov']['nl_choices']},
    }


# ---------------------------------------------------------------------------------------------
def main():
    t0 = time.time()
    sha_before = {f: sha256(os.path.join(DATA, f)) for f in SEALED_INPUTS}
    pr, sha, P, V, Lb, score = load_sealed()
    GS = sealed_era_ground_states(pr)            # FIRST v0-less eigsh calls of the process
    res = {}
    res['0_frozen_rule_check'] = part0_rule_check(pr, sha, P, V)
    res['a_primary_endpoint'] = part_a(pr, V, Lb, score)
    res['b_fpr_sensitivity'] = part_b(pr, sha, P, V, Lb)
    c = part_c(pr, sha, P, V, Lb, GS)
    res['c_truncations_rebuild_and_shared_state'] = c
    res['d_krylov'] = part_d(pr, P, V, Lb)
    res['e_calibrated_scope'] = part_e(pr, Lb, GS, c['per_instance'])
    sha_after = {f: sha256(os.path.join(DATA, f)) for f in SEALED_INPUTS}
    if sha_after != sha_before:
        raise SystemExit('a sealed input changed during the run')
    runtime = time.time() - t0
    res = {'provenance': {
        'script': 'src/blind_addenda.py', 'date': '2026-09-27',
        'seed': SEED, 'seeds_used': {'b': f'default_rng([{SEED}, 2])', 'c': f'default_rng([{SEED}, 3, id])',
                                     'sealed_generator_replay': 'default_rng(master_seed + 1000 + id)'},
        'runtime_s': round(runtime, 2),
        'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__,
        'platform': platform.platform(),
        'prereg_sha256': sha,
        'inputs_sha256': sha_after,
        'code_sha256': {f: sha256(os.path.join(HERE, f)) for f in CODE_INPUTS},
        'sealed_files_modified': False,   # checked: input hashes identical before and after the run
        'sim_only': True}, **res}
    json.dump(res, open(OUT_JSON, 'w'), indent=1)

    a, bb, cc, dd = (res['a_primary_endpoint'], res['b_fpr_sensitivity'],
                     res['c_truncations_rebuild_and_shared_state'], res['d_krylov'])
    print('frozen rule on sealed instances :', res['0_frozen_rule_check']['verbatim_classifier_reproduces_sealed_verdicts'])
    t, f = a['primary_endpoint']['TPR'], a['primary_endpoint']['FPR']
    print(f"TPR {t['k']}/{t['n']} = {t['rate']:.4f} {np.round(t['wilson95'], 4)}   "
          f"FPR {f['k']}/{f['n']} = {f['rate']:.4f} {np.round(f['wilson95'], 4)}")
    for s in a['sealed_secondary_predictions_that_failed']:
        o = s['outcome']
        print(f"  failed sealed prediction: {s['prediction']}: {o['k']}/{o['n']} {np.round(o['wilson95'], 3)}")
    for k, v in bb['scenarios'].items():
        print(f"FPR {k:45s}: {v['FPR_mean']:.4f} +- {v['FPR_mc_se']:.4f}")
    print('  vectorised vs verbatim:', bb['vectorised_rule_vs_verbatim_classifier'])
    rb = cc['rebuild']
    print('rebuild m_bar/var/m_hat:', rb['match_sealed_m_bar_(rel<1e-8)'], rb['match_sealed_var_loc_(rel<1e-8)'],
          rb['match_sealed_m_hat_by_RNG_replay_(rel<1e-12)'], ' naive differs on', rb['naive_regeneration_(fixed-v0 psi0)_differs_from_sealed_on'])
    ss = cc['shared_state_rescoring']
    print('shared-state rescoring:', ss['rejected_sealed_noise'], ' fires/moment', ss['fires_per_moment_[m0,m1,m2]'],
          ' fresh-noise min P(detect)', round(ss['fresh_noise']['min_detection_prob'], 4))
    k1 = dd['per_nl']['nl=1']
    print('krylov nl=1 m2 rel err', np.round([k1['rel_err_m_bar_vs_m_exact_min_[m0,m1,m2]'][2],
                                              k1['rel_err_m_bar_vs_m_exact_max_[m0,m1,m2]'][2]], 4),
          ' max r2', round(k1['max_residual_over_frozen_threshold_[r0,r1,r2]'][2], 4))
    print(f'wrote {OUT_JSON}  ({runtime:.1f} s)')


if __name__ == '__main__':
    main()
