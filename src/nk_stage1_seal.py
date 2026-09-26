# -*- coding: utf-8 -*-
"""STAGE 1 — seal manifest_nk_device_v1.json for the on-device n_k discriminating-falsifier job.

Assembles every frozen quantity mandated by the audit/adjudication chain (wf_348330ab-00a R7-R11,
wf_90a0b575-9ab, wf_e79137d4-be0) from the PASSED Stage-0 record (v3 + v3-R2 QUOTABLE-PASS),
computes the manifest SHA-256, and writes both. The device job may only run against this seal.
"""
import os, sys, json, time, hashlib, io
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from nk_stage0_gate import (L, U, FN, THETA, PHI, FILL, SHOTS, ed_references, ft_matrix,
                            build_circuits)
from nk_stage0_gate_v3 import build_isa, T1_FLOOR, T2_FLOOR
from qiskit import qpy

RES = os.path.normpath(os.path.join(HERE, '..', 'data'))

def sha(b): return hashlib.sha256(b).hexdigest()
def fsha(p): return sha(open(p, 'rb').read())

def main():
    if os.path.exists(os.path.join(RES, 'manifest_nk_device_v1.json')):   # sealed, read-only
        raise SystemExit('REFUSING to overwrite the sealed manifest_nk_device_v1.json; '
                         'run in a scratch copy of the repository to re-seal.')
    refs, _ = ed_references()
    W = ft_matrix().conj().T
    circs = build_circuits(W)
    val_map = {'rung0': 'rung0_nk', 'mirrorFT': 'rung0_nk', 'rungB': 'rungB_nk',
               'calibB2': 'calibB2_nk', 'calibFT': 'rung0_nk'}
    isa, czc, meta, dt = build_isa(circs, refs, val_map)
    qpy_hashes = {}
    for nm, qc in isa.items():
        buf = io.BytesIO(); qpy.dump(qc, buf)
        qpy_hashes[nm] = sha(buf.getvalue())

    v3  = json.load(open(os.path.join(RES, '2026-09-02_nk_stage0_gate_v3.json')))
    r2  = json.load(open(os.path.join(RES, '2026-09-02_nk_stage0_v3R2.json')))
    assert r2['verdict'][0] == 'QUOTABLE-PASS', 'manifest may only be sealed on a passed gate'

    maxdev_B = float(np.max(np.abs(refs['rungB_nk'] - FILL)))       # 0.545717...
    LAM = 0.50
    dn = LAM*maxdev_B                                                # Delta_n(0.50)
    Ttilde = refs['rungB_Ttilde']

    manifest = {
      'name': 'manifest_nk_device_v1', 'version': 1, 'sealed_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
      'sealed_by': ['wf_348330ab-00a (gate audit)', 'wf_90a0b575-9ab (v2 adjudication)',
                    'wf_e79137d4-be0 (threshold adjudication + v3-R2 protocol)',
                    'v3 10/10 PASS via v3-R2 QUOTABLE-PASS (mean 0.023850, UCB 0.026018 < 0.0266667)'],
      'experiment': {
        'system': f'L={L} PBC Hubbard, doped 2up/2dn (filling 1/3), U/t={U:g}, 12 qubits',
        'backend': 'ibm_fez (Heron r2)', 'claim_level_lambda': LAM,
        'claim': 'pre-registered on-device discriminating falsifier: the mitigated device momentum '
                 'distribution accepts the honest ED reference and rejects dispersion-flattening '
                 'corruptions of strength lambda >= 0.50, at calibration-staging size (all references '
                 'ED-exact by design; this is the diagonal channel of the paper, not the off-diagonal '
                 'moment screen, which stays deferred).',
      },
      'state_prep': {'theta': THETA, 'phi': PHI, 'note': 'phi layer is free momentum phases (R2 theorem)',
                     'rungB_overlap2_gs': refs['rungB_overlap2_gs']},
      'frozen_refs_12dig': {k: ([round(float(x), 12) for x in v] if isinstance(v, np.ndarray) else round(float(v), 12))
                            for k, v in refs.items()},
      'kinetic_anchor': {'Ttilde_ED_rungB': round(float(Ttilde), 12),
                         'identity': 'm_{-1} = (1/2)<-T> - D; D (Kohn stiffness) frozen from the '
                                     'committed run_inverse_moment_falsifier.py at this size; the device '
                                     'delivers <T> only ("kinetic f-sum anchor", never "measured m_-1"); '
                                     'derived scalar, NOT an independent falsifier (linear contraction '
                                     'of the same counts)',
                         'tau_T': round(0.5*abs(float(Ttilde))*LAM, 6)},
      'isa': {'recipe': 'FakeFez heavy-hex, sabre, optimization_level=3, seed-scan 0-4 (min CZ), '
                        'ASAP scheduling, compaction to active qubits; measurements added pre-transpile',
              'qpy_sha256': qpy_hashes, 'cz_counts': czc,
              'meta': {k: {kk: vv for kk, vv in v.items()} for k, v in meta.items()},
              'note': 'day-of device transpile re-runs the same recipe onto the accepted physical chain; '
                      'the actual job QPY hashes are recorded at submission time alongside this seal'},
      'noise_floors': {'T1_min_us': 150, 'T2_min_us': 100,
                       'day_of_chain_rule': 'accept the chain iff min over chain qubits of day-of '
                                            'calibrated T1_q >= 150 us AND T2_q >= 100 us AND no chain '
                                            'CZ error > 2x device median; otherwise re-select the chain '
                                            'or forfeit the window'},
      'job_spec': {'circuits': ['rung0', 'mirrorFT', 'calibFT', 'rungB', 'calibB2'],
                   'shots_per_circuit': SHOTS,
                   'twirling': 'GATE-LEVEL Pauli twirling, num_randomizations=32, independently '
                               'compiled (SamplerV2 twirling enable_gates=True) — a hard requirement '
                               '(untwirled coherent ZZ breaks the calibration model: p 0.143 vs 0.060); '
                               'TREX/measure twirling additionally enabled; DD XpXm',
                   'batch': 'one Batch; shot-block interleaving of all circuits in randomized order '
                            '(twin is blind to inter-run drift by construction)',
                   'raw_counts': 'retained and committed', 'qpu_minutes_budget': '3-4 of one 10-min window',
                   'repeat_policy': 'one repeat in the NEXT window only on pre-registered inconsistency'},
      'estimator_freeze': {
        'primary_iii': 'calibB2 twin least-squares flattening fit (sole verdict estimator; clip [0,1])',
        'iv_D': 'p measured on calibFT (pulse-identical acquisition of mirrorFT ISA), applied to '
                'mirrorFT, post-mitigation residual vs ED; threshold Delta_n(0.20)/5 = 0.0266667 '
                'WITH the adjudicated band rule: |B - thr| <= 1.96*sigma_B_hat is INDETERMINATE '
                '(reported, does not void the headline); only above-band = model-chain failure',
        'diagnostics_nonverdict': 'CZ-ratio transfer, mirror survival, halfB scaling — recorded whatever they say',
        'no_post_hoc_swap': 'if (iii) fails on device the pre-committed reading is model-chain failure, '
                            'reported as such — never a fallback estimator',
        'estimator_history_disclosure': 'mirror-survival and calibB (Slater) were tried and rejected '
                                        'pre-seal for stated physical defects; both fail the binding '
                                        'sim config; Stage-0 is a design result — the pre-registered '
                                        'test is this device run alone',
        'twin_disclosure': 'calibB2 is NOT independent: a deliberately maximally-matched noise replica, '
                           'pulse-identical to rungB up to virtual-Z frames; the falsifier is by design '
                           'insensitive to corruption components proportional to (n_ref - fill); '
                           'correctness of the mitigated device leg is certified solely by the residual '
                           'check vs ED (criterion i), and (iii) functions as a residual-coherence detector',
      },
      'criteria': {
        'i':  f'max_k|ntilde_k - n_k^ED(rungB)| < Delta_n(0.50)/5 = {dn/5:.6f}',
        'ii': f'Delta_n(0.50) - B > 1.96*sigma_eff (sigma_shot from actual kept counts; sigma_p=0.02)',
        'iii': '|p_hat(calibB2) - p_fit(rungB vs ED)| < 0.02',
        'iv': 'see estimator_freeze.iv_D (band rule pre-registered)',
        'falsifier_decision': {
          'tau': round(dn/2.0, 6),
          'rule': f'REJECT a reconstruction R iff max_k|ntilde_k - nbar_k(R)| > tau = Delta_n(0.50)/2 = {dn/2:.6f}; '
                  'the honest ED reference must be ACCEPTED (statistic < tau) and every C1 corruption at '
                  'lambda >= 0.50 must be REJECTED; C1: nbar_k(lambda) = (1-lambda)*n_k^ED + lambda*fill',
          'gate_forecast': 'sim: honest statistic ~0.02 << tau ~0.136 << corrupted ~0.25 (margins ~6x both sides)'},
        'kill_rules': ['keep_frac(rungB) < 0.25 -> abort claim, report bias floor',
                       'day-of chain rule fails -> forfeit window',
                       'canary above-band -> model-chain failure branch',
                       'no post-hoc estimator or threshold changes of any kind'],
        'never_quote': 'lambda=0.20 margins are inside shot noise and are never quotable (audit 1)',
      },
      'outcome_paragraphs': {
        'success': 'On ibm_fez we measured the doped-sector momentum distribution n_k and its kinetic '
                   'f-sum anchor of a 12-qubit correlated preparation (overlap^2 = 0.92 with the U/t=8 '
                   'doped ground state) through the number-conserving fermionic Fourier transform, and '
                   'the pre-registered falsifier REJECTED dispersion-flattening corruptions at '
                   'lambda >= 0.50 while ACCEPTING the honest reference — mitigation bias measured, not '
                   'assumed; every reference ED-exact by design (calibration size, stated throughout); '
                   'the preparation one interaction layer beyond matchgates, classically simulable. This '
                   'substantially addresses the dominant objection for the diagonal channel at '
                   'calibration size; the off-diagonal moment screen stays deferred.',
        'negative': 'The pre-registered on-device test did not discriminate at the sealed level: we '
                    'report the measured post-mitigation bias floor of the diagonal channel versus '
                    'corruption strength, the device lambda*_min, and the gate-vs-device comparison, '
                    'as sealed before the run. This is the pre-registered negative branch, published '
                    'as such.'},
      'mandatory_report_line': 'At manifest floors the iv-D canary runs at 74-100% of its threshold '
                               '(headroom <= 26%, within ~1 sigma); on hardware it may trip '
                               'stochastically; sub-shot-noise trips are indeterminate by the sealed '
                               'band rule, and a v4-style claim-level alignment is the pre-registered '
                               'path for any future deployment.',
      'analysis_code_sha256': {os.path.basename(p): fsha(os.path.join(HERE, p)) for p in
                               ('nk_stage0_gate.py', 'nk_stage0_gate_v2.py', 'nk_stage0_gate_v3.py',
                                'nk_stage0_v3R2.py', 'nk_falsifier.py', 'spectral_lanczos.py')},
      'gate_record_sha256': {'v3': fsha(os.path.join(RES, '2026-09-02_nk_stage0_gate_v3.json')),
                             'v3R2': fsha(os.path.join(RES, '2026-09-02_nk_stage0_v3R2.json'))},
      'credentials_note': 'the device run uses the author\'s OWN freshly-created IBM credentials; the '
                          'previously exposed token must be revoked before the run; no token is ever committed',
    }
    mpath = os.path.join(RES, 'manifest_nk_device_v1.json')
    blob = json.dumps(manifest, indent=1, sort_keys=True).encode()
    open(mpath, 'wb').write(blob)
    h = sha(blob)
    open(os.path.join(RES, 'manifest_nk_device_v1.sha256'), 'w').write(h + '\n')
    print('sealed', mpath)
    print('SHA-256:', h)

if __name__ == '__main__':
    main()
