# -*- coding: utf-8 -*-
"""ACTION 3 (1/4) -- PRE-REGISTER the blind falsification harness. SIM-ONLY.

Freezes, BEFORE any instance exists: the battery {m0,m1,m2}+interval-Hankel, the interval
thresholds delta_k (from the ibm_fez shot budget), the error catalogue (three ALGORITHMIC
causes + clean controls), and the master seed. Writes prereg.json and its SHA-256 to
prereg.sha256; blind_generate/classify re-verify that hash and refuse to run if the seal
is broken, so the battery and thresholds cannot be tuned after seeing outcomes.

HONEST SCOPE (must travel with any use of this harness):
  * NOT a field catch: every 'truth' is classically known here (L=6); the corruptions are
    algorithmic errors we injected, not errors of genuinely-unknown-to-science provenance.
  * NOT device-validated: shot noise is modelled from the device-side local-estimator
    variance at the ibm_fez budget; there is no QPU in the loop.
  * NOT robust to unanticipated modes: by construction the battery is BLIND to any
    corruption that preserves m0,m1,m2 (moment-preserving analytic-continuation artefacts
    and >=2-node under-converged Krylov). The harness is designed to EXPOSE that blindness,
    not hide it. A high miss-rate on those modes is the honest result, not a failure.
"""
import os, json, hashlib

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, '..', 'data'))
os.makedirs(OUT, exist_ok=True)

PREREG = {
    'title': 'Blind pre-registered falsification harness for reconstruction-independent moment screen',
    'version': 'battery-v1',
    'sim_only': True,
    'system': {'L': 6, 'U_grid': [3.0, 4.0, 5.0, 6.0], 'filling': 'nup=nd=L//3 (doped)',
               'probe': 'current operator J; T=0 Lehmann measure'},
    'battery': {
        'moments': ['m0', 'm1', 'm2'],
        'rule': 'REJECT iff |m_hat_k - m_bar_k| > delta_k for any k in {0,1,2} '
                'OR interval-Hankel H1=[[m0,m1],[m1,m2]] infeasible (min eig < -tol) '
                'OR m1 interval entirely < -tol.',
        'shot_model': 'm_hat_k ~ Normal(m_exact_k, sigma_k^2); sigma_k = z * local-estimator '
                      'standard error of the DEVICE-SIDE PREPARED state at Ns shots (honesty '
                      'fix: prepared, possibly-truncated, state -- NOT the exact |0>).',
        'Ns': 50000, 'z': 1.96, 'assumed_bias_frac': 0.02, 'hankel_tol': 1e-6,
    },
    'error_catalogue': {
        'clean':  {'prob': 0.40, 'desc': 'faithful reconstruction; m_bar = m_exact (control)'},
        'trunc':  {'prob': 0.20, 'desc': 'SQD/determinant-subspace truncation; m_bar = moments '
                                         'of the truncated-subspace reconstruction; d ~ U[0.15,0.75]*nsup',
                   'd_frac_range': [0.15, 0.75]},
        'krylov': {'prob': 0.20, 'desc': 'under-converged Krylov; m_bar = nl-node Lanczos '
                                         'representation (exact through m_{2nl-1}); nl in {1,2,3,4}',
                   'nl_choices': [1, 2, 3, 4]},
        'ac':     {'prob': 0.20, 'desc': 'analytic-continuation artefact; m_bar = m_exact + a '
                                         'spurious Gaussian atom (random freq/weight), renormalized '
                                         'to preserve m0 with prob 0.5 (the moment-preserving subclass '
                                         'is designed-blind).',
                   'spurious_weight_range': [0.02, 0.25], 'preserve_m0_prob': 0.5},
    },
    'n_instances': 300,
    'master_seed': 20260824,
    'decision_prereg': 'Primary endpoint: TPR (reject | corrupted) and FPR (reject | clean) at '
                       'the frozen thresholds. Secondary: per-mode catch rate, expected to be '
                       'high for trunc / nl=1 krylov / non-preserving ac and LOW (by design) for '
                       'nl>=2 krylov and m0-preserving ac. FPR target ~ nominal joint 1-(1-0.05)^3.',
}


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(',', ':')).encode('utf-8')


def main():
    for _p in ('prereg.json', 'prereg.sha256'):   # the committed seal is a read-only record
        if os.path.exists(os.path.join(OUT, _p)):
            raise SystemExit(f"REFUSING to overwrite the sealed record {_p}; "
                             "run in a scratch copy of the repository to re-seal.")
    blob = canonical(PREREG)
    h = hashlib.sha256(blob).hexdigest()
    with open(os.path.join(OUT, 'prereg.json'), 'wb') as f:
        f.write(blob)
    with open(os.path.join(OUT, 'prereg.sha256'), 'w') as f:
        f.write(h + '\n')
    print("PRE-REGISTRATION SEALED")
    print(f"  n_instances = {PREREG['n_instances']}   master_seed = {PREREG['master_seed']}")
    print(f"  battery = {PREREG['battery']['moments']} + interval-Hankel")
    print(f"  modes   = {list(PREREG['error_catalogue'])}")
    print(f"  SHA-256 = {h}")
    print(f"  wrote {os.path.join(OUT,'prereg.json')}  (+ .sha256)")
    print("\nNext: python blind_generate.py   (uses ONLY the sealed prereg)")


if __name__ == '__main__':
    main()
