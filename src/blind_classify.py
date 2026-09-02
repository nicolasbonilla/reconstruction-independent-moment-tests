# -*- coding: utf-8 -*-
"""ACTION 3 (3/4) -- CLASSIFY blind. Reads the SEALED prereg + the PUBLIC instances (m_hat,
m_bar, var_loc). NEVER opens blind_labels_sealed.json. Applies the frozen battery verbatim."""
import os, json, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, '..', '06_results'))


def main():
    blob = open(os.path.join(OUT, 'prereg.json'), 'rb').read()
    if hashlib.sha256(blob).hexdigest() != open(os.path.join(OUT, 'prereg.sha256')).read().strip():
        raise SystemExit("SEAL BROKEN.")
    pr = json.loads(blob); b = pr['battery']
    Ns, z, bias, tol = b['Ns'], b['z'], b['assumed_bias_frac'], b['hankel_tol']

    pub = json.load(open(os.path.join(OUT, 'blind_instances_public.json')))
    if pub['_prereg_sha256'] != open(os.path.join(OUT, 'prereg.sha256')).read().strip():
        raise SystemExit("instances were generated under a different prereg.")
    verdicts = []
    for inst in pub['instances']:
        m_hat = np.array(inst['m_hat']); m_bar = np.array(inst['m_bar']); var = np.array(inst['var_loc'])
        delta = z * np.sqrt(var / Ns + (bias * np.abs(m_hat)) ** 2)
        g = np.abs(m_hat - m_bar)
        moment_fire = bool(np.any(g > delta))
        # interval-Hankel: is there a feasible (m0,m1,m2) in the boxes making H1 PSD & m1>=0?
        # conservative check on the reconstruction point m_bar (what the consumer shipped):
        H1 = np.array([[m_bar[0], m_bar[1]], [m_bar[1], m_bar[2]]])
        hankel_fire = bool(np.linalg.eigvalsh(H1)[0] < -tol or m_bar[1] < -tol)
        reject = moment_fire or hankel_fire
        verdicts.append({'id': inst['id'], 'reject': reject,
                         'moment_fire': moment_fire, 'hankel_fire': hankel_fire,
                         'g': g.tolist(), 'delta': delta.tolist()})
    json.dump({'_prereg_sha256': pub['_prereg_sha256'], 'verdicts': verdicts},
              open(os.path.join(OUT, 'blind_verdicts.json'), 'w'), indent=1)
    nrej = sum(v['reject'] for v in verdicts)
    print(f"classified {len(verdicts)} instances BLIND: {nrej} REJECT, {len(verdicts)-nrej} corroborate")
    print("  wrote blind_verdicts.json   Next: python blind_score.py")


if __name__ == '__main__':
    main()
