# -*- coding: utf-8 -*-
"""ACTION 3 (2/4) -- GENERATE blind instances from the SEALED prereg. SIM-ONLY.

Real physics per mode: 'trunc' builds an actual truncated determinant subspace and takes its
reconstruction moments; 'krylov' uses an nl-node Lanczos (Haydock) representation (moments
exact through 2nl-1); 'ac' adds a spurious atom to the exact measure. The independent estimate
m_hat_k ~ N(m_exact_k, se_k^2) has its true mean at the EXACT full-operator moment (the
reconstruction-independent quantity); the corruption lives entirely in the reconstruction
back-moment m_bar_k. Shot standard error se_k and the local-estimator variance are taken from
the DEVICE-SIDE PREPARED state (honesty fix).

Writes instances_public.json (m_hat, m_bar, var_loc -- NO labels) and labels_sealed.json
(mode + params + is_corrupted). The classifier reads only the public file.
"""
import os, sys, json, hashlib
import numpy as np
import scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import spectral_lanczos as sl
OUT = os.path.normpath(os.path.join(HERE, '..', 'data'))


def load_sealed_prereg():
    blob = open(os.path.join(OUT, 'prereg.json'), 'rb').read()
    h = hashlib.sha256(blob).hexdigest()
    want = open(os.path.join(OUT, 'prereg.sha256')).read().strip()
    if h != want:
        raise SystemExit("SEAL BROKEN: prereg.json hash != prereg.sha256. Re-run blind_preregister.py.")
    return json.loads(blob)


def build_U(L, U):
    nup = nd = L // 3
    Tu, Su, iu = sl.hop(L, nup); Td, Sd, idd = sl.hop(L, nd)
    Du, Dd = len(Su), len(Sd)
    upocc = sl.occ_matrix(Su, L); dnocc = sl.occ_matrix(Sd, L)
    Hs = (sp.kron(Tu, sp.identity(Dd)) + sp.kron(sp.identity(Du), Td)
          + sp.diags(U * (upocc @ dnocc.T).ravel())).tocsr()
    Ju = sl.current_string(L, nup); Jd = sl.current_string(L, nd)
    Js = (sp.kron(Ju, sp.identity(Dd)) + sp.kron(sp.identity(Du), Jd)).tocsr()
    E0, psi0 = sl._sub_gs(Hs)
    Hc = (Hs - E0 * sp.identity(Hs.shape[0])).tocsr()
    return Hs, Js, Hc, psi0


def moments_of_state(psi, Js, Hc):
    """full-operator moments m_k=<phi|Hc^k|phi>, phi=Js psi, k=0,1,2 (psi need not be normalized)."""
    phi = Js @ psi
    Hphi = Hc @ phi; H2phi = Hc @ Hphi
    return np.array([float(np.vdot(phi, phi).real),
                     float(np.vdot(phi, Hphi).real),
                     float(np.vdot(phi, H2phi).real)])


def local_var(psi, Js, Hc):
    """device-side local-estimator variance of O_k on the PREPARED state psi (honesty fix)."""
    p = np.abs(psi) ** 2; s = p.sum()
    if s < 1e-300:
        return np.array([0.0, 0.0, 0.0])
    p = p / s; mask = np.abs(psi) > 1e-12
    phi = Js @ psi
    Oks = [Js @ phi, Js @ (Hc @ phi), Js @ (Hc @ (Hc @ phi))]   # O_k psi
    out = []
    for Ok in Oks:
        loc = np.zeros_like(psi); loc[mask] = Ok[mask] / psi[mask]
        mean = float(np.sum(p * loc)); out.append(max(float(np.sum(p * loc ** 2) - mean ** 2), 0.0))
    return np.array(out)


def krylov_moments(psi0, Js, Hc, nl):
    phi = Js @ psi0
    theta, wj = sl.haydock_poles(lambda v: Hc @ v, phi, nl)
    return np.array([float(np.sum(wj)), float(np.sum(wj * theta)), float(np.sum(wj * theta ** 2))])


def main():
    for _p in ('blind_instances_public.json', 'blind_labels_sealed.json'):   # sealed, read-only
        if os.path.exists(os.path.join(OUT, _p)):
            raise SystemExit(f"REFUSING to overwrite the sealed record {_p}; "
                             "run in a scratch copy of the repository to regenerate.")
    pr = load_sealed_prereg()
    L = pr['system']['L']; Ugrid = pr['system']['U_grid']
    G = pr['n_instances']; seed0 = pr['master_seed']; Ns = pr['battery']['Ns']
    cat = pr['error_catalogue']
    modes = list(cat); probs = np.array([cat[m]['prob'] for m in modes]); probs /= probs.sum()

    # precompute per-U truth
    cache = {}
    for U in Ugrid:
        Hs, Js, Hc, psi0 = build_U(L, U)
        m_ex = moments_of_state(psi0, Js, Hc)
        var_ex = local_var(psi0, Js, Hc)
        # exact poles (for AC spurious-atom frequency range)
        Ev, Vv = np.linalg.eigh(Hs.toarray())
        om = Ev - float(np.vdot(psi0, Hs @ psi0).real)
        cache[U] = dict(Hs=Hs, Js=Js, Hc=Hc, psi0=psi0, m_ex=m_ex, var_ex=var_ex,
                        om_lo=float(om[om > 1e-6].min()), om_hi=float(om.max()),
                        nsup=int((np.abs(psi0) ** 2 > 1e-14).sum()))
    print(f"truth cached for U in {Ugrid} (L={L})")

    public, sealed = [], []
    for i in range(G):
        rng = np.random.default_rng(seed0 + 1000 + i)
        U = float(rng.choice(Ugrid)); C = cache[U]
        Hs, Js, Hc, psi0, m_ex = C['Hs'], C['Js'], C['Hc'], C['psi0'], C['m_ex']
        mode = str(rng.choice(modes, p=probs))
        params = {'U': U}
        if mode == 'clean':
            m_bar = m_ex.copy(); var_loc = C['var_ex'].copy(); corrupted = False
        elif mode == 'trunc':
            lo, hi = cat['trunc']['d_frac_range']
            d = int(np.clip(rng.uniform(lo, hi) * C['nsup'], 4, C['nsup']))
            prob = np.abs(psi0) ** 2; idx = np.sort(np.argsort(prob)[::-1][:d])
            Hsub = Hs[idx][:, idx]; Jsub = Js[idx][:, idx]
            e0, g = sl._sub_gs(Hsub); Hcs = (Hsub - e0 * sp.identity(d)).tocsr()
            Jg = Jsub @ g
            m_bar = np.array([float(np.vdot(Jg, Jg).real), float(np.vdot(Jg, Hcs @ Jg).real),
                              float(np.vdot(Jg, Hcs @ (Hcs @ Jg)).real)])
            gfull = np.zeros_like(psi0); gfull[idx] = g          # device-side prepared (truncated) state
            var_loc = local_var(gfull, Js, Hc); corrupted = True; params['d'] = d
        elif mode == 'krylov':
            nl = int(rng.choice(cat['krylov']['nl_choices']))
            m_bar = krylov_moments(psi0, Js, Hc, nl)
            var_loc = C['var_ex'].copy(); corrupted = True; params['nl'] = nl
        else:  # ac
            lo, hi = cat['ac']['spurious_weight_range']
            w_s = float(rng.uniform(lo, hi)) * m_ex[0]
            om_s = float(rng.uniform(C['om_lo'], C['om_hi']))
            m_bar = m_ex + w_s * np.array([1.0, om_s, om_s ** 2])
            if rng.random() < cat['ac']['preserve_m0_prob']:
                m_bar = m_bar * (m_ex[0] / m_bar[0])             # preserve m0 (still shifts m1,m2)
                params['preserve_m0'] = True
            var_loc = C['var_ex'].copy(); corrupted = True
            params['w_s'] = w_s; params['om_s'] = om_s
        # independent estimate: true mean at m_exact, device-side shot noise
        se = np.sqrt(var_loc / Ns)
        m_hat = m_ex + rng.normal(0.0, se)
        public.append({'id': i, 'm_hat': m_hat.tolist(), 'm_bar': m_bar.tolist(),
                       'var_loc': var_loc.tolist()})
        sealed.append({'id': i, 'mode': mode, 'corrupted': bool(corrupted), 'params': params,
                       'm_exact': m_ex.tolist()})

    json.dump({'_prereg_sha256': open(os.path.join(OUT, 'prereg.sha256')).read().strip(),
               'instances': public}, open(os.path.join(OUT, 'blind_instances_public.json'), 'w'), indent=1)
    json.dump({'labels': sealed}, open(os.path.join(OUT, 'blind_labels_sealed.json'), 'w'), indent=1)
    ncorr = sum(s['corrupted'] for s in sealed)
    print(f"generated {G} instances: {ncorr} corrupted, {G-ncorr} clean")
    print(f"  wrote blind_instances_public.json (no labels) + blind_labels_sealed.json")
    print("Next: python blind_classify.py")


if __name__ == '__main__':
    main()
