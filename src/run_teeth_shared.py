# -*- coding: utf-8 -*-
r"""
THE TEETH EXPERIMENT, SHARED-STATE VERSION (classical post-processing of the subspace state; no shots).
=========================================================================================================
Addresses a referee-style objection to Fig.~teeth: the teeth demo evaluated the reconstruction-INDEPENDENT
first moment on the FULL EXACT ground state psi0 (of the whole sector), while only the reconstruction used
the truncated subspace. That is more independent than one prepared state allows: here BOTH branches are
computed from ONE state, the subspace ground state.

The quantity is the operator "leak" of the paper's own Eq.(app-leak):
evaluate BOTH moments on the SAME subspace ground state |0_S> (the ground state OF the truncated
subspace Hamiltonian H_S = Pi_S H Pi_S, NOT the full-sector GS), with a SHARED energy reference
e0 = <0_S|H|0_S>:

    m_hat_1  = <0_S| J (H - e0)             J |0_S>     INDEPENDENT estimator: FULL sparse J and H
                                                        act on |0_S> and reach determinants OUTSIDE S
                                                        (classical post-processing of the shared state)
    m_bar_1  = <0_S| J_S (H_S - e0) J_S |0_S>           RECONSTRUCTION moment: within-S Lehmann,
             = \int omega A(omega) domega,              J_S = Pi_S J Pi_S, H_S = Pi_S H Pi_S
    Delta_1^shared = | m_hat_1 - m_bar_1 |              exactly Eq.(app-leak), k=1, on the shared |0_S>.

This is the paper's own appendix definition (Eq. app-leak): |0> there is the subspace ground state,
e0 the subspace ground energy. The earlier demo used psi0 (full GS) for the independent branch; here
both branches use |0_S>. Everything is evaluated exactly and classically on |0_S> (full sparse J and H):
no shots are drawn, no device is involved, and m_hat_1 is not a moment measured from samples. On a device
it would be a ground-state expectation of the prepared state, which needs rotated-basis measurements and
carries device noise.
Truncation order: np.argsort(|psi0|^2); ties in the degenerate |psi0|^2 shells are broken by
floating-point noise (the L=12 teeth were left on this order on purpose; see the interval scripts for the
deterministic lexsort key).
Reproducibility (checked 2026-09-27 in a scratch copy, scipy 1.13.1 / numpy 1.26.4): a rerun reproduces the
two 5% crossings quoted in the paper (d ~ 15977 idealized and d ~ 43554 shared, relative change < 2e-6);
the sweep points agree to <= 2.4e-4 (relative, every column) for d >= 460, but the three smallest d move
(m_bar by 1.4%, 0.4% and 6.2% at d = 337, 246, 180; worst rel_shared 0.9983 -> 0.9982), because the tie order
and the subspace eigensolves are not bit-reproducible against the 2026-08-26 run. The committed data/2026-08-26_teeth_shared_state.json, paper/figs/teeth_shared.dat
and src/cache/teeth_shared_L12.npz are therefore kept as the record and were NOT regenerated; that JSON
still carries the pre-2026-09-27 wording ('device-realizable', 'closes', and a verdict about the threshold
'a device actually enforces'), which the strings written below replace.

We ALSO carry the published (idealized) curve for side-by-side plotting:
    Delta_1^ideal(rel) = | m_bar_1 - m1_op | / | m1_op |,   m1_op = <psi0| J (H - E0) J | psi0 >
with psi0 the FULL-sector exact ground state (what the old demo used).

L=12 doped ring (nup=nd=4, 2/3 filling), U/t=8 -- identical to the teeth figure. REAL numbers only.
"""
import os, sys, json, time
import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import spectral_lanczos as sl

t0 = time.time()
log = lambda *a: print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)

CACHE = os.path.join(HERE, 'cache')
RES = os.path.normpath(os.path.join(HERE, '..', 'data'))
# write the plot data into the paper's figs dir
FIGDIRS = [
    os.path.normpath(os.path.join(HERE, '..', 'paper', 'figs')),
]
DATE = '2026-08-26'


def run(L=12, U=8.0, t=1.0, npts=24, dmin=180):
    nup = nd = L // 3                                     # 2/3 filling, matches fig_teeth
    Tu, Su, iu = sl.hop(L, nup, t); Td, Sd, idd = sl.hop(L, nd, t)
    Du, Dd = len(Su), len(Sd)
    upocc = sl.occ_matrix(Su, L); dnocc = sl.occ_matrix(Sd, L)
    diagU = U * (upocc @ dnocc.T).ravel()
    Hs = (sp.kron(Tu, sp.identity(Dd), format='csr') + sp.kron(sp.identity(Du), Td, format='csr')
          + sp.diags(diagU)).tocsr()
    Ju = sl.current_string(L, nup, t); Jd = sl.current_string(L, nd, t)
    Js = (sp.kron(Ju, sp.identity(Dd), format='csr') + sp.kron(sp.identity(Du), Jd, format='csr')).tocsr()
    log(f"[shared] L={L} U={U} filling={2*nup}/{L}: sector dim={Du*Dd}, H nnz={Hs.nnz}, J nnz={Js.nnz}")

    # full-sector exact GS: only needed for the IDEALIZED (old demo) reference m1_op
    E0, psi0 = sl._sub_gs(Hs)
    Jp = Js @ psi0
    m1_op = float(np.real(np.vdot(Jp, Hs @ Jp)) - E0 * np.real(np.vdot(Jp, Jp)))   # independent on FULL psi0
    prob = np.abs(psi0) ** 2
    order = np.argsort(prob)[::-1]
    nsup = int((prob > 1e-14).sum()); order = order[:nsup]
    log(f"[shared] full-sector E0={E0:.6f}  m1_op(idealized, on psi0)={m1_op:.6f}  n_support={nsup}")

    ds = sl._logspace_d(nsup, npts, dmin)                # log grid over the FULL fig_teeth range (d>=180)
    rows = []
    for d in ds:
        d = int(d); idx = np.sort(order[:d])
        Hsub = Hs[idx][:, idx]; Jsub = Js[idx][:, idx]
        e0, g = sl._sub_gs(Hsub)                          # |0_S>, subspace ground state + energy (SHARED ref)

        # RECONSTRUCTION moment m_bar_1 (within-S Lehmann) == the old m1_trunc
        vsub = Jsub @ g
        m_bar = float(np.real(np.vdot(vsub, Hsub @ vsub)) - e0 * np.real(np.vdot(vsub, vsub)))

        # INDEPENDENT estimator m_hat_1 on the SAME |0_S>: embed g into the full sector, apply FULL J,H
        gfull = np.zeros(Du * Dd, dtype=g.dtype)
        gfull[idx] = g                                   # |0_S> as a vector of the full Hilbert space
        vful = Js @ gfull                                # J|0_S>: reaches OUTSIDE S (the genuine leak)
        m_hat = float(np.real(np.vdot(vful, Hs @ vful)) - e0 * np.real(np.vdot(vful, vful)))

        # leak diagnostic: fraction of |J|0_S>|^2 that lands outside S
        n_in = float(np.real(np.vdot(vsub, vsub)))
        n_tot = float(np.real(np.vdot(vful, vful)))
        leak_frac = max(0.0, 1.0 - n_in / n_tot) if n_tot > 0 else 0.0

        d_shared_abs = abs(m_hat - m_bar)                        # Eq.(app-leak), k=1, shared |0_S>
        rel_shared = d_shared_abs / abs(m_hat) if m_hat != 0 else 0.0
        rel_ideal = abs(m_bar - m1_op) / abs(m1_op)              # OLD published curve (m_bar vs full-GS)

        cov = d / nsup
        rows.append(dict(d=d, cov=cov, m_hat=m_hat, m_bar=m_bar, m1_op=m1_op,
                         delta_shared_abs=d_shared_abs, rel_shared=rel_shared,
                         rel_ideal=rel_ideal, leak_frac=leak_frac))
        log(f"[shared]  d={d:6d} cov={100*cov:5.1f}%  m_hat={m_hat:8.4f}  m_bar={m_bar:8.4f}  "
            f"Dshared={d_shared_abs:7.4f} ({100*rel_shared:6.2f}%)  ideal={100*rel_ideal:6.2f}%  leak={100*leak_frac:5.2f}%")
    return rows, m1_op, nsup, float(E0)


def crossing(rows, key, thr=0.05):
    """log-linear interpolate the determinant count d where the given relative curve first crosses thr."""
    for k in range(1, len(rows)):
        a, b = rows[k - 1][key], rows[k][key]
        if a < thr <= b:
            la, lb = np.log10(rows[k - 1]['d']), np.log10(rows[k]['d'])
            dcx = 10 ** (la + (lb - la) * (thr - a) / (b - a))
            return float(dcx)
    return None


def main():
    rows, m1_op, nsup, E0 = run()
    d_arr = np.array([r['d'] for r in rows])
    cov = np.array([r['cov'] for r in rows])
    rel_shared = np.array([r['rel_shared'] for r in rows])
    rel_ideal = np.array([r['rel_ideal'] for r in rows])
    mhat = np.array([r['m_hat'] for r in rows])
    mbar = np.array([r['m_bar'] for r in rows])
    leak = np.array([r['leak_frac'] for r in rows])

    cx_ideal = crossing(rows, 'rel_ideal')
    cx_shared = crossing(rows, 'rel_shared')

    os.makedirs(CACHE, exist_ok=True)
    np.savez(os.path.join(CACHE, 'teeth_shared_L12.npz'),
             L=12, U=8.0, m1_op=m1_op, nsup=nsup, E0=E0,
             d=d_arr, cov=cov, rel_shared=rel_shared, rel_ideal=rel_ideal,
             m_hat=mhat, m_bar=mbar, leak_frac=leak)

    # --- plot data file (columns match fig_teeth conventions; percentages) ---
    header = "d cov rideal rshared leak\n"
    body = "".join(f"{r['d']} {100*r['cov']:.4f} {100*r['rel_ideal']:.4f} "
                   f"{100*r['rel_shared']:.4f} {100*r['leak_frac']:.4f}\n" for r in rows)
    for fd in FIGDIRS:
        os.makedirs(fd, exist_ok=True)
        with open(os.path.join(fd, 'teeth_shared.dat'), 'w') as f:
            f.write(header); f.write(body)
        log(f"[shared] wrote {os.path.join(fd, 'teeth_shared.dat')}")

    # --- results JSON (provenance + honest verdict) ---
    out = {
        '_provenance': {
            'script': 'src/run_teeth_shared.py', 'date': DATE,
            'params': 'doped Hubbard ring L=12, U/t=8, nup=nd=4 (2/3 filling)',
            'addresses': 'referee-style objection to Fig. teeth: the independent estimator was evaluated '
                         'on the FULL-sector exact GS psi0 while only the reconstruction used the truncated '
                         'subspace',
            'scope': 'classical, exact evaluation on the subspace ground state |0_S>; no shots, no device',
            'quantity': 'Eq.(app-leak), k=1, evaluated on the SAME subspace ground state |0_S> for BOTH '
                        'branches with a shared energy reference e0=<0_S|H|0_S>.',
            'definitions': {
                'm_hat_1_shared': '<0_S| J (H - e0)   J |0_S>   (FULL sparse J,H act on |0_S>; reach outside S)',
                'm_bar_1':        '<0_S| J_S (H_S-e0) J_S|0_S> = int omega A(omega) domega   (within-S Lehmann)',
                'Delta_1_shared': '| m_hat_1_shared - m_bar_1 |   (operator leak on the shared state; classical, no shots)',
                'm1_op_idealized': '<psi0| J (H - E0) J |psi0>   (OLD demo: independent branch on FULL GS)',
            },
        },
        'm1_op_idealized': m1_op, 'E0_full_sector': E0, 'n_support': nsup,
        'first_crossing_5pct_ideal_d': cx_ideal, 'first_crossing_5pct_ideal_cov':
            (cx_ideal / nsup if cx_ideal else None),
        'first_crossing_5pct_shared_d': cx_shared, 'first_crossing_5pct_shared_cov':
            (cx_shared / nsup if cx_shared else None),
        'worst_rel_shared': float(rel_shared.max()), 'worst_rel_ideal': float(rel_ideal.max()),
        'sweep': rows,
        'verdict': (
            'The shared-state leak Delta_1^shared (classical post-processing of |0_S>, no shots) is a nonzero '
            'operator leak that GROWS as the determinant subspace is truncated: the qualitative screen still fires. '
            f'It crosses the 5% falsification line at d~{cx_shared:.0f} '
            f'(cov~{100*cx_shared/nsup:.1f}%), versus the idealized full-GS curve at d~{cx_ideal:.0f} '
            f'(cov~{100*cx_ideal/nsup:.1f}%). The shared-state crossing, not the idealized one, is the relevant '
            'reference when both branches come from one prepared state (5% relative guide line; no shot noise).'
            if (cx_shared and cx_ideal) else
            'See sweep; report crossings directly.'),
    }
    os.makedirs(RES, exist_ok=True)
    outpath = os.path.join(RES, f'{DATE}_teeth_shared_state.json')
    with open(outpath, 'w') as f:
        json.dump(out, f, indent=2)
    log(f"[shared] wrote {outpath}")

    print("\n=== TEETH, SHARED-STATE (classical, no shots) — L=12 doped ring, U/t=8 ===")
    print(f"idealized independent m1_op (on full psi0) = {m1_op:.4f}   E0={E0:.4f}   n_support={nsup}")
    print(f"{'d':>7} {'cov%':>6} {'m_hat':>9} {'m_bar':>9} {'D_sh_abs':>9} {'rel_sh%':>8} {'rel_id%':>8} {'leak%':>7}")
    for r in rows:
        print(f"{r['d']:>7} {100*r['cov']:>6.1f} {r['m_hat']:>9.4f} {r['m_bar']:>9.4f} "
              f"{r['delta_shared_abs']:>9.4f} {100*r['rel_shared']:>8.2f} {100*r['rel_ideal']:>8.2f} "
              f"{100*r['leak_frac']:>7.2f}")
    print(f"\n5% crossing: idealized d~{cx_ideal}  shared d~{cx_shared}")
    print("VERDICT:", out['verdict'])


if __name__ == '__main__':
    main()
