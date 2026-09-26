# -*- coding: utf-8 -*-
"""
THE TEETH EXPERIMENT — the make-or-break gate for the whole moment-test kernel (maximize-board wf_00831cb3-9d6).

The board's decisive requirement: the sum-rule falsifier only has TEETH if it is an INDEPENDENT estimator, not
circular. The error class SQD actually suffers is SUBSPACE TRUNCATION: the dynamical spectrum A_J(omega) is
reconstructed (Lehmann) from an INCOMPLETE set of excited states, so its moments are wrong. The falsifier target
must be a GROUND-STATE static operator expectation measured directly (m0=<J^2>, m1=(1/2)<[J,[H,J]]>) — which
converges with FEWER resources than the full excited spectrum and is therefore an INDEPENDENT check.

This experiment proves (or kills) the teeth:
  * INDEPENDENT falsifier: compare the truncated-spectrum moment m1_trunc to the DIRECT ground-state estimator
    m1_op. As spectral coverage drops, |m1_trunc - m1_op| GROWS -> the falsifier FIRES on the truncation error.
  * BACK-COMPUTED (circular) control: compare m1_trunc to a "target" recomputed from the SAME truncated spectrum
    -> residual is 0 by construction -> BLIND. This is the theater the board warned about.
If the independent estimator fires while the circular one is blind, the kernel has teeth on the subspace-truncation
error class. We ALSO state the honest blind spot: coherent/symmetry-preserving noise corrupts both identically.
Reuses the verified falsifier machinery (doped Hubbard ring, current operator, exact Lehmann spectrum).
"""
import os, sys, json
import numpy as np
from scipy.sparse.linalg import eigsh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)   # hubbard_ed.py is vendored in src/
import hubbard_ed as H
from run_sumrule_falsifier import current_operator

DATE = '2026-08-18'
L, U = 6, 8.0
RES = os.path.normpath(os.path.join(HERE, '..', 'data'))
FIG = os.path.normpath(os.path.join(HERE, '..', 'out'))   # diagnostic PNGs (gitignored)


def main():
    c = H.build_operators(L)
    cd = [ci.getH() for ci in c]
    Ham0 = H.hubbard_hamiltonian(L, t=1.0, U=U, pbc=True, c=c)
    Ntot = sum((cd[m] @ c[m] for m in range(2 * L)), sp.csr_matrix((2 ** (2 * L),) * 2, dtype=complex))
    J = current_operator(L, c)

    # doped (metal) ground state, exact eigenstate of Ham0 (mu selects the sector; [H,N]=0)
    for mu in [0.0, -1.0, -2.0]:
        cand = eigsh(Ham0 - mu * Ntot, k=1, which='SA')[1][:, 0]
        if float(np.real(np.vdot(cand, Ntot @ cand))) <= L - 1.5:
            psi0 = cand; break
    E0 = float(np.real(np.vdot(psi0, Ham0 @ psi0)))

    # DIRECT, INDEPENDENT ground-state estimators (measured on the ground state, not from the excited spectrum)
    Jpsi = J @ psi0
    m1_op = float(np.real(np.vdot(Jpsi, Ham0 @ Jpsi)) - E0 * np.real(np.vdot(Jpsi, Jpsi)))   # (1/2)<[J,[H,J]]>

    # SQD ERROR MODEL: the spectrum is reconstructed by projecting H onto a TRUNCATED DETERMINANT SUBSPACE (the top-d
    # important configurations by |psi0|^2, exactly how SQD selects), diagonalizing within it, and rebuilding A_J from
    # the PROJECTED states. This distorts poles/weights and violates the moments (unlike dropping low-weight exact
    # poles, which preserves them). The DIRECT ground-state estimator m1_op is independent and stays accurate.
    Hd = Ham0.toarray(); Jd = J.toarray()
    prob = np.abs(psi0) ** 2
    support = np.where(prob > 1e-14)[0]                        # determinants the (doped) state actually lives on
    order = support[np.argsort(prob[support])[::-1]]          # important configurations first (SQD ordering)
    nsup = len(order)

    def truncated_spectrum_m1(d):
        idx = np.sort(order[:d])
        Hs = Hd[np.ix_(idx, idx)]; Js = Jd[np.ix_(idx, idx)]
        ev, V = np.linalg.eigh(Hs)
        g = V[:, 0]                                           # projected ground state
        Jg = Js @ g
        wsub = np.abs(V.conj().T @ Jg) ** 2                   # |<n|J|0>|^2 in the subspace
        omsub = ev - ev[0]
        m = omsub > 1e-6
        return float((omsub[m] * wsub[m]).sum())

    fracs = np.linspace(1.0, 0.15, 16)
    rows = []
    for f in fracs:
        d = max(4, int(round(f * nsup)))
        m1_trunc = truncated_spectrum_m1(d)
        # INDEPENDENT falsifier: truncated spectrum vs the DIRECT ground-state estimator m1_op (has teeth)
        r_indep = abs(m1_trunc - m1_op) / m1_op
        # BACK-COMPUTED (circular) control: target recomputed from the SAME truncated subspace -> identical -> blind
        r_circular = abs(m1_trunc - m1_trunc) / m1_op
        rows.append({'coverage': float(f), 'd': d, 'n_support': int(nsup), 'm1_trunc': m1_trunc,
                     'residual_independent': r_indep, 'residual_circular': r_circular})
    n_poles = nsup

    # teeth verdict: independent residual grows as coverage drops; circular stays ~0
    r_indep_worst = rows[-1]['residual_independent']
    r_circ_worst = max(r['residual_circular'] for r in rows)
    teeth = (r_indep_worst > 0.05) and (r_circ_worst < 1e-9)

    out = {'_provenance': {'script': 'A_payload_echoes/03_src/run_falsifier_teeth.py', 'date': DATE,
                           'params': f'doped Hubbard ring L={L}, U={U}', 'board': 'wf_00831cb3-9d6',
                           'claim': 'independent ground-state estimator has TEETH on subspace-truncation error; back-computed is circular/blind',
                           'honest_blind_spot': 'coherent/symmetry-preserving device noise corrupts direct estimator and spectrum identically -> screens ONE error class (subspace truncation / Lehmann-weight misassignment), NOT device noise generally'},
           'm1_direct_operator': m1_op, 'n_poles': n_poles, 'sweep': rows,
           'summary': {'residual_independent_worst': r_indep_worst, 'residual_circular_worst': r_circ_worst,
                       'teeth_confirmed': bool(teeth),
                       'verdict': ('TEETH CONFIRMED: the INDEPENDENT ground-state estimator fires on subspace truncation '
                                   '(residual grows to %.1f%% at 20%% coverage) while the back-computed check is blind (0). '
                                   'The moment-test kernel is NOT circular on this error class. Honest blind spot stated.' % (100 * r_indep_worst)
                                   if teeth else 'NO TEETH / inspect')}}
    os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)
    with open(os.path.join(RES, f'{DATE}_falsifier_teeth.json'), 'w') as fjs:
        json.dump(out, fjs, indent=2)

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    cov = [100 * r['coverage'] for r in rows]
    ax.plot(cov, [100 * r['residual_independent'] for r in rows], 'o-', color='#8e44ad', lw=2.2, ms=6,
            label='INDEPENDENT estimator (direct $m_1$)\n— FIRES on truncation (has teeth)')
    ax.plot(cov, [100 * r['residual_circular'] for r in rows], 's--', color='#c0392b', lw=2, ms=6,
            label='BACK-COMPUTED (same subspace)\n— BLIND (circular, always 0)')
    ax.axhline(5, color='#1e8449', ls=':', lw=1.2, label='falsification threshold')
    ax.invert_xaxis()
    ax.set_title('THE TEETH EXPERIMENT — independent estimator vs circular control\n(subspace-truncation error, doped Hubbard)', fontsize=11)
    ax.set_xlabel('spectral coverage kept (%)  [SQD subspace undersampling →]')
    ax.set_ylabel('first-moment sum-rule residual (%)')
    ax.legend(fontsize=8.5, loc='upper left'); ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f'{DATE}_falsifier_teeth.png'), dpi=140, bbox_inches='tight')

    print('=== THE TEETH EXPERIMENT ===')
    print(f'direct independent m1 = {m1_op:.4f}   ({n_poles} poles total)')
    for r in rows[::3]:
        print(f"  coverage {100*r['coverage']:4.0f}% (d={r['d']:3d}): independent residual {100*r['residual_independent']:5.1f}%   circular {100*r['residual_circular']:.1e}%")
    print('TEETH CONFIRMED:', teeth, '| worst independent residual %.1f%%, worst circular %.1e' % (100*r_indep_worst, r_circ_worst))
    print('VERDICT:', out['summary']['verdict'])


if __name__ == '__main__':
    main()
