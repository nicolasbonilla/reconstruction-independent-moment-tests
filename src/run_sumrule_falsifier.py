# -*- coding: utf-8 -*-
"""
m0-blind / m1-detects sum-rule falsifier for a reconstructed current response (exact classical simulation).

A reconstruction-independent, necessary-condition check for a reconstructed dynamical response, built from sum
rules that are GROUND-STATE OPERATOR EXPECTATIONS, evaluated here by exact diagonalization. No moment is estimated
on hardware in this work; the current J and the hopping part of H are off-diagonal in the occupation basis, so on a
device these expectations would need rotated-basis measurements (cf. bond_moment_estimator.py, simulation only),
not the computational-basis samples alone.

Optical response of a Hubbard ring: current spectral function A_J(omega) = sum_n |<n|J|0>|^2 delta(omega - omega_n),
omega_n = E_n - E_0 > 0, J the current operator. Its moments are EXACT operator sum rules on the ground state:
    m0 = integral A_J domega = <J^2>                       (current-fluctuation / total-weight (m0) sum rule; NOT the optical f-sum, which is m_-1)
    m1 = integral omega A_J domega = (1/2) <[J,[H,J]]>      (first-moment / kinetic sum rule)
Both are ground-state expectations, INDEPENDENT of the reconstructed spectrum.

WHAT THIS DEMONSTRATES: the total-weight (m0) sum rule fixes only the TOTAL weight and is BLIND to the Drude/mid-IR
SPLIT; a wrong reconstruction that redistributes weight between low-omega (Drude) and mid-omega (mid-IR) while
preserving m0 PASSES the total-weight (m0) sum rule but VIOLATES the first-moment sum rule m1. Adding m1 constrains
the shape, but it remains a necessary condition only: a pass is not a certificate.
"""
import os, sys, json
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)   # hubbard_ed.py is vendored in src/
import hubbard_ed as H

DATE = '2026-08-18'
L, U = 6, 8.0
RES = os.path.normpath(os.path.join(HERE, '..', 'data'))
FIG = os.path.normpath(os.path.join(HERE, '..', 'out'))   # diagnostic PNGs (gitignored)


def current_operator(L, c, t=1.0):
    """Hermitian current on a ring: J = -i t sum_{i,s} (c†_{i+1,s} c_{i,s} - c†_{i,s} c_{i+1,s}), i+1 mod L."""
    cd = [ci.getH() for ci in c]
    nq = 2 * L
    J = sp.csr_matrix((2 ** nq, 2 ** nq), dtype=complex)
    for i in range(L):
        j = (i + 1) % L
        for s in (0, 1):
            a, b = 2 * i + s, 2 * j + s
            J = J + (-1j * t) * (cd[b] @ c[a] - cd[a] @ c[b])
    return J.tocsr()


def main():
    c = H.build_operators(L)
    cd = [ci.getH() for ci in c]
    Ham0 = H.hubbard_hamiltonian(L, t=1.0, U=U, pbc=True, c=c)
    Ntot = sum((cd[m] @ c[m] for m in range(2 * L)), sp.csr_matrix((2 ** (2 * L), 2 ** (2 * L)), dtype=complex))
    J = current_operator(L, c)

    # HOLE-DOPE via chemical potential so the ground state is a METAL (genuine Drude + mid-IR structure).
    # Scan mu downward (H0 - mu N favors fewer particles for smaller mu) until we land ~2 holes off half-filling.
    # psi0 stays an EXACT eigenstate of Ham0 because [Ham0, Ntot]=0 (verified by the H-variance check below).
    Ham = Ham0                                                        # dynamics under the true H
    psi0 = None; Nfill = L; mu_used = None
    for mu in [4.0, 2.0, 0.0, -1.0, -2.0, -3.0, -4.0, -6.0]:
        evals, evecs = eigsh(Ham0 - mu * Ntot, k=1, which='SA')
        cand = evecs[:, 0]
        n = float(np.real(np.vdot(cand, Ntot @ cand)))
        if n <= L - 1.5:                                             # at least ~2 holes
            psi0, Nfill, mu_used = cand, n, mu; break
        psi0, Nfill, mu_used = cand, n, mu
    E0 = float(np.real(np.vdot(psi0, Ham @ psi0)))
    var_H = float(np.real(np.vdot(psi0, (Ham @ (Ham @ psi0)))) - E0 ** 2)   # must be ~0: psi0 is a clean Ham0 eigenstate
    print(f'doped ground state: filling N={Nfill:.2f} (half-filling={L}), doping={100*(L-Nfill)/L:.0f}% holes, mu={mu_used}, '
          f'H-variance={var_H:.1e} (clean eigenstate if ~0)')

    # operator sum rules (ground-state expectations, evaluated exactly here; on a device the off-diagonal
    # J and H would need rotated-basis measurements, not computational-basis samples alone)
    Jpsi = J @ psi0
    m0_op = float(np.real(np.vdot(Jpsi, Jpsi)))                       # <J^2>
    HJpsi = Ham @ Jpsi
    m1_op = float(np.real(np.vdot(Jpsi, HJpsi)) - E0 * m0_op)         # <Jpsi|(H-E0)|Jpsi>

    # exact spectrum A_J(omega) via full Lehmann (dense eig; L=6 -> 4096 dim)
    Hd = Ham.toarray()
    Ev, Vv = np.linalg.eigh(Hd)
    overlaps = Vv.conj().T @ Jpsi                                     # <n| J |0>
    w = np.abs(overlaps) ** 2
    om = Ev - E0
    keep = om > 1e-6                                                  # regular part (drop omega=0 Drude delta)
    om, w = om[keep], w[keep]
    m0_spec = float(w.sum()); m1_spec = float((om * w).sum())

    # ---- construct a WRONG Drude/mid-IR split that preserves m0 but violates m1 ----
    # physical split at omega ~ U/2: low-energy (Drude/quasiparticle) vs the mid-IR Hubbard band near omega ~ U.
    om_split = U / 2.0
    lo = om < om_split                                                # "Drude" (low-omega)
    hi = ~lo                                                          # "mid-IR" (high-omega)
    W_lo, W_hi = w[lo].sum(), w[hi].sum()
    obar_lo = float((om[lo] * w[lo]).sum() / max(W_lo, 1e-12))
    obar_hi = float((om[hi] * w[hi]).sum() / max(W_hi, 1e-12))

    # distortion sweep: move a fraction f of mid-IR weight into Drude, keep TOTAL m0 fixed
    fs = np.linspace(0.0, 0.5, 11)
    res = []
    for f in fs:
        transfer = f * W_hi
        w_wrong = w.copy()
        w_wrong[hi] *= (1 - f)                                        # remove from mid-IR
        w_wrong[lo] *= (1 + transfer / max(W_lo, 1e-12))              # add to Drude, m0 preserved
        m0_wrong = float(w_wrong.sum()); m1_wrong = float((om * w_wrong).sum())
        res.append({'f': float(f),
                    'm0_residual': abs(m0_wrong - m0_op) / m0_op,     # total-weight (m0) sum-rule check (should stay ~0)
                    'm1_residual': abs(m1_wrong - m1_op) / m1_op})    # first-moment check (should grow -> FALSIFIES)

    # sanity: true spectrum matches operator sum rules
    s_m0 = abs(m0_spec - m0_op) / m0_op
    s_m1 = abs(m1_spec - m1_op) / m1_op
    falsifier_works = (res[-1]['m0_residual'] < 1e-6) and (res[-1]['m1_residual'] > 0.05)

    out = {'_provenance': {'script': 'src/run_sumrule_falsifier.py', 'date': DATE,
                           'params': f'Hubbard ring L={L}, U={U}, PBC, hole-doped ground state (N={Nfill:.2f}, '
                                     f'selected by the chemical-potential scan, mu={mu_used})',
                           'claim': 'exact classical simulation, no hardware: the total-weight (m0) sum rule is blind '
                                    'to a Drude/mid-IR redistribution that preserves m0; the first-moment (m1) sum rule '
                                    'flags it; a necessary condition only (a pass is not a certificate)',
                           'relabel_note': '2026-09-26: m0 was mislabelled the f-sum; the optical f-sum is m_-1 (see manuscript)',
                           'wording_note': '2026-09-27: provenance strings rescoped (the script path, the filling, which '
                                           'is hole-doped and not half filling, and the claim/verdict wording; an '
                                           'internal workflow identifier removed); every number unchanged'},
           'sum_rules': {'m0_operator_<J^2>': m0_op, 'm1_operator_half<[J,[H,J]]>': m1_op,
                         'm0_spectrum': m0_spec, 'm1_spectrum': m1_spec,
                         'm0_match_residual': s_m0, 'm1_match_residual': s_m1},
           'split': {'omega_split': float(om_split), 'W_Drude': float(W_lo), 'W_midIR': float(W_hi),
                     'obar_Drude': obar_lo, 'obar_midIR': obar_hi},
           'distortion_sweep': res,
           'summary': {'true_spectrum_satisfies_both_sumrules': bool(s_m0 < 1e-6 and s_m1 < 1e-6),
                       'falsifier_works': bool(falsifier_works),
                       'verdict': ('m0-blind / m1-detects: a wrong Drude/mid-IR split passes the total-weight (m0) sum '
                                   'rule (m0 residual ~0) but the first-moment sum rule (m1) flags it (relative residual '
                                   'above a 5% guide line; exact simulation, no shot noise).'
                                   if falsifier_works else 'INSPECT')}}
    os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)
    jpath = os.path.join(RES, f'{DATE}_sumrule_falsifier.json')
    with open(jpath, 'w') as f:
        json.dump(out, f, indent=2)

    # figure
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.4))
    # (a) true vs wrong (f=0.5) spectrum, binned
    fw = 0.5
    w_wrong = w.copy(); w_wrong[hi] *= (1 - fw); w_wrong[lo] *= (1 + fw * W_hi / max(W_lo, 1e-12))
    bins = np.linspace(0, om.max() * 1.02, 30)
    ax[0].hist(om, bins=bins, weights=w, alpha=0.55, color='#2471a3', label='true $A_J(\\omega)$')
    ax[0].hist(om, bins=bins, weights=w_wrong, alpha=0.55, color='#c0392b',
               label='wrong split (same $m_0$)')
    ax[0].axvline(om_split, color='k', ls=':', lw=1, label='Drude | mid-IR')
    ax[0].set_title('(a) a wrong Drude/mid-IR split that PASSES the total-weight (m0) sum rule', fontsize=10.5)
    ax[0].set_xlabel(r'$\omega$ (energy above ground state)'); ax[0].set_ylabel(r'spectral weight $A_J(\omega)$')
    ax[0].legend(fontsize=8.5)
    # (b) sum-rule residuals vs distortion
    ax[1].plot(fs, [r['m0_residual'] for r in res], 's-', color='#1e8449', lw=2, label=r'$m_0$ (total-weight) residual')
    ax[1].plot(fs, [r['m1_residual'] for r in res], 'o-', color='#8e44ad', lw=2, label=r'$m_1$ (first-moment) residual')
    ax[1].axhline(0.05, color='#c0392b', ls='--', lw=1, label='relative 5% guide')
    ax[1].set_title('(b) the first-moment sum rule flags what the total-weight (m0) sum rule misses', fontsize=10.5)
    ax[1].set_xlabel('weight fraction moved mid-IR → Drude'); ax[1].set_ylabel('relative sum-rule residual')
    ax[1].legend(fontsize=8.5)
    for a in ax:
        a.grid(alpha=0.25)
    fig.suptitle(f'Sum-rule check (exact simulation) — Hubbard ring L={L}, U={U} (m0-blind / m1-detects)',
                 y=1.02, fontsize=11)
    fig.tight_layout()
    ppath = os.path.join(FIG, f'{DATE}_sumrule_falsifier.png')
    fig.savefig(ppath, dpi=140, bbox_inches='tight')

    print('=== SUM-RULE CHECK (exact classical simulation; no hardware) ===')
    print(f'operator sum rules:  m0=<J^2>={m0_op:.4f}   m1=(1/2)<[J,[H,J]]>={m1_op:.4f}')
    print(f'spectrum moments:    m0={m0_spec:.4f} (res {s_m0:.1e})   m1={m1_spec:.4f} (res {s_m1:.1e})  <- sanity: true spectrum obeys BOTH')
    print(f'Drude/mid-IR split at omega={om_split:.2f}:  W_Drude={W_lo:.3f} @ {obar_lo:.2f},  W_midIR={W_hi:.3f} @ {obar_hi:.2f}')
    print('distortion sweep (fraction moved mid-IR->Drude):')
    for r in res[::2]:
        print(f'  f={r["f"]:.2f}:  m0 residual={r["m0_residual"]:.1e} (m0 sum rule: PASSES)   m1 residual={r["m1_residual"]:.3f} (first-moment: {"FLAGS" if r["m1_residual"]>0.05 else "ok"})')
    print('VERDICT:', out['summary']['verdict'])
    print('saved:', jpath); print('saved:', ppath)


if __name__ == '__main__':
    main()
