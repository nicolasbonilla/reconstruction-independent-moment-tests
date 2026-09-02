# -*- coding: utf-8 -*-
"""
run_collective_sumrule_falsifier.py -- DEMONSTRATED breadth of the moment screen on the COLLECTIVE channels.

Referee finding closed here: the paper claims the reconstruction-independent moment screen "transports across
observables", but only RUNS it on the optical current A_J(omega) (run_sumrule_falsifier.py / fig_falsifier) and on an
injected flatten of A(k,omega). The collective structure factors displayed in fig_sqw --
  charge S(q,omega)      via rho_q = sum_j e^{-i q j} n_j       (n_j = n_{j up} + n_{j dn})
  spin   S^zz(q,omega)   via Sz_q  = sum_j e^{-i q j} Sz_j      (Sz_j = 1/2 (n_{j up} - n_{j dn}))
are NEVER screened. This script screens them, exactly as for the current probe.

For each channel the two independent spectral moments are GROUND-STATE OPERATOR EXPECTATIONS (no reconstruction read
back), for probe A in {rho_q, Sz_q} defining the T=0 Lehmann measure S(q,omega)=sum_n |<n|A|0>|^2 delta(omega-omega_n):
    m0(q) = integral S(q,omega) domega = <0| A^dag A |0>_c                (static structure factor / f-sum zeroth moment)
    m1(q) = integral omega S(q,omega) domega = <A0|(H-E0)|A0>            (first moment; A0 = A|0> - <A>|0>)
          = (1/2) <0| [A^dag,[H,A]] |0>                                  (double-commutator f-sum form)
The double-commutator equals the one-sided <A0|(H-E0)|A0> because inversion symmetry gives S(q,.)=S(-q,.) so the +q and
-q channels carry equal first moment; we verify this numerically. Both moments are INDEPENDENT of the reconstructed
spectrum -- the whole point of the screen. Only the KINETIC term contributes to m1 (the Hubbard U-term is diagonal in the
occupation basis, so it commutes with rho_q and Sz_q); we do NOT assume this, we compute the full double commutator.

Cross-check: exact Lehmann spectrum via dense diagonalization (L=6 -> 4096-dim) reproduces both operator moments to
~1e-12 for BOTH channels at every representative q. Then an INJECTED corruption that preserves m0 (a spurious low-omega
peak stealing weight from the true band -- an in-gap phantom for the gapped charge channel, a near-zero phantom for the
gapless spin channel) is shown to PASS the zeroth-moment (weight) check yet FAIL the first-moment check, firing the
screen -- exactly as the current probe does in fig_falsifier. Honest scope: this is an EXACT-SIMULATION demonstration
(the same footing as the current-probe screen); it certifies the screen transports to the collective observables, not a
hardware result.

Params match fig_sqw physics: 1D Hubbard ring (PBC), U/t=8, half filling. Uses the L=6 dense JW engine hubbard_ed.py so
the Lehmann cross-check is exact; the L=12 fig_sqw is the same physics one system size up.
"""
import os, sys, json
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.normpath(os.path.join(HERE, '..', '..', '00_shared', 'lib'))
sys.path.insert(0, LIB); sys.path.insert(0, HERE)
import hubbard_ed as H

DATE = '2026-08-26'
L, U, T = 6, 8.0, 1.0
RES = os.path.normpath(os.path.join(HERE, '..', '06_results'))
FIGS = os.path.normpath(os.path.join(HERE, '..', 'paper', 'arxiv-submission', 'figs'))

# representative momenta q = 2*pi*m/L ; m in {1,2,3} -> q/pi in {1/3, 2/3, 1}
M_REP = [1, 2, 3]
# corruption: a spurious low-omega phantom peak that steals weight from the true band, preserving m0 exactly.
OM_SPUR = {'charge': 2.0, 'spin': 0.10}     # in-gap phantom (charge), near-zero phantom (spin)
OM_SPLIT = {'charge': U / 2.0, 'spin': 0.30}  # boundary above which "true band" weight is stolen
F_SWEEP = np.linspace(0.0, 0.5, 11)
F_REPORT = 0.4
FIRE_THRESH = 0.05      # 5% relative first-moment residual = the screen fires (matches fig_falsifier threshold)


def site_number_ops(L, c):
    """n_j (total), and Sz_j = 1/2 (n_up - n_dn), per site j, as sparse operators."""
    cd = [ci.getH() for ci in c]
    n_tot, sz = [], []
    for j in range(L):
        nup = cd[2 * j] @ c[2 * j]
        ndn = cd[2 * j + 1] @ c[2 * j + 1]
        n_tot.append((nup + ndn).tocsr())
        sz.append((0.5 * (nup - ndn)).tocsr())
    return n_tot, sz


def probe_operator(site_ops, q, L):
    """A_q = sum_j e^{-i q j} O_j  (rho_q for O=n_j; Sz_q for O=Sz_j)."""
    nq2 = site_ops[0].shape[0]
    A = sp.csr_matrix((nq2, nq2), dtype=complex)
    for j in range(L):
        A = A + np.exp(-1j * q * j) * site_ops[j]
    return A.tocsr()


def moments_operator(A, Ham, psi0, E0):
    """Independent ground-state operator moments of the Lehmann measure of probe A:
       m0 = <A0|A0>,  m1 = <A0|(H-E0)|A0>  with connected A0 = A|0> - <A>|0>. Also the double-commutator
       dc = (1/2)<0|[A^dag,[H,A]]|0> as an INDEPENDENT operator identity (must equal m1)."""
    Adag = A.getH()
    Apsi = A @ psi0
    expA = np.vdot(psi0, Apsi)                       # <A> (=0 for q!=0 by translation symmetry)
    A0 = Apsi - expA * psi0                          # connected
    m0 = float(np.real(np.vdot(A0, A0)))
    HA0 = Ham @ A0
    m1 = float(np.real(np.vdot(A0, HA0)) - E0 * m0)
    # double commutator (1/2)<[A^dag,[H,A]]> as a pure operator expectation, no A0 shortcut
    HA = Ham @ Apsi
    comm = HA - A @ (Ham @ psi0)                     # [H,A]|0>
    # [A^dag,[H,A]]|0> = A^dag [H,A]|0> - [H,A] A^dag|0>
    term1 = Adag @ comm
    Adpsi = Adag @ psi0
    term2 = (Ham @ (A @ Adpsi)) - (A @ (Ham @ Adpsi))   # [H,A] A^dag |0>
    dc = 0.5 * float(np.real(np.vdot(psi0, term1 - term2)))
    return m0, m1, dc, complex(expA)


def lehmann_moments(A, Evals, Evecs, psi0, E0):
    """Exact spectrum cross-check: S(q,omega)=sum_n |<n|A|0>|^2 delta(omega-(E_n-E0)); return (m0,m1,poles)."""
    Apsi = A @ psi0
    expA = np.vdot(psi0, Apsi)
    A0 = Apsi - expA * psi0
    overlaps = Evecs.conj().T @ A0                   # <n|A0|0>
    w = np.abs(overlaps) ** 2
    om = Evals - E0
    keep = om > 1e-8                                 # regular part, drop the (empty) elastic line
    om, w = om[keep], w[keep]
    m0 = float(w.sum()); m1 = float((om * w).sum())
    return m0, m1, om, w


def corrupt_sweep(om, w, m0_op, m1_op, om_split, om_spur):
    """Inject a spurious phantom peak at om_spur that steals fraction f of the weight above om_split, preserving m0
    EXACTLY (weight conserved) but shifting m1. Return per-f (m0_residual, m1_residual)."""
    hi = om >= om_split
    W_hi = float(w[hi].sum())
    rows = []
    for f in F_SWEEP:
        w_w = w.copy(); w_w[hi] *= (1.0 - f)
        om_ext = np.append(om, om_spur)
        w_ext = np.append(w_w, f * W_hi)             # phantom carries the stolen weight
        m0_w = float(w_ext.sum()); m1_w = float((om_ext * w_ext).sum())
        rows.append({'f': float(f),
                     'm0_residual': abs(m0_w - m0_op) / m0_op,
                     'm1_residual': abs(m1_w - m1_op) / m1_op})
    return rows, W_hi


def main():
    c = H.build_operators(L)
    cd = [ci.getH() for ci in c]
    Ham = H.hubbard_hamiltonian(L, t=T, U=U, pbc=True, c=c)
    nq = 2 * L
    Ntot = sum((cd[m] @ c[m] for m in range(nq)), sp.csr_matrix((2 ** nq, 2 ** nq), dtype=complex))

    # half-filled ground state: PH-symmetric point mu=U/2 pins N=L; moments use the TRUE H (the -mu*N shift cancels
    # in E_n-E0 because rho_q and Sz_q are N-conserving, so all reached states share N=L).
    mu = U / 2.0
    evals, evecs = eigsh((Ham - mu * Ntot).tocsc(), k=1, which='SA')
    psi0 = evecs[:, 0]
    Nfill = float(np.real(np.vdot(psi0, Ntot @ psi0)))
    E0 = float(np.real(np.vdot(psi0, Ham @ psi0)))
    var_H = float(np.real(np.vdot(psi0, Ham @ (Ham @ psi0))) - E0 ** 2)
    assert abs(Nfill - L) < 1e-8, f'not half-filled: N={Nfill}'
    assert var_H < 1e-8, f'psi0 not a clean H eigenstate: var={var_H}'

    # exact dense spectrum for the Lehmann cross-check (4096-dim)
    Hd = Ham.toarray()
    Evals, Evecs = np.linalg.eigh(Hd)

    n_tot, sz = site_number_ops(L, c)
    ops = {'charge': n_tot, 'spin': sz}
    labels = {'charge': r'S(q,\omega) via rho_q', 'spin': r'S^{zz}(q,\omega) via Sz_q'}

    out = {'_provenance': {'script': 'A_payload_echoes/03_src/run_collective_sumrule_falsifier.py', 'date': DATE,
                           'params': f'Hubbard ring L={L}, U={U}, t={T}, PBC, half filling (N={L}), mu=U/2',
                           'closes': 'referee finding: assert->demonstrate that the moment screen transports across '
                                     'observables, now RUN on collective S(q,w) [charge] and S^zz(q,w) [spin]',
                           'scope': 'EXACT-SIMULATION demonstration (same footing as the current-probe screen); '
                                    'certifies transport of the screen to the collective channels, not a hardware result'},
              'ground_state': {'E0': E0, 'N': Nfill, 'H_variance': var_H},
              'channels': {}}

    worst_m0_match = 0.0; worst_m1_match = 0.0
    all_fire = True
    print(f'=== COLLECTIVE-CHANNEL MOMENT SCREEN  (Hubbard ring L={L}, U={U}, PBC, half filling) ===')
    print(f'ground state: E0={E0:.6f}, N={Nfill:.4f}, H-variance={var_H:.1e}\n')

    for chan in ('charge', 'spin'):
        chan_rec = {'probe': labels[chan], 'om_split': OM_SPLIT[chan], 'om_spur': OM_SPUR[chan], 'per_q': {}}
        print(f'--- {chan.upper()} channel  [{labels[chan]}] ---')
        for m in M_REP:
            q = 2 * np.pi * m / L
            A = probe_operator(ops[chan], q, L)
            m0_op, m1_op, dc, expA = moments_operator(A, Ham, psi0, E0)
            m0_sp, m1_sp, om, w = lehmann_moments(A, Evals, Evecs, psi0, E0)
            r0 = abs(m0_sp - m0_op) / abs(m0_op)
            r1 = abs(m1_sp - m1_op) / abs(m1_op)
            rdc = abs(dc - m1_op) / abs(m1_op)
            worst_m0_match = max(worst_m0_match, r0); worst_m1_match = max(worst_m1_match, max(r1, rdc))
            mean_om = m1_op / m0_op

            sweep, W_hi = corrupt_sweep(om, w, m0_op, m1_op, OM_SPLIT[chan], OM_SPUR[chan])
            rep = next(r for r in sweep if abs(r['f'] - F_REPORT) < 1e-9)
            fires = (rep['m0_residual'] < 1e-9) and (rep['m1_residual'] > FIRE_THRESH)
            all_fire = all_fire and fires

            chan_rec['per_q'][f'm={m}'] = {
                'q_over_pi': 2 * m / L,
                'm0_operator': m0_op, 'm1_operator_onesided': m1_op, 'm1_double_commutator': dc,
                'mean_frequency_m1_over_m0': mean_om,
                'm0_lehmann': m0_sp, 'm1_lehmann': m1_sp,
                'xcheck_m0_residual': r0, 'xcheck_m1_residual': r1, 'xcheck_dc_vs_m1_residual': rdc,
                'expA_abs': abs(expA),
                'corruption_at_f0.4': {'m0_residual': rep['m0_residual'], 'm1_residual': rep['m1_residual']},
                'fires': bool(fires),
                'sweep': sweep}
            print(f'  q/pi={2*m/L:.3f}:  m0={m0_op:.6f}  m1={m1_op:.6f}  <w>=m1/m0={mean_om:.4f}  '
                  f'dc(1/2<[A+,[H,A]]>)={dc:.6f}')
            print(f'            Lehmann xcheck: m0 res={r0:.1e}  m1 res={r1:.1e}  dc-vs-m1 res={rdc:.1e}')
            print(f'            corruption f=0.4: m0 res={rep["m0_residual"]:.1e} (PASS weight)  '
                  f'm1 res={rep["m1_residual"]*100:.2f}% ({"FIRES" if fires else "quiet"})')
        out['channels'][chan] = chan_rec
        print()

    out['summary'] = {
        'worst_lehmann_m0_match_residual': worst_m0_match,
        'worst_lehmann_m1_match_residual': worst_m1_match,
        'xcheck_passes_1e-10': bool(worst_m0_match < 1e-10 and worst_m1_match < 1e-10),
        'screen_fires_both_channels_all_q': bool(all_fire),
        'fire_threshold': FIRE_THRESH,
        'verdict': ('DEMONSTRATED: the reconstruction-independent moment screen RUNS on both collective channels -- '
                    'operator m0,m1 match the exact Lehmann spectrum to ~1e-12, and an m0-preserving corruption passes '
                    'the weight check yet fails the first moment, firing the screen for BOTH S(q,w) and S^zz(q,w).'
                    if (worst_m0_match < 1e-10 and worst_m1_match < 1e-10 and all_fire) else 'INSPECT')}

    os.makedirs(RES, exist_ok=True)
    jpath = os.path.join(RES, f'{DATE}_collective_sumrule_falsifier.json')
    with open(jpath, 'w') as f:
        json.dump(out, f, indent=2)

    # ---- figure data: residual sweep for the representative momentum of each channel (q/pi=2/3, m=2) ----
    os.makedirs(FIGS, exist_ok=True)
    m_fig = 2
    for chan in ('charge', 'spin'):
        sweep = out['channels'][chan]['per_q'][f'm={m_fig}']['sweep']
        dpath = os.path.join(FIGS, f'collective_screen_{chan}.dat')
        with open(dpath, 'w') as f:
            f.write('f r0 r1\n')
            for r in sweep:
                f.write(f'{r["f"]*100:.3f} {r["m0_residual"]*100:.6f} {r["m1_residual"]*100:.4f}\n')
        print('wrote', dpath)

    print('\nVERDICT:', out['summary']['verdict'])
    print('worst Lehmann match: m0', f'{worst_m0_match:.1e}', ' m1/dc', f'{worst_m1_match:.1e}')
    print('saved:', jpath)


if __name__ == '__main__':
    main()
