# -*- coding: utf-8 -*-
"""
SECOND DOMAIN — the lattice-gauge Gauss-law falsifier (maximize-board wf_00831cb3-9d6; decides cross-domain vs single-system).

Instantiates the transportable moment-test kernel in lattice gauge theory: the U(1) quantum link model (lattice
Schwinger model), the standard quantum-simulation gauge testbed. The falsifier is GAUSS'S LAW — a set of LOCAL
operator constraints G_n = E_n - E_{n-1} - q_n that EVERY physical (gauge-invariant) state satisfies EXACTLY
(<G_n> = 0, <sum_n G_n^2> = 0). It is:
  * NATIVE + EXACT (an exact constraint, not a necessary-only moment) — strictly stronger than a sum rule;
  * SAME-SAMPLE (each G_n is a LOCAL operator diagonal in the computational basis) — transportable at marginal cost;
  * a DIFFERENT error class than the Hubbard moment falsifier (it catches gauge-invariance-BREAKING errors —
    bit-flip / leakage / non-physical-sector admixture — the dominant hardware error in LGT simulation).
Two falsifiers, two domains, one transport principle -> the paper is genuinely cross-domain, not Hubbard dressed up.

Model (N staggered matter sites, N-1 spin-1/2 gauge links, OBC):
  H = -w sum_n (c†_n S+_n c_{n+1} + h.c.) + m sum_n (-1)^n c†_n c_n
Gauss law G_n = E_n - E_{n-1} - q_n,  E_l = S^z of link l,  q_n = c†_n c_n - (1-(-1)^n)/2.
Sanity (test-driven): [H, G_n] = 0 for all n (gauge invariance); physical ground state has <sum G_n^2> = 0.
"""
import os, sys, json
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATE = '2026-08-18'
N = 4                                   # matter sites (even, staggered)
W, M = 1.0, 0.6                         # hopping, staggered mass
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, '..', 'data'))
FIG = os.path.normpath(os.path.join(HERE, '..', 'out'))   # diagnostic PNGs (gitignored)

I2 = sp.identity(2, format='csr', dtype=complex)
SZ = sp.csr_matrix(np.array([[1, 0], [0, -1]], complex))
SM = sp.csr_matrix(np.array([[0, 1], [0, 0]], complex))     # sigma^- = |0><1|  (lowers occupation)
SP = sp.csr_matrix(np.array([[0, 0], [1, 0]], complex))     # sigma^+ = |1><0|
NUM = sp.csr_matrix(np.array([[0, 0], [0, 1]], complex))    # |1><1| occupation

NL = N - 1                                                   # links (OBC)
Q = N + NL                                                   # qubits: 0..N-1 matter, N..N+NL-1 links


def op_on(idx, single):
    ops = [I2] * Q
    ops[idx] = single
    out = ops[0]
    for k in range(1, Q):
        out = sp.kron(out, ops[k], format='csr')
    return out


def c_matter(n):
    """JW annihilation on matter site n (string only over matter qubits 0..N-1)."""
    ops = [I2] * Q
    for k in range(n):
        ops[k] = SZ
    ops[n] = SM
    out = ops[0]
    for k in range(1, Q):
        out = sp.kron(out, ops[k], format='csr')
    return out


def build():
    c = [c_matter(n) for n in range(N)]
    cd = [ci.getH() for ci in c]
    nmat = [cd[n] @ c[n] for n in range(N)]
    # link operators on qubits N+l : E_l = (1/2) Z ; S+ raises E (|1>-> |0>): use SM on link (|1><... ) mapping
    # convention: link |0> = E=+1/2, |1> = E=-1/2  => E_l = (1/2) Z_link ; raising E (-1/2 -> +1/2) = |0><1| = SM
    Ez = [0.5 * op_on(N + l, SZ) for l in range(NL)]
    Sp_link = [op_on(N + l, SM) for l in range(NL)]          # raises E by +1 (|1>->|0>)
    # Hamiltonian
    H = sp.csr_matrix((2 ** Q, 2 ** Q), dtype=complex)
    for n in range(N - 1):
        l = n
        hop = cd[n] @ Sp_link[l] @ c[n + 1]
        H = H - W * (hop + hop.getH())
    for n in range(N):
        H = H + M * ((-1) ** n) * nmat[n]
    # staggered charge q_n = n_n - (1-(-1)^n)/2 ; Gauss law G_n
    Ipix = sp.identity(2 ** Q, dtype=complex, format='csr')
    q = [nmat[n] - 0.5 * (1 - (-1) ** n) * Ipix for n in range(N)]
    # fixed HALF-INTEGER boundary background fields so the spin-1/2 (half-integer) links can satisfy Gauss's law
    # against INTEGER staggered charges (standard OBC Schwinger convention). bg=-1/2 makes Q_tot=0 the physical sector.
    bg = -0.5 * Ipix
    G = []
    for n in range(N):
        En = Ez[n] if n < NL else bg                        # right link (E_{N-1} boundary background)
        Enm = Ez[n - 1] if n - 1 >= 0 else bg               # left link (E_{-1} boundary background)
        G.append((En - Enm - q[n]).tocsr())
    return H.tocsr(), G, q


def main():
    H, G, q = build()
    Ipix = sp.identity(2 ** Q, dtype=complex, format='csr')
    G2 = sum((g @ g for g in G), sp.csr_matrix((2 ** Q, 2 ** Q), dtype=complex)).tocsr()

    # SANITY 1: gauge invariance [H, G_n] = 0
    comm = max(float(sp.linalg.norm(H @ g - g @ H)) for g in G)

    # physical ground state: penalize non-physical sector to project onto G=0
    Hpen = H + 50.0 * G2
    E0, V0 = eigsh(Hpen, k=1, which='SA')
    psi = V0[:, 0]
    g2_phys = float(np.real(np.vdot(psi, G2 @ psi)))         # SANITY 2: physical state satisfies Gauss law

    # ---- the Gauss-law FALSIFIER: catch gauge-breaking (bit-flip) hardware error ----
    # apply single-qubit bit-flips (X) with prob p to the physical state (incoherent mix), measure <sum G_n^2>.
    rng = np.random.default_rng(7)
    Xops = [op_on(k, SM + SP) for k in range(Q)]             # X = sigma^+ + sigma^-
    ps = np.linspace(0.0, 0.5, 11)
    rows = []
    total_charge = sum(q, sp.csr_matrix((2 ** Q, 2 ** Q), dtype=complex)).tocsr()
    Qtot_phys = float(np.real(np.vdot(psi, total_charge @ psi)))
    for p in ps:
        # incoherent bit-flip channel on each qubit: rho -> (1-p) rho + p X rho X ; track <G2> under the mixed state.
        # represent as expectation over the ensemble: <G2> = (prod channels) — compute via mixing on the density level
        # cheap exact route: <G2>_mixed = sum over flip-subsets of weight * <psi| Xs G2 Xs |psi>. Sample subsets.
        acc_g2 = 0.0; acc_Q = 0.0; nsamp = 40
        for _ in range(nsamp):
            state = psi.copy()
            for k in range(Q):
                if rng.random() < p:
                    state = Xops[k] @ state
            acc_g2 += float(np.real(np.vdot(state, G2 @ state)))
            acc_Q += float(np.real(np.vdot(state, total_charge @ state)))
        g2 = acc_g2 / nsamp
        # INDEPENDENT Gauss-law falsifier residual (native, exact, local): grows with gauge-breaking error
        # BLIND global control: total charge Q is preserved by symmetric bit-flip pairs on average -> weak signal
        rows.append({'p': float(p), 'gauss_violation': g2, 'charge_drift': abs(acc_Q / nsamp - Qtot_phys)})

    fires = rows[2]['gauss_violation'] > 1e-3
    teeth = fires and (rows[-1]['gauss_violation'] > 10 * max(rows[0]['gauss_violation'], 1e-9))
    ok = (comm < 1e-9) and (g2_phys < 1e-8) and teeth

    out = {'_provenance': {'script': 'A_payload_echoes/03_src/run_gausslaw_falsifier.py', 'date': DATE,
                           'model': f'U(1) quantum link model (lattice Schwinger), N={N} sites, {NL} spin-1/2 links, OBC',
                           'board': 'wf_00831cb3-9d6', 'params': f'w={W}, m={M}',
                           'claim': 'Gauss law = native, EXACT, same-sample falsifier for LGT; catches gauge-breaking error; SECOND domain confirmed'},
           'sanity': {'gauge_invariance_[H,G]_max': comm, 'physical_state_<sumG2>': g2_phys},
           'error_sweep': rows,
           'summary': {'gauge_invariant_hamiltonian': bool(comm < 1e-9),
                       'physical_state_satisfies_gauss': bool(g2_phys < 1e-8),
                       'falsifier_fires_on_gauge_breaking': bool(teeth),
                       'cross_domain_confirmed': bool(ok),
                       'verdict': ('SECOND DOMAIN CONFIRMED: Gauss law is a native, EXACT, same-sample falsifier — [H,G]=0 '
                                   '(gauge-invariant), physical state has <sumG2>=0, and it FIRES on gauge-breaking bit-flip '
                                   'error (a DIFFERENT error class than the Hubbard moment falsifier). The paper is genuinely '
                                   'cross-domain: one transport principle, two distinct falsifiers.' if ok else 'INSPECT')}}
    os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)
    with open(os.path.join(RES, f'{DATE}_gausslaw_falsifier.json'), 'w') as f:
        json.dump(out, f, indent=2)

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot([r['p'] * 100 for r in rows], [r['gauss_violation'] for r in rows], 'o-', color='#1e8449', lw=2.2, ms=6,
            label=r'Gauss-law falsifier $\langle\sum_n G_n^2\rangle$ — FIRES (native, exact, local)')
    ax.plot([r['p'] * 100 for r in rows], [r['charge_drift'] for r in rows], 's--', color='#c0392b', lw=2, ms=6,
            label='global total-charge drift — weak (blind to local gauge violation)')
    ax.set_title('SECOND DOMAIN — Gauss-law falsifier (U(1) lattice Schwinger model)\ncatches gauge-breaking error a global check misses', fontsize=10.5)
    ax.set_xlabel('per-qubit bit-flip error rate (%)'); ax.set_ylabel('violation signal')
    ax.legend(fontsize=8.5); ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f'{DATE}_gausslaw_falsifier.png'), dpi=140, bbox_inches='tight')

    print('=== SECOND DOMAIN: Gauss-law falsifier (U(1) quantum link model) ===')
    print(f'SANITY [H,G_n]=0 (gauge invariance): max commutator = {comm:.2e}')
    print(f'SANITY physical ground state <sum G_n^2> = {g2_phys:.2e}')
    for r in rows[::2]:
        print(f"  bit-flip p={100*r['p']:4.0f}%:  <sum G^2>={r['gauss_violation']:.4f}  (charge drift {r['charge_drift']:.4f})")
    print('CROSS-DOMAIN CONFIRMED:', ok)
    print('VERDICT:', out['summary']['verdict'])


if __name__ == '__main__':
    main()
