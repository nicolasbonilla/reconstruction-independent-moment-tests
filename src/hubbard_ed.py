# -*- coding: utf-8 -*-
# Vendored 2026-09-26 from the author's shared ED library (numpy/scipy only); code unchanged, docstring shortened.
"""
hubbard_ed.py  —  exact-diagonalization utilities for the 1D Fermi-Hubbard chain (Jordan-Wigner, numpy/scipy only).

1D Fermi-Hubbard chain, Jordan-Wigner mapped to 2L qubits: fermionic operators, the Hamiltonian,
real-time evolution and a few exact-state diagnostics. This repository uses build_operators()
and hubbard_hamiltonian().

Mode ordering (crucial): mode m = 2*site + spin  ->  qubit m, with qubit 0 the MOST significant.
Therefore a contiguous block of the first L_A sites = qubits 0..2*L_A-1 (contiguous).

No project-specific deps: numpy + scipy only.  Verified by sanity checks in `_selftest()`.
"""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply

I2 = sp.identity(2, format='csr', dtype=complex)
SZ = sp.csr_matrix(np.array([[1, 0], [0, -1]], dtype=complex))
SM = sp.csr_matrix(np.array([[0, 1], [0, 0]], dtype=complex))   # sigma^- = |0><1| : lowers (annihilates)


def _op_on(q, op, nq):
    """Place single-qubit `op` on qubit q (qubit 0 most significant) among nq qubits."""
    ops = [I2] * nq
    ops[q] = op
    out = ops[0]
    for k in range(1, nq):
        out = sp.kron(out, ops[k], format='csr')
    return out


def jw_annihilation(i, nq):
    """Jordan-Wigner annihilation operator c_i on nq qubits: (Z_0...Z_{i-1}) sigma^-_i."""
    ops = [I2] * nq
    for k in range(i):
        ops[k] = SZ
    ops[i] = SM
    out = ops[0]
    for k in range(1, nq):
        out = sp.kron(out, ops[k], format='csr')
    return out


def build_operators(L):
    """Return c[m] for all 2L modes; mode m = 2*site+spin (spin 0=up,1=down)."""
    nq = 2 * L
    return [jw_annihilation(m, nq) for m in range(nq)]


def hubbard_hamiltonian(L, t=1.0, U=8.0, t2=0.0, pbc=False, c=None):
    """1D Hubbard: H = -t sum_{<ij>,s} - t2 sum_{<<ik>>,s} (c†c + h.c.) + U sum_i n_i↑ n_i↓.
    t2 = next-nearest-neighbour hopping: t2=0 is the Bethe-ansatz-INTEGRABLE model (relaxes to a GGE);
    t2!=0 (e.g. 0.3t) breaks integrability -> genuine thermalization/scrambling. [review-board FIX-phys#1]"""
    nq = 2 * L
    if c is None:
        c = build_operators(L)
    cd = [ci.getH() for ci in c]
    H = sp.csr_matrix((2 ** nq, 2 ** nq), dtype=complex)
    # nearest-neighbour hopping (same spin)
    nn = [(s, s + 1) for s in range(L - 1)] + ([(L - 1, 0)] if pbc and L > 2 else [])
    for (a, b) in nn:
        for spin in (0, 1):
            ma, mb = 2 * a + spin, 2 * b + spin
            H = H - t * (cd[ma] @ c[mb] + cd[mb] @ c[ma])
    # next-nearest-neighbour hopping (integrability breaker)
    if t2 != 0.0:
        nnn = [(s, s + 2) for s in range(L - 2)] + ([(L - 2, 0), (L - 1, 1)] if pbc and L > 3 else [])
        for (a, b) in nnn:
            for spin in (0, 1):
                ma, mb = 2 * a + spin, 2 * b + spin
                H = H - t2 * (cd[ma] @ c[mb] + cd[mb] @ c[ma])
    # on-site interaction
    for s in range(L):
        nup = cd[2 * s] @ c[2 * s]
        ndn = cd[2 * s + 1] @ c[2 * s + 1]
        H = H + U * (nup @ ndn)
    return H.tocsr()


def cdw_state(L):
    """Charge-density-wave product state: even sites doubly occupied, odd sites empty. Half filling."""
    nq = 2 * L
    bits = ['0'] * nq
    for s in range(0, L, 2):          # even sites doubly occupied
        bits[2 * s] = '1'             # up
        bits[2 * s + 1] = '1'         # down
    idx = int(''.join(bits), 2)       # qubit 0 = most significant
    psi = np.zeros(2 ** nq, dtype=complex)
    psi[idx] = 1.0
    return psi


def evolve(psi0, H, t):
    """|psi(t)> = exp(-i H t) |psi0>."""
    if t == 0:
        return psi0.astype(complex).copy()
    return expm_multiply(-1j * H * t, psi0.astype(complex))


def subsystem_purity(psi, L, L_A):
    """Tr[rho_A^2] for A = first L_A sites (qubits 0..2L_A-1). Pure global state assumed."""
    nq = 2 * L
    dA = 2 ** (2 * L_A)
    dB = 2 ** (nq - 2 * L_A)
    M = psi.reshape(dA, dB)           # qubit 0 most significant -> A is the leading index block
    rhoA = M @ M.conj().T             # (dA x dA)
    return float(np.real(np.trace(rhoA @ rhoA)))


def renyi2(psi, L, L_A):
    p = subsystem_purity(psi, L, L_A)
    return -np.log(max(p, 1e-300))


def one_rdm(psi, L, c=None):
    """gamma_{mn} = <psi| c†_m c_n |psi>  (2L x 2L one-body reduced density matrix)."""
    nq = 2 * L
    if c is None:
        c = build_operators(L)
    g = np.zeros((nq, nq), dtype=complex)
    for m in range(nq):
        cm_psi_d = (c[m] @ psi)       # c_m |psi>
        for n in range(nq):
            cn_psi = (c[n] @ psi)
            g[m, n] = np.vdot(cm_psi_d, cn_psi)   # <psi|c†_m c_n|psi> = <c_m psi | c_n psi>
    return g


def one_body_magic_F1(psi, L, c=None):
    """F1 = 4 tr[gamma(1-gamma)] = 4(tr gamma - tr gamma^2). Gaussian-invariant one-body magic."""
    g = one_rdm(psi, L, c)
    return float(np.real(4.0 * (np.trace(g) - np.trace(g @ g))))


def _selftest(L=4, U=8.0, t=1.0):
    import numpy.linalg as npl
    c = build_operators(L)
    H = hubbard_hamiltonian(L, t=1.0, U=U, c=c)
    psi0 = cdw_state(L)
    checks = []
    # (1) initial product state: every subsystem pure
    checks.append(('init purity L_A=2 == 1', abs(subsystem_purity(psi0, L, 2) - 1.0) < 1e-10))
    # (2) initial CDW is a Slater determinant -> F1 == 0
    checks.append(('init F1 == 0', abs(one_body_magic_F1(psi0, L, c)) < 1e-9))
    # (3) evolve; global purity stays 1 (pure state)
    psit = evolve(psi0, H, 0.8)
    checks.append(('norm preserved', abs(np.vdot(psit, psit).real - 1.0) < 1e-9))
    checks.append(('global purity == 1', abs(subsystem_purity(psit, L, L) - 1.0) < 1e-9))
    # (4) energy conservation
    e0 = np.vdot(psi0, H @ psi0).real
    et = np.vdot(psit, H @ psit).real
    checks.append(('energy conserved', abs(e0 - et) < 1e-8))
    # (5) entanglement grows away from 0 after the quench
    checks.append(('S2 grows after quench', renyi2(psit, L, 2) > 1e-3))
    # (6) U=0 stays Gaussian: F1 stays ~0 under free evolution of a Slater determinant
    H0 = hubbard_hamiltonian(L, t=1.0, U=0.0, c=c)
    psi_free = evolve(psi0, H0, 0.8)
    checks.append(('U=0 keeps F1 ~ 0 (Gaussian)', abs(one_body_magic_F1(psi_free, L, c)) < 1e-7))
    # (7) interactions generate one-body magic content is state-dependent; just check finiteness
    checks.append(('F1(t) finite & >=0 with U>0', one_body_magic_F1(psit, L, c) >= -1e-9))
    return checks


if __name__ == '__main__':
    print('hubbard_ed self-test (L=4, U=8):')
    ok = True
    for name, passed in _selftest():
        print(f'  [{"PASS" if passed else "FAIL"}] {name}')
        ok = ok and passed
    print('ALL PASS' if ok else 'SOME FAILED')
