# -*- coding: utf-8 -*-
"""RECONSTRUCTION-INDEPENDENT MOMENT ESTIMATOR FROM DEVICE-STYLE MEASUREMENTS (B1 core).

The paper's independent estimator computes spectral moments m_j of a probe's spectral function
as ground-state operator expectations -- from the SAME samples that produced the spectrum. For a
DIAGONAL probe (density rho_q) m_0=<drho^2> is a pure counts estimator; but m_1=<drho(H-E0)drho>
contains the HOPPING (off-diagonal) part of H, which is NOT counts-only. The rigorous, device-
ready way to estimate it -- exactly what a real experiment does -- is Pauli grouping: expand the
Hermitian operator into Pauli strings, group them into qubit-wise-commuting (QWC) sets, and
measure each set in one rotated basis (a few 'bond-basis' circuits). This module builds those
operators and estimators and VERIFIES, by finite-shot sampling, that the counts-based m_0, m_1
recover the exact operator moments within shot noise. On real hardware only the state prep and
the sampler backend change; the estimator is identical.

Model (device demonstration): L=6 Hubbard chain (OBC, clean Jordan-Wigner, no wrap-around sign),
2L=12 qubits, blocked ordering (qubits 0..L-1 spin-up sites, L..2L-1 spin-down). Probe: staggered
density rho_q = sum_i cos(q i)(n_{i,up}+n_{i,dn}), q=pi (diagonal). SIM-ONLY verification here.
"""
import numpy as np
from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit import QuantumCircuit
import scipy.sparse.linalg as sla

import os as _os
L = int(_os.environ.get('BME_L', '6'))    # L=6 (12q) matches the paper's hardware system; chunked
T, U = 1.0, 4.0                            # products + matrix-free exact moments make the full m0,m1,m2 tractable.
NQ = 2 * L
QPROBE = np.pi


def _pstr(letters):
    """length-NQ Pauli label in qiskit convention (qubit 0 = RIGHTMOST char)."""
    s = ['I'] * NQ
    for qb, l in letters.items():
        s[qb] = l
    return ''.join(reversed(s))


def build_hubbard():
    terms = []
    for spin in (0, 1):
        for i in range(L - 1):                      # OBC hopping (adjacent -> no JW string)
            a, b = spin * L + i, spin * L + i + 1
            terms.append((_pstr({a: 'X', b: 'X'}), -T * 0.5))
            terms.append((_pstr({a: 'Y', b: 'Y'}), -T * 0.5))
    for i in range(L):                              # U n_up n_dn ; n_p = 1/2 - 1/2 Z_p
        a, b = i, L + i
        terms += [(_pstr({}), U * 0.25), (_pstr({a: 'Z'}), -U * 0.25),
                  (_pstr({b: 'Z'}), -U * 0.25), (_pstr({a: 'Z', b: 'Z'}), U * 0.25)]
    return SparsePauliOp.from_list(terms).simplify()


def build_rho():
    terms = []
    for i in range(L):
        c = float(np.cos(QPROBE * i))
        for spin in (0, 1):
            p = spin * L + i
            terms += [(_pstr({}), c * 0.5), (_pstr({p: 'Z'}), -c * 0.5)]
    return SparsePauliOp.from_list(terms).simplify()


def pauli_expectation_from_counts(label, counts, ntot):
    """<P> from computational-basis counts already measured in P's eigenbasis. label is the
    qiskit string (qubit0 = rightmost); support = non-I positions; eigenvalue = prod (-1)^bit."""
    supp = [NQ - 1 - k for k, ch in enumerate(label) if ch != 'I']   # qubit indices with non-I
    if not supp:
        return 1.0
    acc = 0.0
    for bitstr, c in counts.items():
        par = sum(int(bitstr[NQ - 1 - j]) for j in supp) & 1        # bit of qubit j
        acc += (-c if par else c)
    return acc / ntot


_LET = {'I': 0, 'X': 1, 'Y': 2, 'Z': 3}
_LETINV = 'IXYZ'


def qwc_group_labels(labels):
    """fast greedy qubit-wise-commuting grouping (vectorized) -> (basis_strings, [label lists]).
    Replaces qiskit's group_commuting, which is O(N^2) and blows memory at ~10^4 Paulis (L=6 M2)."""
    codes = np.array([[_LET[c] for c in lab] for lab in labels], dtype=np.int8)  # (N, NQ)
    bases = np.zeros((0, codes.shape[1]), dtype=np.int8)
    members = []
    for i in range(len(labels)):
        p = codes[i]
        idx = -1
        if len(bases):
            compat = np.all((p == 0) | (bases == 0) | (bases == p), axis=1)
            if compat.any():
                idx = int(np.argmax(compat))
        if idx >= 0:
            b = bases[idx]; m = (b == 0); b[m] = p[m]
            members[idx].append(labels[i])
        else:
            bases = np.vstack([bases, p.copy()]); members.append([labels[i]])
    basis_strs = [''.join(_LETINV[c] for c in b) for b in bases]
    return basis_strs, members


def group_expectations(labels, counts, ntot):
    """vectorized <P> for all Paulis in one QWC group from its rotated-basis counts."""
    keys = list(counts.keys())
    cnt = np.array([counts[k] for k in keys], dtype=float)
    # bit matrix bm[i, j] = measured bit of qubit j in shot-key i  (qiskit key: leftmost = highest qubit)
    bm = np.array([[1 if ch == '1' else 0 for ch in k[::-1]] for k in keys], dtype=np.int8)
    out = {}
    for lab in labels:
        supp = [NQ - 1 - k for k, ch in enumerate(lab) if ch != 'I']
        if not supp:
            out[lab] = 1.0
        else:
            sign = 1 - 2 * (bm[:, supp].sum(axis=1) & 1)
            out[lab] = float((cnt * sign).sum() / ntot)
    return out


def group_basis(labels):
    """for a QWC group of Pauli labels, the per-qubit measurement letter (X/Y/Z or I)."""
    basis = ['I'] * NQ
    for lab in labels:
        for k, ch in enumerate(lab):
            if ch != 'I':
                qb = NQ - 1 - k
                if basis[NQ - 1 - qb] == 'I':
                    basis[k] = ch
    return basis


def rotation_circuit(basis_label):
    """circuit that rotates each qubit's measurement basis to Z (X->H, Y->Sdg;H)."""
    qc = QuantumCircuit(NQ)
    for k, ch in enumerate(basis_label):
        qb = NQ - 1 - k
        if ch == 'X':
            qc.h(qb)
        elif ch == 'Y':
            qc.sdg(qb); qc.h(qb)
    return qc


def estimate_moments(psi, Mops, ns=50000, seed=1, reps=1):
    """finite-shot counts estimate of <M> for each operator in Mops via ONE shared QWC grouping
    (all bases measured once per rep, reused across operators). psi: Statevector.
    Returns array (reps, len(Mops)) and the number of measurement bases."""
    groups, bases, coeffs = build_measurement_plan(Mops)
    rng = np.random.default_rng(seed)
    out = np.zeros((reps, len(Mops)))
    for r in range(reps):
        pexp = {}
        for labs, basis in zip(groups, bases):
            rotated = psi.evolve(rotation_circuit(basis))
            rotated.seed(int(rng.integers(1, 2**31)))
            counts = rotated.sample_counts(ns)
            pexp.update(group_expectations(labs, counts, ns))
        for j, c in enumerate(coeffs):
            out[r, j] = sum(c[l] * pexp[l] for l in c)
    return out, len(groups)


def chunked_dot(A, B, chunk=48):
    """(A.dot(B)).simplify() with peak memory bounded by chunk*|B| (avoids the L=6 M2 blowup)."""
    labels = A.to_list()
    acc = None
    for i in range(0, len(labels), chunk):
        piece = SparsePauliOp.from_list(labels[i:i + chunk]).dot(B).simplify()
        acc = piece if acc is None else (acc + piece).simplify()
    return acc


def build_operators(nmom=3):
    """return ([M0..M_{nmom-1}], E0, psi, m_exact). Operators via memory-bounded chunked products;
    exact moments computed matrix-free (phi=drho|psi>, m_k=<phi|(H-E0)^k|phi>) so they are cheap and
    exact even where the M_k SparsePauliOp is large (works at L=6, 12 qubits)."""
    H = build_hubbard()
    rho = build_rho()
    Hmat = H.to_matrix(sparse=True)
    E0v, V = sla.eigsh(Hmat, k=1, which='SA')
    E0 = float(E0v[0]); psivec = V[:, 0]; psi = Statevector(psivec)
    Iop = SparsePauliOp.from_list([('I' * NQ, 1.0)])
    rho_ev = float(psi.expectation_value(rho).real)
    drho = (rho - rho_ev * Iop).simplify()
    HmE = (H - E0 * Iop).simplify()
    # operators for the measurement plan (chunked to bound memory)
    Mops = [chunked_dot(drho, drho)]
    v = drho
    for _ in range(1, nmom):
        v = chunked_dot(HmE, v)
        Mops.append(chunked_dot(drho, v))
    # exact moments matrix-free: phi = drho psi ; m_k = <phi|(H-E0)^k|phi>  (cheap, exact, any L)
    import scipy.sparse as sp
    HmE_mat = (Hmat - E0 * sp.identity(Hmat.shape[0], format='csr')).tocsr()
    phi = drho.to_matrix(sparse=True) @ psivec
    m_exact, w = [], phi.copy()
    for _ in range(nmom):
        m_exact.append(float(np.vdot(phi, w).real))
        w = HmE_mat @ w
    return Mops, E0, psi, m_exact


def build_measurement_plan(Mops):
    """QWC grouping shared across all Mops -> (groups[labels], bases[str], coeffs[dict per Mop])."""
    coeffs = [{lab: complex(co).real for lab, co in M.to_list()} for M in Mops]
    all_labels = sorted(set().union(*[set(c) for c in coeffs]))
    bases, groups = qwc_group_labels(all_labels)
    return groups, bases, coeffs


def moments_from_group_counts(group_counts, groups, coeffs, ns):
    """combine retained per-group device counts into [m0,m1,...]. group_counts[i] matches groups[i]."""
    pexp = {}
    for labs, counts in zip(groups, group_counts):
        pexp.update(group_expectations(labs, counts, sum(counts.values())))
    return [sum(c[l] * pexp[l] for l in c) for c in coeffs]


def main():
    Mops, E0, psi, m_exact = build_operators(nmom=3)
    M0, M1, M2 = Mops

    print("=" * 70)
    print(f"DEVICE-STYLE MOMENT ESTIMATOR  (Hubbard L={L} chain, {NQ} qubits, U/t={U})")
    print("=" * 70)
    print(f"ground energy E0 = {E0:.4f}")
    print(f"operator sizes (Paulis): M0 {len(M0)}, M1 {len(M1)}, M2 {len(M2)}")
    print(f"EXACT moments (matrix-free):  m0={m_exact[0]:.5f}  m1={m_exact[1]:.5f}  m2={m_exact[2]:.5f}\n")

    est, ngroups = estimate_moments(psi, Mops, ns=50000, seed=1, reps=8)
    print(f"measurement bases (QWC groups / 'bond-basis' circuits): {ngroups}")
    print(f"COUNTS-BASED (Ns=50000, {est.shape[0]} reps):")
    oks = []
    for j, nm in enumerate(['m0', 'm1', 'm2']):
        mu, sd = est[:, j].mean(), est[:, j].std()
        se = sd / np.sqrt(est.shape[0])
        ok = abs(mu - m_exact[j]) < 4 * se + 1e-6
        oks.append(ok)
        print(f"  {nm}_hat = {mu:.5f} +/- {sd:.5f}   (exact {m_exact[j]:.5f}, "
              f"bias {mu-m_exact[j]:+.5f} = {abs(mu-m_exact[j])/max(se,1e-9):.1f} SE-of-mean)")
    print(f"\nVERIFIED: counts-based estimator recovers exact moments within shot noise? "
          f"m0:{oks[0]}  m1:{oks[1]}  m2:{oks[2]}")
    print("=> the reconstruction-independent m0 (counts-only) and m1,m2 (bond-basis, off-diagonal")
    print("   hopping via QWC rotated measurements) are DEVICE-MEASURABLE. On hardware only the")
    print("   state prep and the sampler backend change; this estimator is identical.")
    res = {'_provenance': {'script': 'bond_moment_estimator.py', 'sim_only': True,
                           'model': f'Hubbard L={L} chain OBC, {NQ} qubits, U/t={U}, density probe rho_q q=pi',
                           'claim': 'reconstruction-independent m0,m1,m2 (incl off-diagonal hopping) '
                                    'device-measurable via counts + bond-basis (QWC) circuits'},
           'E0': E0, 'ns': 50000, 'reps': int(est.shape[0]), 'n_measurement_bases': ngroups,
           'op_sizes_paulis': {'M0': len(M0), 'M1': len(M1), 'M2': len(M2)},
           'm_exact': {'m0': m_exact[0], 'm1': m_exact[1], 'm2': m_exact[2]},
           'm_hat': {f'm{j}': float(est[:, j].mean()) for j in range(3)},
           'm_std': {f'm{j}': float(est[:, j].std()) for j in range(3)}}
    import json
    RES = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '06_results'))
    _os.makedirs(RES, exist_ok=True)
    json.dump(res, open(_os.path.join(RES, f'2026-08-25_bond_moment_estimator_L{L}.json'), 'w'), indent=2)
    print(f"\nwrote 06_results/2026-08-25_bond_moment_estimator_L{L}.json")
    return res


if __name__ == '__main__':
    main()
