# -*- coding: utf-8 -*-
"""COMMITTED script for the estimator circuit-count scaling table (fixes the hand-written,
inflated table caught by self-audit S-1/S-3).  SIM-ONLY, structure-only (no exponential statevector).

The Pauli SET of M0,M1,M2 = drho (H-E0)^k drho depends on the OPERATOR structure only, not on the
numeric values of E0 or <rho> (those scale the identity term). E0 shifts identity only. For <rho>:
at half filling with q=pi and even L, <rho_q>=0 by particle-hole/staggering symmetry, so drho=rho
EXACTLY (no identity shift) -- this is why we may skip the exponential eigsh. (A NONZERO <rho> would
enlarge the support; we do NOT assume a generic nonzero placeholder -- that was the S-3 bug.)

Verified against the full eigsh-based build_operators at L=4,6,8,10 (where eigsh is affordable):
the label-set counts match exactly.
"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import bond_moment_estimator as bme
from qiskit.quantum_info import SparsePauliOp


def operators_structure(L):
    """build M0,M1,M2 SparsePauliOp with E0=0 and drho=rho (i.e. <rho>=0, exact by symmetry).
    No eigsh, no statevector."""
    old = bme.L, bme.NQ, bme.T, bme.U
    bme.L, bme.NQ = L, 2 * L
    H = bme.build_hubbard(); rho = bme.build_rho()
    # <rho>=0 by symmetry -> drho = rho (no identity subtraction)
    drho = rho.simplify()
    HmE = H.simplify()                        # E0=0 shifts identity only; irrelevant to the Pauli SET
    M0 = bme.chunked_dot(drho, drho)
    v1 = bme.chunked_dot(HmE, drho); M1 = bme.chunked_dot(drho, v1)
    v2 = bme.chunked_dot(HmE, v1);   M2 = bme.chunked_dot(drho, v2)
    bme.L, bme.NQ, bme.T, bme.U = old
    return M0, M1, M2


def main():
    rows = []
    for L in [4, 6, 8, 10, 12]:
        t = time.time()
        M0, M1, M2 = operators_structure(L)
        _, b_full, _ = bme.build_measurement_plan([M0, M1, M2])
        _, b_01, _ = bme.build_measurement_plan([M0, M1])
        rows.append({'L': L, 'qubits': 2 * L, 'M2_paulis': len(M2),
                     'full_battery_bases': len(b_full), 'm0m1_bases': len(b_01),
                     'build_s': round(time.time() - t, 1)})
        print(f"L={L:2d} ({2*L:2d}q): M2={len(M2):7d} Paulis | full battery {len(b_full):5d} bases | "
              f"(m0,m1) {len(b_01):3d} bases | {rows[-1]['build_s']}s", flush=True)

    # cross-check against the eigsh-based build_operators where affordable (L<=10)
    print("\ncross-check vs full eigsh build_operators (structure must match):")
    for L in [4, 6, 8]:
        os.environ['BME_L'] = str(L)
        import importlib; importlib.reload(bme)
        Mops, E0, psi, mex = bme.build_operators(nmom=3)
        _, bf, _ = bme.build_measurement_plan(Mops)
        ref = next(r for r in rows if r['L'] == L)
        match = (len(Mops[2]) == ref['M2_paulis']) and (len(bf) == ref['full_battery_bases'])
        print(f"  L={L}: eigsh M2={len(Mops[2])} bases={len(bf)}  vs structure M2={ref['M2_paulis']} "
              f"bases={ref['full_battery_bases']}  -> {'MATCH' if match else 'MISMATCH'}")

    out = {'_provenance': {'script': 'estimator_scaling.py', 'sim_only': True,
                           'method': 'structure-only Pauli label-set counts; <rho>=0 exact by half-filling '
                                     'symmetry (q=pi, even L), so drho=rho; verified vs eigsh build_operators at L<=8'},
           'columns': ['L', 'qubits', 'M2_paulis', 'full_battery_bases', 'm0m1_bases'],
           'scaling': [[r['L'], r['qubits'], r['M2_paulis'], r['full_battery_bases'], r['m0m1_bases']] for r in rows]}
    RES = os.path.normpath(os.path.join(HERE, '..', '06_results'))
    json.dump(out, open(os.path.join(RES, '2026-08-25_estimator_circuit_scaling.json'), 'w'), indent=2)
    print(f"\nwrote 06_results/2026-08-25_estimator_circuit_scaling.json (CORRECTED)")


if __name__ == '__main__':
    main()
