r"""Verify: <-That> needs 2 QWC settings, and those settings are a SUBSET of the
bases the M0/M1 measurement plan already requires (paper Sec. IV claim)."""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ["BME_L"] = "6"
import bond_moment_estimator as bme
from qiskit.quantum_info import SparsePauliOp

L, NQ = bme.L, bme.NQ
T = bme.T

# -That = +T * sum_{<ij>,s} (c^dag_i c_j + h.c.)  ->  JW (adjacent, blocked): +T/2 (XX + YY)
terms = []
for spin in (0, 1):
    for i in range(L - 1):
        a, b = spin * L + i, spin * L + i + 1
        terms.append((bme._pstr({a: 'X', b: 'X'}), +T * 0.5))
        terms.append((bme._pstr({a: 'Y', b: 'Y'}), +T * 0.5))
negT = SparsePauliOp.from_list(terms).simplify()
labels_negT = sorted({lab for lab, _ in negT.to_list()})
bases_T, groups_T = bme.qwc_group_labels(labels_negT)
print(f"<-T> Pauli strings: {len(labels_negT)}")
print(f"<-T> QWC settings : {len(bases_T)}  -> bases {bases_T}")

# measurement plan for M0, M1 (the deployed primary battery)
Mops, E0, psi, m_exact = bme.build_operators(nmom=2)
groups, bases, coeffs = bme.build_measurement_plan(Mops)
print(f"M0,M1 QWC settings: {len(bases)}")

def covers(basis_big, basis_small):
    """a setting `basis_big` can measure any Pauli whose per-qubit letters match where non-I."""
    return all(bs == 'I' or bs == bb for bb, bs in zip(basis_big, basis_small))

covered = []
for bt in bases_T:
    hit = [b for b in bases if covers(b, bt)]
    covered.append(bool(hit))
    print(f"  setting {bt}  covered by existing plan? {bool(hit)}  ({len(hit)} matching)")

# stronger check: is every -T Pauli string measurable in some existing plan basis?
def measurable(label, basis):
    return all(ch == 'I' or ch == basis[k] for k, ch in enumerate(label))
allm = all(any(measurable(l, b) for b in bases) for l in labels_negT)
print(f"\nEvery <-T> Pauli string measurable inside an existing M0/M1 setting: {allm}")
print(f"Numeric check  <-That> = {float(psi.expectation_value(negT).real):.10f}")
