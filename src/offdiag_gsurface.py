# -*- coding: utf-8 -*-
"""G-surface: exact local-depol bias floor of the density-probe m1 vs circuit depth G (CZ count),
at the ibm_fez 2q anchor eps=2.5e-3, with the global-depolarizing folklore overlay. Establishes the
depth dependence needed to state the paper's off-diagonal floor HONESTLY as a bracket
[exact static-carrier lower bound, folklore upper bound] at the ACTUAL companion depth (~293 CZ),
not the paper's understated '~10^2 CZ'. SIM-ONLY, zero QPU."""
import numpy as np, json, os, time
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error
import bond_moment_estimator as bme

EPS = 2.5e-3            # ibm_fez 2q anchor
SEED = 20260830
G_LIST = [40, 60, 100, 150, 200, 292]
Mops, E0, psi, m_exact = bme.build_operators(nmom=2)
M1d = Mops[1].to_matrix(sparse=True).toarray()
gs = np.asarray(psi.data, dtype=complex); NQ = bme.NQ
m1e = m_exact[1]
trM1 = next(float(np.real(co)) for lab, co in Mops[1].to_list() if set(lab) == {'I'})
Bmax = trM1 - m1e
SIM = AerSimulator(method='density_matrix'); SIM.set_options(fusion_enable=False)
PAIRS = [(i, i + 1) for i in range(NQ - 1)]

def floor_at(ncz):
    qc = QuantumCircuit(NQ); qc.set_statevector(gs)
    for k in range(ncz // 2):
        a, b = PAIRS[k % len(PAIRS)]; qc.cz(a, b); qc.cz(a, b)   # identity carrier
    qc.save_density_matrix()
    nm = NoiseModel(); nm.add_all_qubit_quantum_error(depolarizing_error(EPS, 2), 'cz')
    rho = np.asarray(SIM.run(qc, noise_model=nm, seed_simulator=SEED).result().data()['density_matrix'])
    return float(np.real(np.tensordot(rho, M1d.T, axes=([0, 1], [0, 1]))))

print(f"m1_exact={m1e:.4f}  Tr(M1)/dim={trM1:.4f}  Bmax={Bmax:.4f}  (eps={EPS:.1e})")
print(f"{'G':>5} {'B_exact':>9} {'%m1':>7} {'B_folklore':>11} {'%m1':>7} {'folk/exact':>10}")
rows = []
t0 = time.time()
for G in G_LIST:
    Be = floor_at(G) - m1e
    peff = 1.0 - (1.0 - EPS) ** G
    Bf = peff * Bmax
    rows.append({'G': G, 'B_exact': Be, 'pct_exact': 100 * Be / m1e,
                 'B_folklore': Bf, 'pct_folklore': 100 * Bf / m1e, 'ratio_folk_exact': Bf / Be})
    print(f"{G:>5} {Be:>9.4f} {100*Be/m1e:>6.1f}% {Bf:>11.4f} {100*Bf/m1e:>6.1f}% {Bf/Be:>10.2f}  [{time.time()-t0:.0f}s]")

RES = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data'))
json.dump({'eps': EPS, 'seed': SEED, 'm1_exact': m1e, 'TrM1_dim': trM1, 'Bmax': Bmax,
           'note': 'exact static-carrier local-depol floor is a LOWER bound (no operator scrambling, '
                   'no leakage, no readout); folklore global-depol is the upper edge; real Krylov-Trotter '
                   'circuit scrambles M1 to high weight -> true floor between exact and folklore',
           'gsurface': rows},
          open(os.path.join(RES, f'{time.strftime("%Y-%m-%d")}_offdiag_gsurface.json'), 'w'), indent=2)
print(f"\nwrote {time.strftime('%Y-%m-%d')}_offdiag_gsurface.json  runtime={time.time()-t0:.0f}s")
