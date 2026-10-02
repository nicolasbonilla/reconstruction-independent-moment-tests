# -*- coding: utf-8 -*-
"""G-surface: first-moment bias floor of the density-probe m1 versus circuit depth G (number of CZ
gates) at the ibm_fez two-qubit anchor eps = 2.5e-3, for three noise models:
  - exact local depolarizing (density matrix): each CZ followed by a two-qubit depolarizing channel;
  - global depolarizing ("folklore"): m1 -> (1-p) m1 + p Tr(M1)/dim, p = 1-(1-eps)^G;
  - coherent ZZ over-rotation theta = sqrt(eps) per CZ (statevector, no decoherence).
Surrogate probe: q = pi density probe of the open L = 6 chain at U/t = 4 (bond_moment_estimator.py).

Reference depths (read from data/sampler_cz_count.json, generator src/sampler_cz_count.py):
  - G_ref  = median transpiled CZ count of the 24 deployed ISA circuits (k >= 1) of our four L = 8
             Sampler jobs on ibm_fez (headline.G_L8_deployed_median);
  - G_L6   = the same median for the 12 deployed ISA circuits (k >= 1) of the two 12-qubit ibm_fez
             Sampler jobs of 2026-08-21, identified as build_circuit(k) of hardware_matched_job.py (the
             companion's L = 6 circuit family) at U/t = 4 and 3.5 (headline.G_L6_deployed_median);
             FakeFez transpilation of that family gives the same 213 for every circuit and seed.
Neither is the depth of the ground-state preparation the bond-basis estimator would need; they are
reference depths of the sampler circuits. 292 is kept as a grid point (one of the deployed counts).

Carrier model: the CZ unitaries are applied in cancelling pairs (cz;cz = I) on nearest-neighbour
pairs, so the state is unchanged and only the noise acts. Because the depolarizing channel is
unitarily covariant and RZZ commutes with CZ, this equals G applications of the error channel; an
odd G is completed by one extra error channel (depolarizing or RZZ(theta)) on the next pair, with no
unpaired CZ unitary. Even G is identical to the cz;cz carrier.

The local value is exact for this static noise model, which omits operator spreading by a real
circuit, leakage and readout error. SIM-ONLY, zero QPU.
Environment of the deposited output: qiskit 1.4.6, qiskit-aer 0.17.2, numpy 2.4.6, scipy 1.17.1.
Output: data/2026-08-31_offdiag_gsurface.json (fixed name; run_date inside).
"""
import numpy as np, json, os, time
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error
import bond_moment_estimator as bme

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, '..', 'data'))
OUT = os.path.join(RES, '2026-08-31_offdiag_gsurface.json')
CZJ = os.path.join(RES, 'sampler_cz_count.json')

EPS = 2.5e-3            # ibm_fez 2q anchor
SEED = 20260830
_cz = json.load(open(CZJ))['headline']
G_REF = int(round(_cz['G_L8_deployed_median']))
G_L6 = int(round(_cz['G_L6_deployed_median']))
G_LIST = sorted({40, 60, 100, 150, 200, 292, G_L6, G_REF})

Mops, E0, psi, m_exact = bme.build_operators(nmom=2)
M1d = Mops[1].to_matrix(sparse=True).toarray()
gs = np.asarray(psi.data, dtype=complex); NQ = bme.NQ
m1e = m_exact[1]
trM1 = next(float(np.real(co)) for lab, co in Mops[1].to_list() if set(lab) == {'I'})
Bmax = trM1 - m1e
SIM = AerSimulator(method='density_matrix'); SIM.set_options(fusion_enable=False)
SIM_SV = AerSimulator(method='statevector'); SIM_SV.set_options(fusion_enable=False)
PAIRS = [(i, i + 1) for i in range(NQ - 1)]


def floor_at(ncz):
    qc = QuantumCircuit(NQ); qc.set_statevector(gs)
    for k in range(ncz // 2):
        a, b = PAIRS[k % len(PAIRS)]; qc.cz(a, b); qc.cz(a, b)   # identity carrier
    if ncz % 2:
        a, b = PAIRS[(ncz // 2) % len(PAIRS)]
        qc.append(depolarizing_error(EPS, 2).to_instruction(), [a, b])  # error channel of the odd CZ
    qc.save_density_matrix()
    nm = NoiseModel(); nm.add_all_qubit_quantum_error(depolarizing_error(EPS, 2), 'cz')
    rho = np.asarray(SIM.run(qc, noise_model=nm, seed_simulator=SEED).result().data()['density_matrix'])
    return float(np.real(np.tensordot(rho, M1d.T, axes=([0, 1], [0, 1]))))


def coherent_at(ncz):
    th = np.sqrt(EPS)
    qc = QuantumCircuit(NQ); qc.set_statevector(gs)
    for k in range(ncz // 2):
        a, b = PAIRS[k % len(PAIRS)]
        qc.cz(a, b); qc.rzz(th, a, b); qc.cz(a, b); qc.rzz(th, a, b)
    if ncz % 2:
        a, b = PAIRS[(ncz // 2) % len(PAIRS)]
        qc.rzz(th, a, b)
    qc.save_statevector()
    sv = np.asarray(SIM_SV.run(qc).result().data()['statevector'])
    return float(np.real(sv.conj() @ (M1d @ sv)))


print(f"m1_exact={m1e:.4f}  Tr(M1)/dim={trM1:.4f}  Bmax={Bmax:.4f}  (eps={EPS:.1e})  G_ref={G_REF}  G_L6={G_L6}")
print(f"{'G':>5} {'B_exact':>9} {'%m1':>7} {'B_folklore':>11} {'%m1':>7} {'B_coher':>9} {'%m1':>7} {'folk/exact':>10}")
rows = []
t0 = time.time()
for G in G_LIST:
    Be = floor_at(G) - m1e
    Bc = coherent_at(G) - m1e
    peff = 1.0 - (1.0 - EPS) ** G
    Bf = peff * Bmax
    rows.append({'G': G, 'B_exact': Be, 'pct_exact': 100 * Be / m1e,
                 'B_folklore': Bf, 'pct_folklore': 100 * Bf / m1e, 'ratio_folk_exact': Bf / Be,
                 'B_coherent': Bc, 'pct_coherent': 100 * Bc / m1e})
    print(f"{G:>5} {Be:>9.4f} {100*Be/m1e:>6.1f}% {Bf:>11.4f} {100*Bf/m1e:>6.1f}% {Bc:>9.4f} {100*Bc/m1e:>6.1f}% "
          f"{Bf/Be:>10.2f}  [{time.time()-t0:.0f}s]")

json.dump({'eps': EPS, 'seed': SEED, 'm1_exact': m1e, 'TrM1_dim': trM1, 'Bmax': Bmax,
           'generator': 'src/offdiag_gsurface.py', 'run_date': time.strftime('%Y-%m-%d'),
           'G_reference': G_REF,
           'G_reference_source': 'data/sampler_cz_count.json headline.G_L8_deployed_median '
                                 '(median CZ of the 24 deployed L=8 ISA circuits with k>=1)',
           'G_L6': G_L6,
           'G_L6_source': 'data/sampler_cz_count.json headline.G_L6_deployed_median',
           'note': 'exact static-carrier local-depolarizing floor: no operator spreading, no leakage, no '
                   'readout error; global depolarizing is the folklore overlay; coherent = ZZ over-rotation '
                   'theta=sqrt(eps) per CZ, pure state. Odd G: one extra error channel without an unpaired CZ.',
           'gsurface': rows},
          open(OUT, 'w'), indent=2)
print(f"\nwrote {OUT}  runtime={time.time()-t0:.0f}s")
