# -*- coding: utf-8 -*-
"""PREPARED, NEVER RUN: a minimal blinded hardware job for the bond-basis moment estimator (ibm_fez / Heron).

STATUS: this job has never been submitted to hardware, and its state preparation is a PLACEHOLDER (below),
not a ground-state preparation. No moment in this work was estimated from device counts. What the script
does today is a local Aer dry run (default) of the same circuits and estimator.

If it were run with a real ground-state preparation, it would estimate the reconstruction-independent
moments m0 (computational-basis counts only) and m1 (off-diagonal hopping, via qubit-wise-commuting
'bond-basis' rotations) of bond_moment_estimator.py from retained per-basis counts, with shot-noise error
bars, so that Delta_k could be formed against a reconstruction. The estimator itself has been checked only
in noiseless simulation (bond_moment_estimator.py); on hardware the state prep and the sampler backend
change, and at ibm_fez depth the off-diagonal moments are forecast to carry a depolarizing bias of about
37-116% of the signal (offdiag_noise_forecast.py, offdiag_gsurface.py).

============================  READ BEFORE RUNNING  ============================
0. ENVIRONMENT (blocking): this account is on the NEW IBM Quantum Platform (quantum.cloud.ibm.com,
   IBM Cloud), NOT the retired quantum-computing.ibm.com. It needs qiskit>=1.2 with a matching
   qiskit-ibm-runtime>=0.34 (the pair installed as of writing, qiskit 1.0.2 + runtime 0.24.0, is
   INCOMPATIBLE: runtime import fails on SamplerPubResult). Use an ISOLATED venv so other projects
   are untouched:
       py -m venv .venv-ibmqpu && .venv-ibmqpu\Scripts\activate
       pip install -U "qiskit>=1.2" "qiskit-ibm-runtime>=0.34" qiskit-aer numpy scipy
1. CREDENTIALS (USER-ONLY): from IBM Quantum Platform save the account ONCE in a private shell
   (NOT this repo, NOT chat). New-platform form (API key + instance CRN, both on your dashboard):
       from qiskit_ibm_runtime import QiskitRuntimeService
       QiskitRuntimeService.save_account(
           channel="ibm_quantum_platform",   # use "ibm_cloud" if your runtime lacks this channel
           token="<YOUR_API_KEY>",           # Dashboard -> API key -> Create
           instance="<YOUR_CRN>",            # Dashboard -> Instances -> open-instance -> CRN
           name="open-instance", overwrite=True)
   This script loads the saved account by name (IBM_ACCOUNT env, default 'open-instance') and never
   sees the raw key. Budget note: NMOM=2 (m0,m1) ~35 circuits ~3-4 min QPU fits the Open 10 min/cycle.
2. No QPU contact unless RUN=1 AND the saved account exists. Unset RUN -> local Aer dry-run that
   exercises the exact same circuits + estimator so you can inspect the numbers first.
3. Probe = density rho_q (diagonal) so m0 is counts-only; m1's off-diagonal hopping is read from
   the QWC bond-basis circuits this script builds. Truth is classically known at L=6 (a run would give
   a shot-noisy Delta_k, not certification of a classically-hard line shape). Call the
   estimator "reconstruction-independent", NEVER "sample-free".
4. The PREP below is a DOCUMENTED PLACEHOLDER (one Trotter-like layer from a product state), not a
   ground-state preparation; its dry-run moments are those of the prepared state, not of the Hubbard
   ground state. A meaningful run needs a ground-state preparation of the L-site chain. The companion's
   sampling circuits are NOT a substitute: they prepare (N+1)-particle states for the spectrum (a
   product determinant followed by on-site R_zz and R_xx R_yy hopping layers without Jordan-Wigner
   strings), not the ground state whose moments are wanted here.
5. Options matched to the companion: dynamical decoupling (XpXm), Pauli gate twirling
   (num_randomizations=32), measurement twirling. No readout-error mitigation: measurement twirling
   returns twirled raw counts (it mitigates nothing by itself), and the resilience line in qpu_run fails on SamplerV2,
   which has no resilience options. Never run on hardware. Submit as ONE Batch. RETAIN raw per-basis
   counts -> data/heron_counts_<jobid>.json.
==============================================================================
"""
import os, sys, json
import numpy as np
from qiskit import QuantumCircuit, transpile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bond_moment_estimator as bme     # L defaults to 6 (12 qubits) via BME_L; see that module

NS = 50000                              # do NOT change: interval thresholds are fixed to this budget
NMOM = int(os.environ.get("NMOM", 3))   # 3 = full battery m0,m1,m2 (~357 bond-basis circuits, ~12-24 min QPU);
                                        # 2 = primary m0,m1 falsifier (~35 circuits, ~1.5-3 min QPU -> fits Open Plan).
                                        # NS stays 50000 either way, so the pre-registered thresholds remain valid.
BACKEND = os.environ.get("IBM_BACKEND", "ibm_fez")        # 156q Heron r2, confirmed available on your instance
ACCOUNT = os.environ.get("IBM_ACCOUNT", "open-instance")  # saved-account name (your dashboard: 'open-instance')
NQ = bme.NQ
OUT = os.path.normpath(os.path.join(HERE, '..', 'data'))


def build_prep_circuit():
    """PLACEHOLDER prep: one first-order Trotter-like layer of the L Hubbard chain on 2L qubits, from a
    staggered occupation. Not a ground-state preparation; replace it with one before any real run."""
    L = bme.L
    qc = QuantumCircuit(NQ, name=f"prep_L{L}")
    for i in range(L):                          # staggered initial occupation (half filling-ish)
        if i % 2 == 0:
            qc.x(i)                              # up on even sites
        else:
            qc.x(L + i)                          # dn on odd sites
    dt = 0.35
    for spin in (0, 1):
        for i in range(L - 1):
            a, b = spin * L + i, spin * L + i + 1
            qc.rxx(dt, a, b); qc.ryy(dt, a, b)
    for i in range(L):
        qc.cp(dt * bme.U, i, L + i)
    return qc


def measurement_circuits(prep):
    """for each QWC group build prep + basis-rotation + measure; return (circuits, groups, coeffs, Mops)."""
    Mops, E0, _psi, _mex = bme.build_operators(nmom=NMOM)
    groups, bases, coeffs = bme.build_measurement_plan(Mops)
    circuits = []
    for basis in bases:
        qc = prep.copy()
        qc.compose(bme.rotation_circuit(basis), inplace=True)
        qc.measure_all()
        circuits.append(qc)
    return circuits, groups, coeffs, Mops, E0


def local_dry_run():
    """No QPU: Aer-sample the exact same circuits + estimator; compare counts-based m0,m1 to exact."""
    from qiskit.quantum_info import Statevector
    prep = build_prep_circuit()
    circuits, groups, coeffs, Mops, E0 = measurement_circuits(prep)
    # exact reference on the PREPARED state (what the device actually prepares here)
    psi_prep = Statevector(prep)
    m_exact = [float(psi_prep.expectation_value(M).real) for M in Mops]
    try:
        from qiskit_aer import AerSimulator
        sim = AerSimulator()
    except Exception:
        print("Aer unavailable; transpile check only."); return
    group_counts = []
    for qc in circuits:
        tqc = transpile(qc, sim)
        res = sim.run(tqc, shots=NS, seed_simulator=1).result()
        group_counts.append(res.get_counts())
    m_hat = bme.moments_from_group_counts(group_counts, groups, coeffs, NS)
    print(f"[dry-run] L={bme.L}, {NQ} qubits, {len(circuits)} bond-basis circuits, Ns={NS}")
    for j in range(NMOM):
        print(f"  m{j}: counts-based {m_hat[j]:.4f}  vs exact(prepared state) {m_exact[j]:.4f}  "
              f"Delta={abs(m_hat[j]-m_exact[j]):.4f}")
    json.dump({'_dry_run': True, 'backend': 'aer', 'shots': NS, 'n_bases': len(circuits),
               'm_hat': m_hat, 'm_exact': m_exact, 'group_counts': group_counts},
              open(os.path.join(OUT, 'heron_counts_DRYRUN.json'), 'w'))
    print(f"[dry-run] RETAINED per-basis counts -> heron_counts_DRYRUN.json")
    print("Before any real run: replace build_prep_circuit (a placeholder, not a ground-state prep); then set RUN=1.")


def qpu_run():
    from qiskit_ibm_runtime import QiskitRuntimeService, Batch, SamplerV2 as Sampler
    service = QiskitRuntimeService(name=ACCOUNT)
    backend = service.backend(BACKEND)
    prep = build_prep_circuit()
    circuits, groups, coeffs, Mops, E0 = measurement_circuits(prep)
    tcircs = transpile(circuits, backend, optimization_level=3)
    with Batch(backend=backend) as batch:
        sampler = Sampler(mode=batch)
        o = sampler.options
        o.default_shots = NS
        o.dynamical_decoupling.enable = True
        o.dynamical_decoupling.sequence_type = "XpXm"
        o.twirling.enable_gates = True
        o.twirling.enable_measure = True
        o.twirling.num_randomizations = 32
        o.resilience.measure_mitigation = True   # kept as written; fails on SamplerV2 (no resilience options)
        job = sampler.run(tcircs)
        jid = job.job_id()
        print(f"submitted job {jid} on {BACKEND} ({len(tcircs)} bond-basis circuits, Batch); RETAINING counts")
        res = job.result()
    group_counts = [res[i].data.meas.get_counts() for i in range(len(tcircs))]
    m_hat = bme.moments_from_group_counts(group_counts, groups, coeffs, NS)
    json.dump({'_dry_run': False, 'backend': BACKEND, 'job_id': jid, 'shots': NS,
               'n_bases': len(tcircs), 'm_hat': m_hat, 'group_counts': group_counts},
              open(os.path.join(OUT, f'heron_counts_{jid}.json'), 'w'))
    print(f"counts-based moments m_hat = {[round(x,4) for x in m_hat]};  RETAINED counts -> heron_counts_{jid}.json")
    print("Next: compare m_hat with the reconstruction's back-moments to form a shot-noisy Delta_k (off-diagonal")
    print("      moments carry a forecast depolarizing bias of ~37-116% at this depth); then unblind.")


def main():
    os.makedirs(OUT, exist_ok=True)
    if os.environ.get("RUN") == "1":
        print("RUN=1 -> QPU. Use your own saved account; revoke any API key that was ever exposed before running.")
        qpu_run()
    else:
        print("RUN unset -> LOCAL Aer DRY-RUN (no QPU). Set RUN=1 to submit to ibm_fez.")
        local_dry_run()


if __name__ == '__main__':
    main()
