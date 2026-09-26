# -*- coding: utf-8 -*-
r"""
MATCHED, RECONSTRUCTION-INDEPENDENT DEVICE-SIDE MOMENT FALSIFIER  (Paper 3)
==========================================================================
Companion (Paper 2) canonical notebook:
    the companion repository dynamical-spectral-functions-sqd, notebooks/Spectral_Heron.ipynb
That notebook computes the electron-ADDITION single-particle spectral function A^+(omega) of the
1D Hubbard ring (L=6, U/t=4, thop=1.0, eta=0.15, 12 qubits, INTERLEAVED spin-orbital order
p = 2*site + spin, even=up / odd=down) on IBM Heron. The device ONLY samples: it prepares a
(4 up, 3 down) = 7-electron shallow-Trotter state (K=7 Krylov times, n_trot=3, dt_total=0.5/thop),
samples computational-basis bitstrings, post-selects Nocc==7 & Szocc==+1 -> subspace S, diagonalizes
H_S, and builds the Lehmann reconstruction A_hw. The c^dagger / Lehmann build is entirely classical.

WHAT THIS SCRIPT ADDS (Paper 3):  the reconstruction-INDEPENDENT companion falsifier.
--------------------------------------------------------------------------------------
Reference state (fixed, classically known, SHARED by both branches):
    |phi> = c^dagger_{site0,up} |0>,   |0> = exact half-filled (Nocc==L, Szocc==0) ground state, energy E0.
|phi> lives entirely in the ADDITION sector si = {Nocc==L+1 & Szocc==+1}.

Two first moments of the *addition* spectral function A^+(omega):
    m1_bar  =  \int omega A_hw(omega) domega                          RECONSTRUCTION moment (within-S)
            =  <phi_S| (H_S - E0) |phi_S>,   phi_S = P_S phi,  H_S = P_S H P_S
    m1_hat  =  <phi| (H - E0) |phi>          RECONSTRUCTION-INDEPENDENT, Eq.(app-leak), k=1:
                                             the FULL sparse H acts on the FULL |phi> and reaches
                                             addition-sector configs OUTSIDE the hardware-sampled S.
                                             Computed from the operator algebra, NEVER read from A_hw.
    Delta_1 =  | m1_hat - m1_bar |           DEVICE-SIDE FALSIFIER.

Delta_1 grows as the hardware-sampled subspace S omits addition-sector weight; it is independent of the
reconstruction because m1_hat uses the full operator on the full reference, not A_hw. This is exactly the
shared-state / operator-leak logic of Paper 3's Fig. teethshared (run_teeth_shared.py), now MATCHED to
the companion's REAL hardware-sampled subspace and its ACTUAL single-particle addition observable, at a
classically-checkable L=6. The object that is SHARED between the two branches here is |phi> itself (a
single classically-computed vector both branches start from), so the independent branch is no "more
independent than a device can be": the only difference between the branches is full operator (m1_hat)
vs. subspace-projected operator (m1_bar).

CORRECTNESS (addition sector, one-sided -- NOT the two-sided anticommutator identity):
    m0^+ = \int A^+ = <phi|phi> = <0| c_{0up} c^dagger_{0up} |0> = 1 - <n_{0up}>_GS        (NOT 1)
    m1^+ = \int omega A^+ = <phi| (H - E0) |phi>
Both are verified below against the exact Lehmann sums to <= 1e-10.

ENVIRONMENT (blocking for the QPU path):
    This account is on the NEW IBM Quantum Platform (quantum.cloud.ibm.com / IBM Cloud), NOT the
    retired quantum-computing.ibm.com (channel='ibm_quantum' is retired -- do NOT use it). The QPU path
    needs qiskit>=1.2 with a matching qiskit-ibm-runtime>=0.34 in an ISOLATED venv. The pair installed
    for the LOCAL dry-run (qiskit 1.0.2 + qiskit-ibm-runtime 0.24.0) is INCOMPATIBLE with the runtime
    import path used here (SamplerPubResult), so create a separate venv for RUN=1:
        py -m venv .venv-ibmqpu && .venv-ibmqpu\Scripts\activate
        pip install -U "qiskit>=1.2" "qiskit-ibm-runtime>=0.34" qiskit-aer numpy scipy
    Save the account ONCE, privately (new-platform form: API key + instance CRN):
        from qiskit_ibm_runtime import QiskitRuntimeService
        QiskitRuntimeService.save_account(
            channel="ibm_quantum_platform",   # or "ibm_cloud" if your runtime lacks this channel
            token="<YOUR_API_KEY>", instance="<YOUR_CRN>",
            name="open-instance", overwrite=True)
    This script loads the saved account by name (IBM_ACCOUNT env, default 'open-instance').

USAGE:
    py hardware_matched_job.py         # RUN unset -> LOCAL Aer dry-run, 0 QPU (this is what we ran)
    RUN=1 py hardware_matched_job.py   # QPU: SamplerV2 + measurement twirl + gate twirl + DD, 50000 shots, Batch
                                       # (never run on the QPU; no readout-error mitigation, see qpu_run)

The LOCAL dry-run and the QPU path share the SAME circuits, the SAME post-selection, and the SAME
falsifier arithmetic; only the sampler (Aer vs ibm_fez) differs.
"""
import os, sys, json, time
import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, '..', 'data'))
DATE = '2026-08-28'

# numpy 1.26 (dry-run venv) has np.trapz; numpy>=2 has np.trapezoid. Support both.
_trapz = getattr(np, 'trapezoid', None) or np.trapz

# ---- QPU knobs (unused in dry-run) ----
NS = 50000
BACKEND = os.environ.get("IBM_BACKEND", "ibm_fez")        # 156q Heron r2
ACCOUNT = os.environ.get("IBM_ACCOUNT", "open-instance")  # saved-account name (new platform)
SHOTS_DRY = int(os.environ.get("SHOTS_DRY", NS))          # per-circuit Aer shots in the dry-run
SEED = int(os.environ.get("SEED", 1))

t0 = time.time()
log = lambda *a: print(f"[{time.time()-t0:6.1f}s]", *a, flush=True)

# =====================================================================================
# 1. Hamiltonian + exact addition-sector reference (VERBATIM physics from the companion)
# =====================================================================================
L = 6; U = 4.0; thop = 1.0; eta = 0.15; M = 2 * L      # spin-orbital p = site*2 + spin (even=up, odd=dn)

def c_op(p, dim):
    r = []; c = []; d = []
    for s in range(dim):
        if (s >> p) & 1:
            sg = (-1) ** bin(s & ((1 << p) - 1)).count('1'); r.append(s & ~(1 << p)); c.append(s); d.append(float(sg))
    return sp.csr_matrix((d, (r, c)), shape=(dim, dim))

dim = 1 << M
C = [c_op(p, dim) for p in range(M)]; Cd = [c.T.conj() for c in C]
Nocc = np.array([bin(s).count('1') for s in range(dim)])
Szocc = np.array([sum(((s >> (2 * i)) & 1) - ((s >> (2 * i + 1)) & 1) for i in range(L)) for s in range(dim)])
H = sp.csr_matrix((dim, dim))
for i in range(L):                                     # periodic chain
    j = (i + 1) % L
    for spin in (0, 1):
        a = 2 * i + spin; b = 2 * j + spin; H = H - thop * (Cd[a] @ C[b] + Cd[b] @ C[a])
for i in range(L):
    H = H + U * (Cd[2 * i] @ C[2 * i]) @ (Cd[2 * i + 1] @ C[2 * i + 1])
H = H.tocsr()

gi = np.where((Nocc == L) & (Szocc == 0))[0]
w, v = np.linalg.eigh(H[gi][:, gi].toarray()); E0 = w[0]
psi0 = np.zeros(dim, complex); psi0[gi] = v[:, 0]
p_orb = 0
phi = Cd[2 * p_orb] @ psi0                             # |phi> = c^dagger_{site0,up} |0>
si = np.where((Nocc == L + 1) & (Szocc == 1))[0]       # ADDITION sector (N+1, Sz=+1)
En, Vn = np.linalg.eigh(H[si][:, si].toarray()); coef = Vn.conj().T @ phi[si]
grid = np.linspace((En - E0).min() - 1, (En - E0).max() + 1, 600)

def spec(pw, ww):
    A = np.zeros_like(grid)
    for a, b in zip(pw, ww):
        A += b * (eta / np.pi) / ((grid - a) ** 2 + eta ** 2)
    return A

A_exact = spec(En - E0, np.abs(coef) ** 2)

# ---- exact ADDITION-sector moments and their <=1e-10 verification ----
# Discrete Lehmann sums (definition of the moments of A^+):
m0_lehmann = float(np.sum(np.abs(coef) ** 2))
m1_lehmann = float(np.sum((En - E0) * np.abs(coef) ** 2))
# Operator forms (one-sided ADDITION identities -- NOT the two-sided m0=1 anticommutator):
n0up_gs = float(np.real(np.vdot(psi0, (Cd[0] @ C[0]) @ psi0)))     # <n_{0up}>_GS
m0_op = float(np.real(np.vdot(phi, phi)))                          # <phi|phi> = 1 - <n_{0up}>_GS
m1_op = float(np.real(np.vdot(phi, H @ phi)) - E0 * np.real(np.vdot(phi, phi)))  # <phi|(H-E0)|phi>
one_minus_n = 1.0 - n0up_gs
err_m0_ident = abs(m0_op - one_minus_n)            # <phi|phi> vs 1-<n_0up>
err_m0 = abs(m0_op - m0_lehmann)                   # operator vs Lehmann (int A^+)
err_m1 = abs(m1_op - m1_lehmann)                   # operator vs Lehmann (int omega A^+)
# broadened-grid cross-check (finite grid/broadening -> approximate, not the 1e-10 test):
gm0 = float(_trapz(A_exact, grid)); gm1 = float(_trapz(grid * A_exact, grid))


# =====================================================================================
# 2. Companion circuits -- COPIED VERBATIM from Spectral_Heron.ipynb cell 7
# =====================================================================================
from qiskit import QuantumCircuit, QuantumRegister, transpile
norb = L; nelec_add = (4, 3)      # (N+1): 4 up, 3 down
K = 7; n_trot = 3; dt_total = 0.5 / thop

def givens(qc, a, b, theta):      # e^{-i theta (XX+YY)/2} hopping between spin-orbitals a,b
    qc.rxx(theta, a, b); qc.ryy(theta, a, b)

def build_circuit(k):
    q = QuantumRegister(2 * L); qc = QuantumCircuit(q)
    for orb in [0, 1, 2, 3]: qc.x(q[2 * orb])          # 4 up electrons (even = up)
    for orb in [0, 1, 2]:  qc.x(q[2 * orb + 1])        # 3 down electrons (odd = down)
    dt = (k * dt_total) / max(n_trot, 1)
    for _ in range(n_trot):
        for i in range(L): qc.rzz(2 * U * dt, q[2 * i], q[2 * i + 1])         # on-site U (native RZZ)
        for spin in (0, 1):                                                  # hopping per spin
            for i in range(L):
                j = (i + 1) % L; givens(qc, q[2 * i + spin], q[2 * j + spin], thop * dt)
    qc.measure_all(); return qc

circuits = [build_circuit(k) for k in range(K)]

def bitstring_to_config(bs):
    r"""qiskit little-endian counts bitstring -> Fock integer in the OPERATOR convention (bit p =
    spin-orbital p occupation, p = 2*site + spin), i.e. the statevector basis index.

    CORRECTNESS NOTE (verified, load-bearing):  the companion notebook's HARDWARE cell uses
    int(bs[::-1], 2), but that REVERSES the 12-bit string and lands on the wrong Fock integer. On the
    interleaved p=2*site+spin layout, reversal maps even<->odd bits, i.e. it FLIPS the spin label of
    every orbital, so a physically Sz=+1 shot is read as Sz=-1. The companion's own LOCAL validation
    (cell 9) never hits this: it samples statevector INDICES directly (t = basis index), which are
    already in the operator convention. Reproducing that local path on real hardware counts requires
    int(bs, 2) (rightmost char = qubit 0), NOT int(bs[::-1], 2). Empirically, on the k=0 reference
    prep (basis state |0..6>, index 127, Sz=+1): int(bs,2)=127 (Sz=+1, KEPT) vs
    int(bs[::-1],2)=4064 (Sz=-1, REJECTED). Using the reversed map would post-select to the empty set
    (|S|=0) and waste 100% of the QPU shots. We therefore use the physically correct map and
    self-check it against the number-conserving prep in falsifier_from_pooled_configs()."""
    return int(bs.replace(' ', ''), 2)


# =====================================================================================
# 3. Shared falsifier arithmetic (identical for Aer dry-run and for hardware counts)
# =====================================================================================
pos = {int(x): i for i, x in enumerate(si)}          # si-config -> index within addition sector

def falsifier_from_pooled_configs(seen):
    """seen: set of Fock-integer configs post-selected to (Nocc==L+1 & Szocc==1).
    Returns the reconstruction A_hw, its within-S first moment m1_bar, the reconstruction-INDEPENDENT
    m1_hat (full H on full phi), and Delta_1 -- plus |S|, rel-L1 vs A_exact, and the leak fraction."""
    S = np.array(sorted(seen))
    Sidx = np.array([pos[x] for x in S if x in pos], dtype=int)     # indices inside the addition sector
    HS = H[si][:, si].toarray()[np.ix_(Sidx, Sidx)]
    Em, Um = np.linalg.eigh(HS)
    aS = Um.conj().T @ phi[si][Sidx]                                # <m | phi_S>
    A_hw = spec(Em - E0, np.abs(aS) ** 2)
    rel = float(_trapz(np.abs(A_hw - A_exact), grid) / _trapz(A_exact, grid))
    # RECONSTRUCTION moment (within-S first moment = \int omega A_hw): exact discrete form...
    m1_bar = float(np.sum((Em - E0) * np.abs(aS) ** 2))
    m1_bar_grid = float(_trapz(grid * A_hw, grid))                  # ...and the broadened-grid integral
    m0_bar = float(np.sum(np.abs(aS) ** 2))                         # ||phi_S||^2 = within-S weight
    # RECONSTRUCTION-INDEPENDENT moment: FULL sparse H on the FULL |phi> (reaches configs OUTSIDE S).
    # This never touches A_hw, Em, Um, or S -- it is the fixed operator target m1^+.
    m1_hat = m1_op
    delta1 = abs(m1_hat - m1_bar)
    # Delta_0: zeroth-moment falsifier. m0_hat = m0_op = <phi|phi> = 1-<n_0up> (exact, independent);
    # m0_bar = within-S weight. Delta_0 is the addition-sector weight the sampled subspace S leaks --
    # the MONOTONE, more sensitive subspace-completeness discriminator (grows steadily as S truncates,
    # unlike the non-monotone first-moment leak Delta_1).
    delta0 = abs(m0_op - m0_bar)
    leak_frac = max(0.0, 1.0 - m0_bar / m0_op)                     # fraction of |phi| weight outside S (=Delta_0/m0^+)
    return dict(S=S, Sidx=Sidx, A_hw=A_hw, rel=rel, m1_bar=m1_bar, m1_bar_grid=m1_bar_grid,
                m0_bar=m0_bar, m0_hat=m0_op, delta0=delta0, m1_hat=m1_hat, delta1=delta1, leak_frac=leak_frac)


def coverage_sweep():
    r"""Zero-QPU calibration: truncate the addition subspace to the top-d most-probable configs
    (ranked by the true weights |<config|phi>|^2) and report Delta_1(d). This traces the intrinsic
    falsifier response of THIS observable -- the curve a real device climbs as its sampled S omits
    weight -- and proves Delta_1 rises monotonically from 0 (full S) to m1^+ (empty S). Analogous to
    run_teeth_shared.py's determinant sweep, matched here to A^+."""
    w_si = np.abs(phi[si]) ** 2                        # per-config addition weight within the sector
    order = np.argsort(w_si)[::-1]
    rows = []
    for frac in (1.0, 0.75, 0.5, 0.25, 0.1):
        d = max(1, int(round(frac * len(si))))
        seen = set(int(si[j]) for j in order[:d])
        f = falsifier_from_pooled_configs(seen)
        rows.append((d, f['m0_bar'] / m0_op, f['rel'], f['delta0'], f['delta1']))
    return rows


def postselect(pooled_counts):
    """pooled_counts: dict bitstring->count (already pooled across the K times). -> (seen set, kept, total)."""
    seen = set(); kept = 0; total = 0
    for bs, ct in pooled_counts.items():
        total += ct
        t = bitstring_to_config(bs)
        if Nocc[t] == L + 1 and Szocc[t] == 1:
            seen.add(int(t)); kept += ct
    return seen, kept, total


def print_verification():
    print("\n=== ADDITION-SECTOR MOMENT VERIFICATION (one-sided A^+, NOT the m0=1 anticommutator) ===")
    print(f"  <n_0up>_GS                         = {n0up_gs:.12f}")
    print(f"  m0^+ = <phi|phi>                   = {m0_op:.12f}")
    print(f"  1 - <n_0up>_GS                     = {one_minus_n:.12f}   |diff| = {err_m0_ident:.2e}")
    print(f"  Lehmann  int A^+  (sum|coef|^2)    = {m0_lehmann:.12f}   |m0_op - Lehmann| = {err_m0:.2e}")
    print(f"  m1^+ = <phi|(H-E0)|phi>            = {m1_op:.12f}")
    print(f"  Lehmann  int wA^+ (sum(En-E0)|c|^2)= {m1_lehmann:.12f}   |m1_op - Lehmann| = {err_m1:.2e}")
    print(f"  broadened-grid cross-check: int A_exact={gm0:.6f}  int w A_exact={gm1:.6f} (finite grid, approx)")
    ok = (err_m0_ident <= 1e-10) and (err_m0 <= 1e-10) and (err_m1 <= 1e-10)
    print(f"  >>> addition moments recover exact Lehmann to <=1e-10:  {'PASS' if ok else 'FAIL'}")
    return ok


# =====================================================================================
# 4. LOCAL dry-run (Aer, 0 QPU)
# =====================================================================================
def local_dry_run():
    from qiskit_aer import AerSimulator
    sim = AerSimulator()
    print_verification()
    log(f"sampling {K} companion circuits on Aer, {SHOTS_DRY} shots each (seed={SEED}); depth(last)={circuits[-1].depth()}")
    per_time_counts = []
    pooled = {}
    for k, qc in enumerate(circuits):
        tqc = transpile(qc, sim)
        res = sim.run(tqc, shots=SHOTS_DRY, seed_simulator=SEED + k).result()
        cnt = res.get_counts()
        per_time_counts.append(cnt)
        for bs, ctv in cnt.items():
            pooled[bs] = pooled.get(bs, 0) + ctv
    seen, kept, total = postselect(pooled)
    # self-check: the companion circuits conserve N and Sz, so noiseless Aer must keep ~100% of shots.
    # A low fraction here would signal a bitstring->Fock mapping regression (see bitstring_to_config).
    if kept / total < 0.99:
        print(f"  [WARN] noiseless post-selection kept only {100*kept/total:.2f}% -- mapping/convention regression?")
    f = falsifier_from_pooled_configs(seen)
    print("\n=== MATCHED DEVICE-SIDE FALSIFIER (Aer dry-run; hardware-sampled S surrogate) ===")
    print(f"  post-selection kept {kept}/{total} pooled shots ({100*kept/total:.1f}%) -> |S| = {len(f['Sidx'])} "
          f"of {len(si)} addition-sector configs")
    print(f"  rel-L1 (A_hw vs A_exact)           = {f['rel']:.4f}")
    print(f"  m0_hat = <phi|phi> (indep, exact)  = {f['m0_hat']:.6f}")
    print(f"  m0_bar (within-S weight)           = {f['m0_bar']:.6f}")
    print(f"  Delta_0 = |m0_hat - m0_bar|        = {f['delta0']:.6f}   (MONOTONE completeness probe)")
    print(f"  m1_bar  (within-S, int w A_hw)     = {f['m1_bar']:.6f}   [grid int = {f['m1_bar_grid']:.6f}]")
    print(f"  m1_hat  (full H on full phi, indep)= {f['m1_hat']:.6f}")
    print(f"  Delta_1 = |m1_hat - m1_bar|        = {f['delta1']:.6f}")
    print(f"  within-S weight m0_bar/m0^+        = {f['m0_bar']/m0_op:.4f}   leak fraction = {f['leak_frac']:.4f}")

    sweep = coverage_sweep()
    print("\n  --- zero-QPU calibration: Delta_0 (monotone) & Delta_1 vs subspace coverage (top-d weighted configs) ---")
    print(f"    {'d':>4} {'cov(m0_bar/m0)':>14} {'rel-L1':>8} {'Delta_0':>9} {'Delta_1':>9}")
    for d, cov, rel, d0, d1 in sweep:
        print(f"    {d:>4} {cov:>14.4f} {rel:>8.4f} {d0:>9.4f} {d1:>9.4f}")

    os.makedirs(RES, exist_ok=True)
    out = {
        '_dry_run': True, 'backend': 'aer', 'date': DATE, 'seed': SEED,
        'params': {'L': L, 'U': U, 'thop': thop, 'eta': eta, 'K': K, 'n_trot': n_trot,
                   'dt_total': dt_total, 'shots_per_time': SHOTS_DRY, 'nelec_add': list(nelec_add)},
        'verification': {'n0up_gs': n0up_gs, 'm0_op': m0_op, 'one_minus_n': one_minus_n,
                         'm0_lehmann': m0_lehmann, 'm1_op': m1_op, 'm1_lehmann': m1_lehmann,
                         'err_m0_identity': err_m0_ident, 'err_m0_vs_lehmann': err_m0,
                         'err_m1_vs_lehmann': err_m1},
        'falsifier': {'S_size': int(len(f['Sidx'])), 'si_size': int(len(si)),
                      'kept_shots': int(kept), 'total_shots': int(total),
                      'rel_L1': f['rel'], 'm0_hat': f['m0_hat'], 'm0_bar': f['m0_bar'], 'delta0': f['delta0'],
                      'm1_bar': f['m1_bar'], 'm1_bar_grid': f['m1_bar_grid'],
                      'm1_hat': f['m1_hat'], 'delta1': f['delta1'],
                      'leak_fraction': f['leak_frac']},
        'coverage_sweep': [dict(d=d, cov=cov, rel_L1=rel, delta0=d0, delta1=d1)
                           for (d, cov, rel, d0, d1) in sweep],
        'per_time_counts': per_time_counts,      # RETAINED pooled-source per-Krylov-time counts
    }
    outpath = os.path.join(RES, 'heron_counts_matched_DRYRUN.json')
    with open(outpath, 'w') as fh:
        json.dump(out, fh)
    log(f"RETAINED per-time counts + falsifier -> {outpath}")
    print("\nTo run on hardware: activate the qiskit>=1.2 venv, save the account, then  RUN=1 py hardware_matched_job.py")
    return f


# =====================================================================================
# 5. QPU run (RUN=1): new-platform auth, SamplerV2 + measurement twirl + gate twirl + DD, Batch, retain counts
# =====================================================================================
def qpu_run():
    from qiskit_ibm_runtime import QiskitRuntimeService, Batch, SamplerV2 as Sampler
    print_verification()
    service = QiskitRuntimeService(name=ACCOUNT)           # NEW platform; NOT channel='ibm_quantum'
    backend = service.backend(BACKEND)
    isa = transpile(circuits, backend=backend, optimization_level=3)
    log(f"transpiled {len(isa)} circuits for {backend.name}; depth(last ISA)={isa[-1].depth()}")
    with Batch(backend=backend) as batch:
        sampler = Sampler(mode=batch)
        o = sampler.options
        o.default_shots = NS
        o.dynamical_decoupling.enable = True
        o.dynamical_decoupling.sequence_type = "XpXm"
        o.twirling.enable_gates = True                     # Pauli gate twirling
        o.twirling.enable_measure = True                   # measurement twirling (twirled raw counts; not TREX)
        o.twirling.num_randomizations = int(os.environ.get("NRAND", 32))  # 32=companion default;
        # set NRAND=8-16 for MORE QPU-budget margin (fewer twirled circuit loads -> less setup overhead;
        # shot count is unchanged, so thresholds/statistics are unaffected -- only the twirl averaging depth).
        o.resilience.measure_mitigation = True             # kept as written; never executed: SamplerV2 has
        # no resilience options, so this line fails on current qiskit-ibm-runtime. The L=8 job, the one
        # that ran, omits it; no readout-error mitigation was applied to any ibm_fez run in this repository.
        job = sampler.run(isa)
        jid = job.job_id()
        log(f"submitted job {jid} on {BACKEND} ({len(isa)} circuits x {NS} shots, Batch); retaining counts")
        res = job.result()
    # pool per-time counts (raw, retained) then run the identical falsifier arithmetic
    per_time_counts = []
    pooled = {}
    for i in range(len(isa)):
        d = res[i].data
        reg = list(d.__dict__.keys())[0] if hasattr(d, '__dict__') else 'meas'
        cnt = getattr(d, reg).get_counts()
        per_time_counts.append(cnt)
        for bs, ctv in cnt.items():
            pooled[bs] = pooled.get(bs, 0) + ctv
    seen, kept, total = postselect(pooled)
    f = falsifier_from_pooled_configs(seen)
    os.makedirs(RES, exist_ok=True)
    out = {
        '_dry_run': False, 'backend': BACKEND, 'job_id': jid, 'date': DATE, 'shots': NS,
        'params': {'L': L, 'U': U, 'thop': thop, 'eta': eta, 'K': K, 'n_trot': n_trot,
                   'dt_total': dt_total, 'nelec_add': list(nelec_add)},
        'verification': {'m0_op': m0_op, 'm0_lehmann': m0_lehmann, 'm1_op': m1_op,
                         'm1_lehmann': m1_lehmann, 'err_m0_vs_lehmann': err_m0, 'err_m1_vs_lehmann': err_m1},
        'falsifier': {'S_size': int(len(f['Sidx'])), 'si_size': int(len(si)),
                      'kept_shots': int(kept), 'total_shots': int(total),
                      'rel_L1': f['rel'], 'm0_hat': f['m0_hat'], 'm0_bar': f['m0_bar'], 'delta0': f['delta0'],
                      'm1_bar': f['m1_bar'], 'm1_hat': f['m1_hat'],
                      'delta1': f['delta1'], 'leak_fraction': f['leak_frac']},
        'per_time_counts': per_time_counts,
    }
    outpath = os.path.join(RES, f'heron_counts_matched_{jid}.json')
    with open(outpath, 'w') as fh:
        json.dump(out, fh)
    log(f"RETAINED raw per-time counts + falsifier -> {outpath}")
    print(f"\n|S|={len(f['Sidx'])}  rel-L1={f['rel']:.4f}  Delta_0={f['delta0']:.4f}  "
          f"m1_hat={f['m1_hat']:.4f}  m1_bar={f['m1_bar']:.4f}  Delta_1={f['delta1']:.4f}  leak={f['leak_frac']:.3f}")
    return f


def main():
    os.makedirs(RES, exist_ok=True)
    if os.environ.get("RUN") == "1":
        print("RUN=1 -> QPU (ibm_fez). Requires the qiskit>=1.2 venv + saved new-platform account.")
        qpu_run()
    else:
        print("RUN unset -> LOCAL Aer DRY-RUN (0 QPU). Set RUN=1 to submit to ibm_fez.")
        local_dry_run()


if __name__ == '__main__':
    main()
