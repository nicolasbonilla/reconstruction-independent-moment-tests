# -*- coding: utf-8 -*-
r"""
MATCHED, RECONSTRUCTION-INDEPENDENT DEVICE-SIDE MOMENT FALSIFIER  (Paper 3) -- L=8 (16 qubits)
=============================================================================================
This is the L=8 sibling of hardware_matched_job.py (L=6). SAME physics, SAME circuit family,
SAME post-selection, SAME reconstruction-independent falsifier arithmetic; only the system size
(and hence the prep occupation) changes:

    L = 8, U/t = 4, thop = 1.0, eta = 0.15, M = 16 spin-orbitals (16 qubits), interleaved
    spin-orbital order p = 2*site + spin (even = up, odd = down).

Reference state (fixed, classically known, SHARED by both branches):
    |phi> = c^dagger_{site0,up} |0>,   |0> = exact half-filled (Nocc==L==8, Szocc==0) GS, energy E0.
|phi> lives entirely in the ADDITION sector si = {Nocc==L+1==9 & Szocc==+1}.

Device prep (number/Sz conserving shallow Trotter): nelec_add = (5 up, 4 down) = 9 electrons,
K=7 Krylov-Trotter circuits (n_trot=3, dt_total=0.5/thop), pooled counts, post-select
Nocc==9 & Szocc==+1 -> subspace S, then the SAME falsifier:

    m0^+ = <phi|phi> = 1 - <n_{0up}>_GS          (ADDITION-only zeroth moment; NOT the m0=1 anticommutator)
    m1^+ = <phi|(H-E0)|phi>                       (ADDITION-only first moment)
    m0_bar = within-S weight ||P_S phi||^2,  m1_bar = within-S first moment
    Delta_0 = |m0_hat - m0_bar|  (MONOTONE completeness probe; m0_hat = m0^+, exact, reconstruction-indep)
    Delta_1 = |m1_hat - m1_bar|  (m1_hat = m1^+, full H on full phi, reconstruction-indep)

WHY L=8 (the strategic point):
    At L=6 the device-sampled subspace covers the ENTIRE addition sector (300/300 configs), so the
    screen returns Delta_0 = 0 -- it only proves executability. At L=8 the addition sector has
    |si| = 3920 configs and the same K=7 shallow-Trotter prep, sampled at finite shots, genuinely
    LEAKS: noiseless Aer at 50k shots gives Delta_0 = 0.0238 (leak 4.76%); on the real device the
    noise BROADENS the sampled support so Delta_0 = 0.0044 at 50k (94% coverage), rising to 0.3726
    at 4k (35% coverage). The sampled S omits real addition-sector weight, so Delta_0 FIRES while the
    3920-dim sector is still EXACTLY diagonalizable -- coverage calibration, validated not asserted.
    L=10 saturates the leak (avoid); L=8 is the sweet spot at the same ~2-3 min budget.

KEEPS the correct int(bs, 2) bitstring->Fock mapping (rightmost char = qubit 0), NOT the reversed
companion bug, and the ADDITION-only moments (m0^+ = 1 - <n_0up>, NOT the two-sided m0 = 1).

USAGE (LOCAL, 0 QPU):
    C:\tmp\ibmqpu\Scripts\python.exe hardware_matched_job_L8.py
Optional env overrides: L=8 NELEC_UP=5 NELEC_DN=4 SHOTS_DRY=50000 SEED=1.
"""
import os, sys, json, time
import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, '..', '06_results'))
DATE = '2026-08-28'

_trapz = getattr(np, 'trapezoid', None) or np.trapz

NS = int(os.environ.get("SHOTS", 50000))   # real-QPU shots/circuit (SHOTS=4000 -> low-shot firing run)
BACKEND = os.environ.get("IBM_BACKEND", "ibm_fez")
ACCOUNT = os.environ.get("IBM_ACCOUNT", "open-instance")
SHOTS_DRY = int(os.environ.get("SHOTS_DRY", NS))
SEED = int(os.environ.get("SEED", 1))

t0 = time.time()
log = lambda *a: print(f"[{time.time()-t0:6.1f}s]", *a, flush=True)

# =====================================================================================
# 1. Hamiltonian + exact addition-sector reference  (L=8)
# =====================================================================================
L = int(os.environ.get("L", 8)); U = 4.0; thop = 1.0; eta = 0.15; M = 2 * L
N_UP = int(os.environ.get("NELEC_UP", (L // 2) + 1))   # addition sector: 5 up  (default L=8)
N_DN = int(os.environ.get("NELEC_DN", (L // 2)))       #                  4 down
nelec_add = (N_UP, N_DN)

def c_op(p, dim):
    r = []; c = []; d = []
    for s in range(dim):
        if (s >> p) & 1:
            sg = (-1) ** bin(s & ((1 << p) - 1)).count('1'); r.append(s & ~(1 << p)); c.append(s); d.append(float(sg))
    return sp.csr_matrix((d, (r, c)), shape=(dim, dim))

dim = 1 << M
log(f"L={L}  M={M} spin-orbitals  Fock dim=2^{M}={dim}  building c-operators + H ...")
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
log(f"H built (nnz={H.nnz}); diagonalizing half-filled GS sector (Nocc=={L}, Szocc==0) ...")

gi = np.where((Nocc == L) & (Szocc == 0))[0]
w, v = np.linalg.eigh(H[gi][:, gi].toarray()); E0 = w[0]
psi0 = np.zeros(dim, complex); psi0[gi] = v[:, 0]
p_orb = 0
phi = Cd[2 * p_orb] @ psi0                             # |phi> = c^dagger_{site0,up} |0>
si = np.where((Nocc == L + 1) & (Szocc == 1))[0]       # ADDITION sector (N+1, Sz=+1)
log(f"GS sector dim={len(gi)}  E0={E0:.6f};  addition sector |si|={len(si)}; diagonalizing H_si ({len(si)}x{len(si)}) ...")
En, Vn = np.linalg.eigh(H[si][:, si].toarray()); coef = Vn.conj().T @ phi[si]
grid = np.linspace((En - E0).min() - 1, (En - E0).max() + 1, 600)

def spec(pw, ww):
    A = np.zeros_like(grid)
    for a, b in zip(pw, ww):
        A += b * (eta / np.pi) / ((grid - a) ** 2 + eta ** 2)
    return A

A_exact = spec(En - E0, np.abs(coef) ** 2)

# ---- exact ADDITION-sector moments and their <=1e-10 verification ----
m0_lehmann = float(np.sum(np.abs(coef) ** 2))
m1_lehmann = float(np.sum((En - E0) * np.abs(coef) ** 2))
n0up_gs = float(np.real(np.vdot(psi0, (Cd[0] @ C[0]) @ psi0)))
m0_op = float(np.real(np.vdot(phi, phi)))
m1_op = float(np.real(np.vdot(phi, H @ phi)) - E0 * np.real(np.vdot(phi, phi)))
one_minus_n = 1.0 - n0up_gs
err_m0_ident = abs(m0_op - one_minus_n)
err_m0 = abs(m0_op - m0_lehmann)
err_m1 = abs(m1_op - m1_lehmann)
gm0 = float(_trapz(A_exact, grid)); gm1 = float(_trapz(grid * A_exact, grid))


# =====================================================================================
# 2. Companion circuits (L=8 prep: 5 up on even qubits, 4 down on odd qubits)
# =====================================================================================
from qiskit import QuantumCircuit, QuantumRegister, transpile
norb = L
K = 7; n_trot = 3; dt_total = 0.5 / thop

def givens(qc, a, b, theta):      # e^{-i theta (XX+YY)/2} hopping between spin-orbitals a,b
    qc.rxx(theta, a, b); qc.ryy(theta, a, b)

def build_circuit(k):
    q = QuantumRegister(2 * L); qc = QuantumCircuit(q)
    for orb in range(N_UP): qc.x(q[2 * orb])           # N_UP up electrons (even = up)
    for orb in range(N_DN): qc.x(q[2 * orb + 1])       # N_DN down electrons (odd = down)
    dt = (k * dt_total) / max(n_trot, 1)
    for _ in range(n_trot):
        for i in range(L): qc.rzz(2 * U * dt, q[2 * i], q[2 * i + 1])         # on-site U (native RZZ)
        for spin in (0, 1):                                                  # hopping per spin
            for i in range(L):
                j = (i + 1) % L; givens(qc, q[2 * i + spin], q[2 * j + spin], thop * dt)
    qc.measure_all(); return qc

circuits = [build_circuit(k) for k in range(K)]

def bitstring_to_config(bs):
    r"""qiskit little-endian counts bitstring -> Fock integer in the OPERATOR convention.
    CORRECT map is int(bs, 2) (rightmost char = qubit 0), NOT int(bs[::-1], 2) (the reversed
    companion bug, which flips even<->odd = every spin label and post-selects the empty set)."""
    return int(bs.replace(' ', ''), 2)


# =====================================================================================
# 3. Shared falsifier arithmetic
# =====================================================================================
pos = {int(x): i for i, x in enumerate(si)}

def falsifier_from_pooled_configs(seen):
    S = np.array(sorted(seen))
    Sidx = np.array([pos[x] for x in S if x in pos], dtype=int)
    HS = H[si][:, si].toarray()[np.ix_(Sidx, Sidx)]
    Em, Um = np.linalg.eigh(HS)
    aS = Um.conj().T @ phi[si][Sidx]
    A_hw = spec(Em - E0, np.abs(aS) ** 2)
    rel = float(_trapz(np.abs(A_hw - A_exact), grid) / _trapz(A_exact, grid))
    m1_bar = float(np.sum((Em - E0) * np.abs(aS) ** 2))
    m1_bar_grid = float(_trapz(grid * A_hw, grid))
    m0_bar = float(np.sum(np.abs(aS) ** 2))
    m1_hat = m1_op
    delta1 = abs(m1_hat - m1_bar)
    delta0 = abs(m0_op - m0_bar)
    leak_frac = max(0.0, 1.0 - m0_bar / m0_op)
    return dict(S=S, Sidx=Sidx, A_hw=A_hw, rel=rel, m1_bar=m1_bar, m1_bar_grid=m1_bar_grid,
                m0_bar=m0_bar, m0_hat=m0_op, delta0=delta0, m1_hat=m1_hat, delta1=delta1, leak_frac=leak_frac)


def coverage_sweep():
    w_si = np.abs(phi[si]) ** 2
    order = np.argsort(w_si)[::-1]
    rows = []
    for frac in (1.0, 0.9, 0.75, 0.5, 0.25, 0.1):
        d = max(1, int(round(frac * len(si))))
        seen = set(int(si[j]) for j in order[:d])
        f = falsifier_from_pooled_configs(seen)
        rows.append((d, f['m0_bar'] / m0_op, f['rel'], f['delta0'], f['delta1']))
    return rows


def postselect(pooled_counts):
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
    print(f"    {'d':>5} {'cov(m0_bar/m0)':>14} {'rel-L1':>8} {'Delta_0':>9} {'Delta_1':>9}")
    for d, cov, rel, d0, d1 in sweep:
        print(f"    {d:>5} {cov:>14.4f} {rel:>8.4f} {d0:>9.4f} {d1:>9.4f}")

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
        'per_time_counts': per_time_counts,
    }
    outpath = os.path.join(RES, 'heron_counts_matched_L8_DRYRUN.json')
    with open(outpath, 'w') as fh:
        json.dump(out, fh)
    log(f"RETAINED per-time counts + falsifier -> {outpath}")
    return f


# =====================================================================================
# 5. QPU run (RUN=1): identical to L=6 job (new-platform auth, SamplerV2 + TREX + twirl + DD)
# =====================================================================================
def qpu_run():
    from qiskit_ibm_runtime import QiskitRuntimeService, Batch, SamplerV2 as Sampler
    print_verification()
    service = QiskitRuntimeService(name=ACCOUNT)
    backend = service.backend(BACKEND)
    # translation_method='translator' bypasses the backend's 'ibm_dynamic_circuits' translation
    # plugin (unavailable in this qiskit/runtime pair); our circuits are static, so the standard
    # BasisTranslator is correct and sufficient.
    isa = transpile(circuits, backend=backend, optimization_level=3, translation_method='translator')
    log(f"transpiled {len(isa)} circuits for {backend.name}; depth(last ISA)={isa[-1].depth()}")
    with Batch(backend=backend) as batch:
        sampler = Sampler(mode=batch)
        o = sampler.options
        o.default_shots = NS
        o.dynamical_decoupling.enable = True
        o.dynamical_decoupling.sequence_type = "XpXm"
        o.twirling.enable_gates = True
        o.twirling.enable_measure = True
        o.twirling.num_randomizations = int(os.environ.get("NRAND", 32))
        # measurement-error mitigation = TREX, already enabled above via twirling.enable_measure;
        # SamplerV2 in qiskit-ibm-runtime 0.44 has no o.resilience attribute (that was the old API).
        job = sampler.run(isa)
        jid = job.job_id()
        log(f"submitted job {jid} on {BACKEND} ({len(isa)} circuits x {NS} shots, Batch); retaining counts")
        res = job.result()
    per_time_counts = []
    pooled = {}
    for i in range(len(isa)):
        d = res[i].data
        ba = getattr(d, 'meas', None)                      # measure_all() names the register 'meas'
        if ba is None:                                     # robust fallback: single measured register
            ba = getattr(d, list(d.keys())[0])
        cnt = ba.get_counts()
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
    outpath = os.path.join(RES, f'heron_counts_matched_L8_{jid}.json')
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
