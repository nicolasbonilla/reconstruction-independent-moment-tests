# -*- coding: utf-8 -*-
r"""
Re-analysis of the two-arm ibm_fez run of 2026-08-21 (jobs da44t8u1vhnc73fkie30 and da459juaa69c739jcb80).

WHAT WAS RUN. Two SamplerV2 jobs on ibm_fez, submitted from the Colab notebook GATE_TP_on_Heron.ipynb
(work/notebooks/ of this project): L=6 Hubbard ring, 12 qubits (qubit 2i = site i spin up, 2i+1 = spin down),
an X-gate product determinant (4 up, 3 down) followed by n_trot=3 layers of on-site R_zz(2 U dt) and
nearest-neighbour R_xx R_yy hopping (periodic wrap, no Jordan-Wigner strings), K=7 circuits with
dt = k * 0.5 / 3, k = 0..6, 35 000 shots per circuit, twirling (gates + measurement) and dynamical
decoupling on. Arm "correct": U=4.0 in the R_zz angle (da44t8u1vhnc73fkie30). Arm "detuned": U=3.5
(da459juaa69c739jcb80). The raw results, the submitted (transpiled) circuits and the usage metrics are
deposited unchanged in data/heron_gatetp_2026-08-21/ (no credential in them; user_id is an account id).

WHAT THIS SCRIPT DOES
  1. Decodes the deposited raw results with numpy only (base64 -> zlib -> .npy BitArray, big-endian bytes,
     rightmost character = clbit 0), keeping the first-appearance order of qiskit's BitArray.get_counts.
     If qiskit_ibm_runtime is importable, the decode is cross-checked against its RuntimeDecoder.
  2. Reproduces the notebook's analysis pipeline EXACTLY (cells 2-4, 6, 8, 13, 15, 17, 19, 20, 22; same
     functions, same rng = default_rng(0) consumed in the same order, the local rehearsal of cell 8
     included because it draws from that rng), once per bit order:
        as_run  : bits_to_cfg(bs) = int(bs[::-1], 2)   (the notebook, cell 3)
        correct : bits_to_cfg(bs) = int(bs, 2)         (qiskit's order: rightmost char = qubit 0)
     The as-run outputs are compared with the values printed by the executed notebook on 2026-08-21.
  3. Verifies the bit order the way src/delta0_reference_mc.py does, with a noiseless simulation of the
     same circuits: Statevector index i <-> bitstring format(i,'012b') <-> config int(bs,2); the in-sector
     probability of every circuit is 1 under int(bs,2) and 0 under int(bs[::-1],2). If
     qiskit_ibm_runtime is importable, the SUBMITTED transpiled circuits (decoded from the job JSON,
     reduced to their 12 active physical qubits, measured clbits read through their measure map) are
     simulated too and compared with the notebook circuits at U=4 and U=3.5.
  4. Determines what the screen tests: the reconstruction diagonalises H(U_rec) on the support S and
     uses the EXACT ED vector phi(U_rec) and E0(U_rec); the device enters only through S. When S is the
     whole sector, mbar(U_rec) = exact moments of H(U_rec), so the screen residual
     |mbar - m(U=4)| / m(U=4) is a function of U_rec alone. Checked by (a) a data-free control
     (S = the full sector, no device data), (b) a 2x2 cross-assignment (each arm's support reconstructed
     with U_rec = 4 and 3.5) under both bit orders, (c) a coverage scan at small retained-shot numbers
     (where S is incomplete and the data do enter), and (d) a direct, non-moment check of whether the
     device counts themselves distinguish U=4 from U=3.5 (total variation and log-likelihood ratio of the
     retained in-sector counts against the noiseless distributions of both schedules).

Output: data/2026-10-02_gatetp_reanalysis.json (new file).
USAGE (0 QPU, 9-12 min on one BLAS thread under load):  python src/gate_tp_analysis.py
Needs numpy, scipy, qiskit (Statevector, circuit building); qiskit_ibm_runtime is optional (cross-checks).
"""
import os, sys, json, time, base64, zlib, io, hashlib, platform
import numpy as np
import scipy.sparse as sp

T_START = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "data"))
RAW = os.path.join(RES, "heron_gatetp_2026-08-21")
OUT = os.environ.get("OUT_JSON", os.path.join(RES, "2026-10-02_gatetp_reanalysis.json"))
JOBS = {"correct": ("da44t8u1vhnc73fkie30", 4.0), "detuned": ("da459juaa69c739jcb80", 3.5)}
SEED_SCAN = 20261002
log = lambda *a: print(f"[{time.time() - T_START:7.1f}s]", *a, flush=True)


def sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


# =============================================================================================
# Notebook cells 2-4, verbatim apart from plotting (GATE_TP_on_Heron.ipynb)
# =============================================================================================
if not hasattr(np, 'trapezoid'):                       # notebook: np.trapezoid = getattr(np,'trapezoid',np.trapz);
    np.trapezoid = np.trapz                            # same effect, but numpy >= 2.4 has no np.trapz attribute
L = 6; thop = 1.0; eta = 0.15; M = 2 * L; U_TRUE = 4.0


def c_op(p, dim):
    r = []; cc = []; d = []
    for s in range(dim):
        if (s >> p) & 1:
            sg = (-1) ** bin(s & ((1 << p) - 1)).count('1'); r.append(s & ~(1 << p)); cc.append(s); d.append(float(sg))
    return sp.csr_matrix((d, (r, cc)), shape=(dim, dim))


dim = 1 << M; Cc = [c_op(p, dim) for p in range(M)]; Cd = [c.T.conj() for c in Cc]
Nocc = np.array([bin(s).count('1') for s in range(dim)])
Szocc = np.array([sum(((s >> (2 * i)) & 1) - ((s >> (2 * i + 1)) & 1) for i in range(L)) for s in range(dim)])
Hhop = sp.csr_matrix((dim, dim))
for i in range(L):
    j = (i + 1) % L
    for spin in (0, 1):
        a = 2 * i + spin; b = 2 * j + spin; Hhop = Hhop - thop * (Cd[a] @ Cc[b] + Cd[b] @ Cc[a])
Hint = sp.csr_matrix((dim, dim))
for i in range(L): Hint = Hint + (Cd[2 * i] @ Cc[2 * i]) @ (Cd[2 * i + 1] @ Cc[2 * i + 1])
gi = np.where((Nocc == L) & (Szocc == 0))[0]; si = np.where((Nocc == L + 1) & (Szocc == 1))[0]
pos = {int(x): k for k, x in enumerate(si)}
GRID = np.linspace(-1.0, 12.0, 600)


def spec(om, wt):
    A = np.zeros_like(GRID)
    for a, b in zip(om, wt): A += b * (eta / np.pi) / ((GRID - a) ** 2 + eta ** 2)
    return A


_cache = {}


def prep(U):
    if U in _cache: return _cache[U]
    H = (Hhop + U * Hint).tocsr()
    w, v = np.linalg.eigh(H[gi][:, gi].toarray()); E0 = w[0]
    psi0 = np.zeros(dim, complex); psi0[gi] = v[:, 0]
    phi = (Cd[0] @ psi0)[si]
    Hs = H[si][:, si].toarray()
    En, Vn = np.linalg.eigh(Hs); coef = Vn.conj().T @ phi
    om = En - E0; wt = np.abs(coef) ** 2
    m = np.array([wt.sum(), (om * wt).sum(), (om ** 2 * wt).sum()])
    d = dict(U=U, E0=E0, Hs=Hs, phi=phi, A_exact=spec(om, wt), m_ind=m); _cache[U] = d; return d


P4 = prep(U_TRUE)


def reconstruct(configs, U_rec):
    Pr = prep(U_rec)
    Sidx = np.array(sorted({pos[int(x)] for x in configs if int(x) in pos}), dtype=int)
    if len(Sidx) == 0: return None
    Em, Um = np.linalg.eigh(Pr['Hs'][np.ix_(Sidx, Sidx)])
    aS = Um.conj().T @ Pr['phi'][Sidx]; wS = np.abs(aS) ** 2; omS = Em - Pr['E0']
    A = spec(omS, wS); mbar = np.array([wS.sum(), (omS * wS).sum(), (omS ** 2 * wS).sum()])
    return dict(A=A, mbar=mbar, nS=int(len(Sidx)))


def screen_resid(mbar):
    return np.abs(mbar - P4['m_ind']) / (np.abs(P4['m_ind']) + 1e-12)


def spectral_distance(A):
    return float(np.trapezoid(np.abs(A - P4['A_exact']), GRID) / np.trapezoid(P4['A_exact'], GRID))


def haydock_A(configs, U_rec, nlanc=200):
    Pr = prep(U_rec)
    Sidx = np.array(sorted({pos[int(x)] for x in configs if int(x) in pos}), dtype=int)
    Hs = Pr['Hs'][np.ix_(Sidx, Sidx)]; f = Pr['phi'][Sidx].astype(complex)
    n0 = np.linalg.norm(f)
    if n0 < 1e-14: return np.zeros_like(GRID)
    q = f / n0; a = []; b = []; qprev = np.zeros_like(q); bprev = 0.0
    for _ in range(min(nlanc, len(Sidx))):
        Hq = Hs @ q; al = np.real(np.vdot(q, Hq)); a.append(al)
        r = Hq - al * q - bprev * qprev; bn = np.linalg.norm(r)
        if bn < 1e-10: break
        qprev = q.copy(); q = r / bn; b.append(bn); bprev = bn
    a = np.array(a); b = np.array(b)
    z = (GRID + Pr['E0']) + 1j * eta

    def cf_eval(zz):
        cf = 0.0
        for k in range(len(a) - 1, -1, -1):
            num = (b[k - 1] ** 2) if k > 0 else 1.0
            cf = num / (zz - a[k] - cf)
        return cf
    return np.array([(n0 ** 2) * (-1 / np.pi) * np.imag(cf_eval(zz)) for zz in z])


MAPS = {"as_run": lambda bs: int(bs[::-1], 2),      # notebook cell 3: def bits_to_cfg(bs): return int(bs[::-1],2)
        "correct": lambda bs: int(bs, 2)}           # qiskit order: rightmost character = qubit 0

from qiskit import QuantumCircuit, QuantumRegister
from qiskit.quantum_info import Statevector
K = 7; n_trot = 3; dt_total = 0.5 / thop


def givens(qc, a, b, th): qc.rxx(th, a, b); qc.ryy(th, a, b)


def build_circuit(k, U):
    q = QuantumRegister(2 * L); qc = QuantumCircuit(q)
    for orb in [0, 1, 2, 3]: qc.x(q[2 * orb])
    for orb in [0, 1, 2]: qc.x(q[2 * orb + 1])
    dt = (k * dt_total) / max(n_trot, 1)
    for _ in range(n_trot):
        for i in range(L): qc.rzz(2 * U * dt, q[2 * i], q[2 * i + 1])
        for spin in (0, 1):
            for i in range(L):
                j = (i + 1) % L; givens(qc, q[2 * i + spin], q[2 * j + spin], thop * dt)
    qc.measure_all(); return qc


DETUNE_U = 3.5
arms = {U_TRUE: 'correct', DETUNE_U: 'detuned'}
SHOTS = 35000


# =============================================================================================
# 1. Raw results -> per-circuit counts (numpy-only decode, first-appearance order)
# =============================================================================================
def decode_counts(path):
    d = json.load(open(path, encoding="utf-8"))
    out = []
    for pr in d["__value__"]["pub_results"]:
        meas = pr["__value__"]["data"]["__value__"]["fields"]["meas"]["__value__"]
        nb = int(meas["num_bits"])
        arr = np.load(io.BytesIO(zlib.decompress(base64.b64decode(meas["array"]["__value__"]))), allow_pickle=False)
        ints = [int.from_bytes(row.tobytes(), "big") & ((1 << nb) - 1) for row in arr]
        c = {}
        for x in ints:
            s = format(x, f"0{nb}b"); c[s] = c.get(s, 0) + 1
        out.append(c)
    return out


raw_files = {}
counts = {}
for arm, (jid, U) in JOBS.items():
    for kind in ("job", "results", "metrics"):
        p = os.path.join(RAW, f"{jid}_{kind}.json")
        raw_files[f"{jid}_{kind}.json"] = sha256(p)
    counts[U] = decode_counts(os.path.join(RAW, f"{jid}_results.json"))
    assert len(counts[U]) == K and all(sum(c.values()) == SHOTS for c in counts[U])
log("decoded raw results:", {arms[U]: [len(c) for c in counts[U]] for U in counts})

runtime_crosscheck = {"available": False}
isa_check = {"available": False}
try:
    from qiskit_ibm_runtime import RuntimeDecoder  # optional
    runtime_crosscheck["available"] = True
    same = True
    jobmeta = {}
    for arm, (jid, U) in JOBS.items():
        r = json.load(open(os.path.join(RAW, f"{jid}_results.json"), encoding="utf-8"), cls=RuntimeDecoder)
        for k in range(K):
            ck = r[k].data.meas.get_counts()
            same &= (list(ck.items()) == list(counts[U][k].items()))
    runtime_crosscheck["numpy_decode_identical_to_RuntimeDecoder_get_counts_incl_order"] = bool(same)
    assert same
except ImportError:
    RuntimeDecoder = None

metrics = {}
jobinfo = {}
for arm, (jid, U) in JOBS.items():
    m = json.load(open(os.path.join(RAW, f"{jid}_metrics.json"), encoding="utf-8"))
    jj = json.load(open(os.path.join(RAW, f"{jid}_job.json"), encoding="utf-8"))
    metrics[arm] = dict(quantum_seconds=m["usage"]["quantum_seconds"], created=m["timestamps"]["created"],
                        running=m["timestamps"]["running"], finished=m["timestamps"]["finished"],
                        caller_versions=m.get("qiskit_version"))
    jobinfo[arm] = dict(job_id=jj["id"], backend=jj["backend"], program=jj["program"]["id"], status=jj["status"],
                        options=jj["params"].get("options"), n_pubs=len(jj["params"]["pubs"]),
                        shots_per_pub=[pb[2] for pb in jj["params"]["pubs"]])

# =============================================================================================
# 2. Bit-order verification with a noiseless simulation of the same circuits
# =============================================================================================
P = {U: np.array([np.abs(Statevector(build_circuit(k, U).remove_final_measurements(inplace=False)).data) ** 2
                  for k in range(K)]) for U in (U_TRUE, DETUNE_U)}
all_bs = [format(i, f"0{M}b") for i in range(dim)]
insec = np.zeros(dim, bool); insec[si] = True
rev_idx = np.array([MAPS["as_run"](bs) for bs in all_bs])
assert all(MAPS["correct"](all_bs[i]) == i for i in range(dim))
noiseless = {}
for U in P:
    noiseless[arms[U]] = dict(
        insector_mass_correct_map=[float(P[U][k][insec].sum()) for k in range(K)],
        insector_mass_as_run_map=[float(P[U][k][insec[rev_idx]].sum()) for k in range(K)],
        support_size_per_circuit=[int((P[U][k] > 1e-12).sum()) for k in range(K)],
        union_support_size=int((P[U] > 1e-12).any(axis=0).sum()))
assert all(abs(x - 1) < 1e-10 for U in P for x in noiseless[arms[U]]["insector_mass_correct_map"])
assert all(x < 1e-20 for U in P for x in noiseless[arms[U]]["insector_mass_as_run_map"])
# mirror: int(bs[::-1],2) sends qubit p to bit 11-p, i.e. up<->down and site i <-> site 5-i: Sz -> -Sz
mirror_sector = np.where((Nocc == L + 1) & (Szocc == -1))[0]
assert set(int(x) for x in rev_idx[si]) == set(int(x) for x in mirror_sector)
log("noiseless: in-sector mass = 1 (int(bs,2)), 0 (int(bs[::-1],2)) for all 14 circuits")

# submitted ISA circuits (optional): reduce to the active qubits, simulate, read clbits through measure map
if RuntimeDecoder is not None:
    isa_check = {"available": True, "per_arm": {}}
    for arm, (jid, Uarm) in JOBS.items():
        jj = json.load(open(os.path.join(RAW, f"{jid}_job.json"), encoding="utf-8"), cls=RuntimeDecoder)
        rows = []
        for k, pub in enumerate(jj["params"]["pubs"]):
            isa = pub[0]
            act = []
            for ins in isa.data:
                if ins.operation.name == "barrier": continue
                for qb in ins.qubits:
                    qi = isa.find_bit(qb).index
                    if qi not in act: act.append(qi)
            red = QuantumCircuit(len(act)); mmap = {}
            for ins in isa.data:
                nm = ins.operation.name
                if nm == "barrier": continue
                qs = [act.index(isa.find_bit(qb).index) for qb in ins.qubits]
                if nm == "measure":
                    mmap[isa.find_bit(ins.clbits[0]).index] = qs[0]; continue
                red.append(ins.operation, qs)
            assert sorted(mmap) == list(range(M)) and len(act) == M
            pr = np.abs(Statevector(red).data) ** 2              # index over reduced qubits
            # probability of clbit integer x: bit c of x = bit mmap[c] of reduced index
            ridx = np.arange(dim)
            xint = np.zeros(dim, dtype=np.int64)
            for c, rq in mmap.items():
                xint |= ((ridx >> rq) & 1) << c
            px = np.zeros(dim); np.add.at(px, xint, pr)
            tv4 = 0.5 * float(np.abs(px - P[U_TRUE][k]).sum()); tv35 = 0.5 * float(np.abs(px - P[DETUNE_U][k]).sum())
            rows.append(dict(k=k, n_active_qubits=len(act), physical_qubits_by_clbit=[act[mmap[c]] for c in range(M)],
                             tv_vs_notebook_U4=tv4, tv_vs_notebook_U3p5=tv35,
                             insector_mass_correct_map=float(px[insec].sum()),
                             insector_mass_as_run_map=float(px[insec[rev_idx]].sum())))
        isa_check["per_arm"][arm] = dict(
            rows=rows,
            max_tv_vs_own_schedule=max(r["tv_vs_notebook_U4" if Uarm == 4.0 else "tv_vs_notebook_U3p5"] for r in rows),
            min_tv_vs_other_schedule_k_ge_1=min(r["tv_vs_notebook_U3p5" if Uarm == 4.0 else "tv_vs_notebook_U4"]
                                                for r in rows if r["k"] >= 1))
    isa_check["note"] = ("the transpiled circuits decoded from each job JSON, reduced to their 12 active physical "
                         "qubits and read through their measure instructions, reproduce the notebook circuits "
                         "of their own schedule (TV below 1e-13) and differ from the other schedule for k >= 1; "
                         "k = 0 is the bare determinant (no U dependence)")
    log("ISA check:", {a: (v["max_tv_vs_own_schedule"], v["min_tv_vs_other_schedule_k_ge_1"])
                       for a, v in isa_check["per_arm"].items()})

# device: retention under both maps, support agreement, Sz=-1 explanation of the as-run retention
device_bitorder = {}
for U in (U_TRUE, DETUNE_U):
    rows = []
    for k in range(K):
        c = counts[U][k]
        xs = np.array([int(bs, 2) for bs in c]); ns = np.array(list(c.values()))
        kept_ok = int(ns[insec[xs]].sum()); kept_rev = int(ns[insec[rev_idx[xs]]].sum())
        mirror = int(ns[(Nocc[xs] == L + 1) & (Szocc[xs] == -1)].sum())
        on_supp = int(ns[insec[xs] & (P[U][k][xs] > 1e-12)].sum())
        rows.append(dict(k=k, raw=int(ns.sum()), retained_correct=kept_ok, retained_as_run=kept_rev,
                         shots_with_N7_Szm1_correct_reading=mirror,
                         retained_correct_on_noiseless_support=on_supp))
        assert kept_rev == mirror
    tot = lambda key: int(sum(r[key] for r in rows))
    device_bitorder[arms[U]] = dict(per_circuit=rows, raw=tot("raw"), retained_correct=tot("retained_correct"),
                                    retained_as_run=tot("retained_as_run"),
                                    retention_correct=tot("retained_correct") / tot("raw"),
                                    retention_as_run=tot("retained_as_run") / tot("raw"),
                                    retained_correct_on_noiseless_support=tot("retained_correct_on_noiseless_support"))
log("retention:", {a: (v["retained_correct"], v["retained_as_run"]) for a, v in device_bitorder.items()})


# =============================================================================================
# 3. The notebook pipeline, reproduced exactly, once per bit order
# =============================================================================================
def run_notebook_pipeline(bits_to_cfg):
    rng = np.random.default_rng(0)                         # cell 2
    out = {}
    # cell 6 (pre-flight gate; no rng)
    gate = []
    for Ud in [3.0, DETUNE_U, 4.5, 5.0]:
        r = np.abs(prep(Ud)['m_ind'] - P4['m_ind']) / np.abs(P4['m_ind'])
        gate.append(dict(U_det=Ud, resid_pct=(100 * r).tolist()))
    out["cell6_preflight_gate"] = gate
    # cell 8 (local rehearsal; consumes rng exactly as in the notebook)
    def sample_arm_local(U, shots):
        o = []
        for k in range(K):
            sv = Statevector(build_circuit(k, U).remove_final_measurements(inplace=False))
            p = np.abs(sv.data) ** 2; idx = rng.choice(len(p), size=shots, p=p / p.sum())
            o += [int(t) for t in idx if Nocc[t] == L + 1 and Szocc[t] == 1]
        return np.array(o)
    demo = {U_TRUE: sample_arm_local(U_TRUE, 35000), DETUNE_U: sample_arm_local(DETUNE_U, 35000)}
    out["cell8_local_rehearsal"] = {}
    for U in demo:
        v = demo[U]; r = reconstruct(v, U); res = screen_resid(r['mbar'])
        out["cell8_local_rehearsal"][arms[U]] = dict(retained=int(len(v)), nS=r['nS'],
                                                     spectral_dist=spectral_distance(r['A']),
                                                     screen_pct=(100 * res).tolist())
    # cell 12 (bits list as the notebook built it)
    raw = {U: dict(bits=[bs for c in counts[U] for bs, ct in c.items() for _ in range(ct)]) for U in arms}
    # cell 13
    def to_configs(bits, clean=True):
        cfg = np.array([bits_to_cfg(b) for b in bits])
        if clean: cfg = cfg[(Nocc[cfg] == L + 1) & (Szocc[cfg] == 1)]
        return cfg
    data = {}
    out["cell13_retention"] = {}
    for U in arms:
        cln = to_configs(raw[U]['bits'], clean=True); rw = to_configs(raw[U]['bits'], clean=False)
        rng.shuffle(cln)
        data[U] = dict(clean=cln, raw=rw)
        out["cell13_retention"][arms[U]] = dict(total=len(raw[U]['bits']), cleaned_sector=int(len(cln)),
                                                unique=int(len(np.unique(cln))), sector_dim=int(len(si)))
    # cell 15
    vc = data[U_TRUE]['clean']; full = reconstruct(vc, U_TRUE); res1 = screen_resid(full['mbar'])
    sd1 = spectral_distance(full['A'])
    B = 200; m2b = np.array([reconstruct(vc[rng.integers(0, len(vc), len(vc))], U_TRUE)['mbar'][2] for _ in range(B)])
    lo, hi = np.percentile(m2b, [16, 84]); inside = lo <= P4['m_ind'][2] <= hi
    out["cell15_result1"] = dict(nS=full['nS'], spectral_dist=sd1, screen_pct=(100 * res1).tolist(),
                                 verdict='PASS' if sd1 < 0.05 else 'FIRES', m2_hw=float(full['mbar'][2]),
                                 m2_ED=float(P4['m_ind'][2]), m2_bias=float(full['mbar'][2] - P4['m_ind'][2]),
                                 band16_84=[float(lo), float(hi)], ED_inside_band=bool(inside))
    # cell 17
    B = 200
    pts = []
    for U in arms:
        v = data[U]['clean']
        for _ in range(B):
            vs = v[rng.integers(0, len(v), len(v))]; r = reconstruct(vs, U)
            pts.append((*screen_resid(r['mbar']), r['nS'], spectral_distance(r['A']), arms[U]))
    res_m0 = np.array([p[0] for p in pts]); res_m1 = np.array([p[1] for p in pts]); res_m2 = np.array([p[2] for p in pts])
    ucount = np.array([p[3] for p in pts]); sdv = np.array([p[4] for p in pts]); wrong = sdv > 0.05

    def auc(feat):
        p, q = feat[wrong], feat[~wrong]
        return sum((a > b) + 0.5 * (a == b) for a in p for b in q) / (len(p) * len(q)) if len(p) and len(q) else float('nan')
    arm_lab = np.array([p[5] for p in pts])
    out["cell17_roc"] = dict(auc_m0=float(auc(res_m0)), auc_m1=float(auc(res_m1)), auc_m2=float(auc(res_m2)),
                             auc_unique_count=float(auc(-ucount)), n_points=len(pts), n_wrong=int(wrong.sum()),
                             label_wrong_equals_arm_detuned=bool(np.array_equal(wrong, arm_lab == 'detuned')),
                             distinct_values=dict(m0=int(len(np.unique(np.round(res_m0, 12)))),
                                                  m1=int(len(np.unique(np.round(res_m1, 12)))),
                                                  m2=int(len(np.unique(np.round(res_m2, 12)))),
                                                  unique_count=int(len(np.unique(ucount))),
                                                  spectral_dist=int(len(np.unique(np.round(sdv, 12))))))
    # cell 19
    vc = data[U_TRUE]['clean']
    def swap_sites(s, a, b):
        for spin in (0, 1):
            pa = 2 * a + spin; pb = 2 * b + spin; ba = (s >> pa) & 1; bb = (s >> pb) & 1
            s = (s & ~((1 << pa) | (1 << pb))) | (ba << pb) | (bb << pa)
        return s
    batt = []
    for name, Ur in [('correct (U=4)', 4.0), ('wrong U=3.5', 3.5), ('wrong U=3.0', 3.0), ('wrong U=4.5', 4.5)]:
        rc = reconstruct(vc, Ur); r = screen_resid(rc['mbar'])
        batt.append(dict(name=name, U_rec=Ur, nS=rc['nS'], screen_pct=(100 * r).tolist(),
                         caught_by=[n for n, x in zip(['m0', 'm1', 'm2'], r) if x > 0.05]))
    sub = vc[:int(0.40 * len(vc))]
    good = reconstruct(sub, 4.0); bad = reconstruct(np.array([swap_sites(int(s), 0, 3) for s in sub]), 4.0)
    rg = screen_resid(good['mbar']); rb = screen_resid(bad['mbar'])
    out["cell19_battery"] = dict(full_coverage=batt, partial40_n_shots=int(len(sub)), partial40_nS_good=good['nS'],
                                 partial40_nS_swapped=bad['nS'], partial40_good_pct=(100 * rg).tolist(),
                                 partial40_swapped_pct=(100 * rb).tolist(),
                                 partial40_caught=[n for n, x in zip(['m0', 'm1', 'm2'], rb) if x > 0.05])
    # cell 20
    vd = data[DETUNE_U]['clean']
    A_lehmann = reconstruct(vd, DETUNE_U)['A']; A_haydock = haydock_A(vd, DETUNE_U)
    m2_leh = np.trapezoid(GRID ** 2 * A_lehmann, GRID); m2_hay = np.trapezoid(GRID ** 2 * A_haydock, GRID)
    out["cell20_transport"] = dict(
        m2_resid_lehmann_grid_pct=float(100 * abs(m2_leh - P4['m_ind'][2]) / P4['m_ind'][2]),
        m2_resid_haydock_grid_pct=float(100 * abs(m2_hay - P4['m_ind'][2]) / P4['m_ind'][2]),
        max_abs_A_lehmann_minus_haydock=float(np.abs(A_lehmann - A_haydock).max()),
        exact_m2_resid_pct=float(100 * abs(prep(DETUNE_U)['m_ind'][2] - P4['m_ind'][2]) / P4['m_ind'][2]),
        grid_m2_of_exact_A_over_exact_m2=dict(
            U4=float(np.trapezoid(GRID ** 2 * P4['A_exact'], GRID) / P4['m_ind'][2]),
            U3p5=float(np.trapezoid(GRID ** 2 * prep(DETUNE_U)['A_exact'], GRID) / prep(DETUNE_U)['m_ind'][2])),
        note=("full support: Lehmann and Haydock evaluate the same 300-dim resolvent, hence identical A; the "
              "grid value differs from the exact 19.44% because omega^2 A is integrated over the finite "
              "window [-1, 12] of a Lorentzian-broadened A, with omega measured from E0(U_rec)"))
    # cell 22
    def m2_resid_of(cfg, U): return screen_resid(reconstruct(cfg, U)['mbar'])[2]
    B = 300
    out["cell22_noise_survival"] = {}
    for U in arms:
        v = data[U]['clean']
        vals = np.array([m2_resid_of(v[rng.integers(0, len(v), len(v))], U) for _ in range(B)])
        valsT = np.array([m2_resid_of(v[rng.integers(0, len(v), int(0.97 * len(v)))], U) for _ in range(B)])
        out["cell22_noise_survival"][arms[U]] = dict(
            m2_resid_mean_pct=float(100 * vals.mean()), band16_84_pct=[float(100 * np.percentile(vals, 16)),
                                                                        float(100 * np.percentile(vals, 84))],
            readout_pert_mean_pct=float(100 * valsT.mean()), clears_5pct=bool(np.percentile(vals, 16) > 0.05))
    return out, data


pipeline = {}
pdata = {}
for name, f in MAPS.items():
    t0 = time.time()
    pipeline[name], pdata[name] = run_notebook_pipeline(f)
    log(f"pipeline [{name}] done in {time.time() - t0:.1f}s:", pipeline[name]["cell13_retention"],
        pipeline[name]["cell17_roc"]["auc_m2"])

# values printed by the executed notebook on 2026-08-21 (Colab outputs; transcript of that day)
reported_2026_08_21 = dict(
    m_ind=[0.5, 2.0836, 9.3344], cleaned_sector=dict(correct=13527, detuned=12969), unique=300,
    result1_screen_pct=[0.0, 0.0, 0.0], m2_hw=9.334, band=[9.334, 9.334],
    detuned_screen_pct=[0.0, 10.60, 19.44], auc=dict(m0=1.0, m1=1.0, m2=1.0, unique_count=0.5),
    transport_pct=16.62, noise_survival_detuned_pct=19.44)
ar = pipeline["as_run"]
reproduces = dict(
    cleaned_sector=(ar["cell13_retention"]["correct"]["cleaned_sector"] == 13527
                    and ar["cell13_retention"]["detuned"]["cleaned_sector"] == 12969),
    unique_300=(ar["cell13_retention"]["correct"]["unique"] == 300 and ar["cell13_retention"]["detuned"]["unique"] == 300),
    result1=(max(ar["cell15_result1"]["screen_pct"]) < 0.005 and abs(ar["cell15_result1"]["m2_hw"] - 9.334) < 5e-4),
    detuned_screen=(abs(ar["cell19_battery"]["full_coverage"][1]["screen_pct"][1] - 10.60) < 5e-3
                    and abs(ar["cell19_battery"]["full_coverage"][1]["screen_pct"][2] - 19.44) < 5e-3),
    auc=(ar["cell17_roc"]["auc_m0"] == 1.0 and ar["cell17_roc"]["auc_m2"] == 1.0 and ar["cell17_roc"]["auc_unique_count"] == 0.5),
    transport=abs(ar["cell20_transport"]["m2_resid_lehmann_grid_pct"] - 16.62) < 5e-3)
log("as-run reproduces the 2026-08-21 notebook outputs:", reproduces)

# =============================================================================================
# 4. What the screen tests
# =============================================================================================
m_exact = {U: prep(U)['m_ind'] for U in (U_TRUE, DETUNE_U)}
data_free = {}
for Ur in (U_TRUE, DETUNE_U):
    r = reconstruct(si, Ur)
    data_free[f"U_rec={Ur}"] = dict(nS=r['nS'], mbar=r['mbar'].tolist(), screen_pct=(100 * screen_resid(r['mbar'])).tolist(),
                                    spectral_dist=spectral_distance(r['A']),
                                    max_abs_mbar_minus_exact_moments_of_H_Urec=float(np.abs(r['mbar'] - m_exact[Ur]).max()))
cross = {}
max_dev_from_datafree = 0.0
for name in MAPS:
    for U in arms:
        v = pdata[name][U]['clean']
        for Ur in (U_TRUE, DETUNE_U):
            r = reconstruct(v, Ur); s = screen_resid(r['mbar'])
            df = np.array(data_free[f"U_rec={Ur}"]["screen_pct"])
            dev = float(np.abs(100 * s - df).max()); max_dev_from_datafree = max(max_dev_from_datafree, dev)
            cross[f"{name}|arm={arms[U]}|U_rec={Ur}"] = dict(nS=r['nS'], screen_pct=(100 * s).tolist(),
                                                             spectral_dist=spectral_distance(r['A']),
                                                             max_abs_diff_vs_data_free_pct=dev)
log("cross-assignment max |screen - data-free screen| (pct points):", max_dev_from_datafree)

# coverage scan: where the device data CAN enter (incomplete S); correct reading, U_rec = arm U
scan_rng = np.random.default_rng(SEED_SCAN)
NS_SCAN = (10, 20, 50, 100, 200, 400, 800, 1600, 3200)
R_SCAN = 200
coverage_scan = {}
for name in MAPS:
    for U in arms:
        v = pdata[name][U]['clean']
        rows = []
        for n in NS_SCAN:
            nS = []; res = []
            for _ in range(R_SCAN):
                sub = v[scan_rng.choice(len(v), size=n, replace=False)]
                r = reconstruct(sub, U); nS.append(r['nS']); res.append(screen_resid(r['mbar']))
            res = np.array(res); nS = np.array(nS)
            rows.append(dict(n_retained=n, mean_nS=float(nS.mean()), frac_full_support=float((nS == len(si)).mean()),
                             mean_screen_pct=(100 * res.mean(axis=0)).tolist(),
                             median_screen_pct=(100 * np.median(res, axis=0)).tolist()))
        # shots needed to complete the support in the shuffled order the notebook used
        seen = set(); n_full = None
        for t, x in enumerate(v):
            seen.add(int(x))
            if len(seen) == len(si): n_full = t + 1; break
        coverage_scan[f"{name}|arm={arms[U]}"] = dict(rows=rows, retained_shots_to_full_support_in_notebook_order=n_full)
log("coverage scan done")

# do the device counts themselves distinguish the schedules? (correct reading; not a moment test)
eps = 1e-6
discrim = {}
for U in arms:
    rows = []; llr_tot = 0.0; n_tot = 0
    for k in range(K):
        c = counts[U][k]
        xs = np.array([int(bs, 2) for bs in c]); ns = np.array(list(c.values()))
        keep = insec[xs]; xs = xs[keep]; ns = ns[keep]
        emp = np.zeros(dim); np.add.at(emp, xs, ns); emp /= emp.sum()
        q4 = P[U_TRUE][k] * insec; q4 /= q4.sum(); q35 = P[DETUNE_U][k] * insec; q35 /= q35.sum()
        tv_own = 0.5 * float(np.abs(emp - (q4 if U == U_TRUE else q35)).sum())
        tv_other = 0.5 * float(np.abs(emp - (q35 if U == U_TRUE else q4)).sum())
        llr = float((ns * (np.log(q4[xs] + eps) - np.log(q35[xs] + eps))).sum())
        llr_tot += llr; n_tot += int(ns.sum())
        rows.append(dict(k=k, retained=int(ns.sum()), tv_device_vs_own_schedule=tv_own,
                         tv_device_vs_other_schedule=tv_other,
                         tv_noiseless_U4_vs_U3p5=0.5 * float(np.abs(q4 - q35).sum()),
                         llr_U4_vs_U3p5_per_shot=llr / max(int(ns.sum()), 1)))
    discrim[arms[U]] = dict(per_circuit=rows, llr_U4_vs_U3p5_per_shot_pooled=llr_tot / n_tot,
                            n_circuits_closer_to_own_schedule=int(sum(r["tv_device_vs_own_schedule"] <
                                                                      r["tv_device_vs_other_schedule"] for r in rows
                                                                      if r["k"] >= 1)))
log("device-count discrimination:", {a: (v["llr_U4_vs_U3p5_per_shot_pooled"], v["n_circuits_closer_to_own_schedule"])
                                     for a, v in discrim.items()})

# =============================================================================================
# 5. Verdict and output
# =============================================================================================
cr = pipeline["correct"]
verdict = dict(
    what_the_screen_tests=(
        "Both arms have a complete support (|S| = 300/300) under both bit orders. The reconstruction "
        "diagonalises H(U_rec) on the whole sector, using the exact ED vector phi(U_rec) and E0(U_rec), so mbar "
        "equals the exact moments of H(U_rec) (difference 0.0, key data_free_control). The screen residual "
        "|mbar - m(U=4)|/m(U=4) is then a function of U_rec alone. Either arm's counts, read in either bit "
        "order, give the same screen as a data-free control (S = the full sector, no device data); the maximum "
        "difference is max_abs_diff_any_device_screen_vs_data_free_pct = 0.0. Reconstructing the correct arm "
        "with U_rec = 3.5 makes the screen fire, and reconstructing the detuned arm with U_rec = 4 makes it "
        "pass (key cross_assignment). The 'detuned arm caught' outcome therefore detects which Hamiltonian the "
        "analyst used, not the detuned device schedule. The ROC label (spectral distance > 0.05 from the U=4 "
        "ED spectrum) is also a function of U_rec, and it coincides with the arm label. All 200 bootstrap "
        "replicates of an arm give the same value (one distinct value per arm). AUC = 1 for m0 too, because "
        "its residuals are 0 and 1e-15 (rounding), not because m0 discriminates. The AUCs therefore carry no "
        "information about the device. Device data enter only through S, and only when S is incomplete (key "
        "coverage_scan): below a few thousand retained shots, the screen then responds to missing support."),
    bit_order=(
        "The notebook read bitstrings with int(bs[::-1],2). That map sends qubit p to bit 11-p (up<->down, "
        "site i<->5-i), so it maps the (N=7, Sz=+1) sector onto (N=7, Sz=-1), and every noiseless shot is "
        "rejected: the noiseless in-sector mass is 0 for all 14 circuits. Every shot retained by the as-run "
        "analysis (13527 and 12969) is therefore a device error with (N=7, Sz=-1) in the true reading, a set "
        "equality checked shot by shot. With the correct map int(bs,2), retention is 73364/245000 = 29.9% and "
        "73475/245000 = 30.0%. The support is again complete, so every screen and AUC output is unchanged."),
    device_data_carry_the_schedule=(
        "This is not a moment test. Read in the correct order, the retained counts favour their own schedule "
        "in pooled log-likelihood: +0.215 per shot (towards U=4) in the correct arm and -0.096 (towards U=3.5) "
        "in the detuned arm. Among circuits k = 1..6, the counts are closer in TV to their own schedule in "
        "5/6 (correct arm) and 4/6 (detuned arm). The device-vs-noiseless TV is 0.20-0.39 per circuit. The "
        "device output therefore carries a weak, noisy imprint of the detuning, and the moment screen as "
        "built does not use it."),
    paper=(
        "Not as a result. The run is not a hardware true positive or true negative of the moment screen, "
        "because the verdict does not depend on device data at complete support, and the numbers printed on "
        "2026-08-21 came from mis-ordered bitstrings. If the paper lists every device job, the run fits as one "
        "execution-record line: two ibm_fez Sampler jobs (ids), 2026-08-21, L=6, 7 circuits x 35 000 shots, "
        "69 QPU-s each. One sentence of what it shows would accompany it: the support was complete in both "
        "arms (correct-order retention 29.9% / 30.0%), so a support-based reconstruction with exact amplitudes "
        "returns the exact moments of whichever Hamiltonian is assumed, and the screen cannot see the device. "
        "Otherwise the run belongs only in the repository (data/heron_gatetp_2026-08-21/, this script, this "
        "JSON). No AUC, true-positive or noise-survival claim from this run should appear anywhere."))

out = dict(
    generator="src/gate_tp_analysis.py",
    date="2026-10-02",
    qpu_spent=0,
    provenance=dict(
        notebook="work/notebooks/GATE_TP_on_Heron.ipynb (saved without outputs; save cell 23 never executed)",
        raw_dir="data/heron_gatetp_2026-08-21/", raw_sha256=raw_files,
        credentials_check="no token, API key, CRN or instance field; job JSON carries user_id (an IBMid account id)",
        environment=dict(python=platform.python_version(), numpy=np.__version__,
                         qiskit=__import__("qiskit").__version__)),
    jobs=jobinfo, metrics=metrics,
    model=dict(L=L, t=thop, U_true=U_TRUE, U_detuned=DETUNE_U, sector="N=7, Sz=+1 (4 up, 3 down)",
               sector_dim=int(len(si)), K=K, n_trot=n_trot, dt_total=dt_total, shots_per_circuit=SHOTS,
               m_exact_U4=m_exact[U_TRUE].tolist(), m_exact_U3p5=m_exact[DETUNE_U].tolist(),
               E0_U4=float(prep(U_TRUE)['E0']), E0_U3p5=float(prep(DETUNE_U)['E0'])),
    decode=dict(method="numpy: base64 -> zlib -> .npy BitArray (shots x 2 bytes, big-endian), format(x,'012b')",
                runtime_crosscheck=runtime_crosscheck),
    bitorder=dict(noiseless=noiseless, submitted_isa_circuits=isa_check, device=device_bitorder,
                  mirror_note="int(bs[::-1],2) maps the (N=7,Sz=+1) sector onto (N=7,Sz=-1); verified set-equal"),
    notebook_pipeline=dict(as_run=pipeline["as_run"], correct=pipeline["correct"]),
    reported_2026_08_21=reported_2026_08_21, as_run_reproduces_reported=reproduces,
    screen_dependence=dict(data_free_control=data_free, cross_assignment=cross,
                           max_abs_diff_any_device_screen_vs_data_free_pct=max_dev_from_datafree,
                           coverage_scan=dict(seed=SEED_SCAN, replicas=R_SCAN, rows=coverage_scan,
                                              note="random subsets (without replacement) of the retained configs; "
                                                   "U_rec = the arm's own U, as in the notebook")),
    device_count_discrimination=dict(per_arm=discrim, eps_floor=eps,
                                     note="correct bit order; retained in-sector counts vs the normalised in-sector "
                                          "noiseless distributions of U=4 and U=3.5; llr > 0 favours U=4"),
    summary=dict(
        retention=dict(as_run=dict(correct=device_bitorder["correct"]["retained_as_run"],
                                   detuned=device_bitorder["detuned"]["retained_as_run"]),
                       correct_order=dict(correct=device_bitorder["correct"]["retained_correct"],
                                          detuned=device_bitorder["detuned"]["retained_correct"]),
                       raw_per_arm=device_bitorder["correct"]["raw"]),
        support=dict(as_run=dict(correct=ar["cell13_retention"]["correct"]["unique"],
                                 detuned=ar["cell13_retention"]["detuned"]["unique"]),
                     correct_order=dict(correct=cr["cell13_retention"]["correct"]["unique"],
                                        detuned=cr["cell13_retention"]["detuned"]["unique"])),
        screen_pct=dict(as_run=dict(correct=ar["cell15_result1"]["screen_pct"],
                                    detuned=ar["cell19_battery"]["full_coverage"][1]["screen_pct"]),
                        correct_order=dict(correct=cr["cell15_result1"]["screen_pct"],
                                           detuned=cr["cell19_battery"]["full_coverage"][1]["screen_pct"])),
        auc=dict(as_run={k: ar["cell17_roc"][k] for k in ("auc_m0", "auc_m1", "auc_m2", "auc_unique_count")},
                 correct_order={k: cr["cell17_roc"][k] for k in ("auc_m0", "auc_m1", "auc_m2", "auc_unique_count")})),
    verdict=verdict,
    runtime_s=time.time() - T_START)


def _clean(o):
    if isinstance(o, dict): return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [_clean(v) for v in o]
    if isinstance(o, (np.floating,)): return float(o)
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, np.bool_): return bool(o)
    return o


with open(OUT, "w", encoding="utf-8") as fh:
    json.dump(_clean(out), fh, indent=1, ensure_ascii=False)
log("wrote", OUT)
log("summary:", json.dumps(_clean(out["summary"])))
