# -*- coding: utf-8 -*-
r"""
Reference distributions for the device coverage residual Delta_0 (plan R1; blocking items B2, B3, B4).

What the ibm_fez L=8 run supplies is ONLY the post-selected support S of 7 Krylov--Trotter circuits
(hardware_matched_job_L8.py). Delta_0 = |m0 - m0_bar(S)| is the exact weight of phi = c^dag_{0,up}|0>
missed by S. This script asks what Delta_0 (and |S|, Delta_1, rel_L1) look like for reference samplers
evaluated with the SAME post-selection and the SAME falsifier arithmetic, so the device values can be
read against distributions instead of against a single Aer seed or an analytic sigma.

Reference samplers (all use the exact noiseless output distributions P_k(x) = |<x|U_k|det>|^2 of the
7 circuits, k = 0..6, from qiskit Statevector; nothing is fitted to the device):

  (i)   NOISELESS, RAW-MATCHED: per circuit a multinomial draw of the device's raw shots n_k
        (50k/30k/16k/4k), post-selected (noiseless retention is 100% up to 1e-14).
  (ii)  NOISELESS, RETAINED-MATCHED: per circuit a multinomial draw of the device's RETAINED
        (post-selected) count for that circuit, i.e. the same usable statistics as the device.
  (iii) UNIFORM-NOISE MIXTURE, RAW-MATCHED: per circuit (1-p) P_k + p * uniform(2^16), n_k raw shots,
        then the same post-selection, for p in {0, .1, .25, .5, .75, 1}; plus (extra) the p whose
        expected retention equals the device's pooled retention at that budget.
  (iv)  UNIFORM IN-SECTOR SAMPLER at the device's pooled retained count N: E[Delta_0] = m0 (1-1/3920)^N
        exactly; Monte Carlo check.

Exact sampling identity used for speed: drawing n shots from a distribution q over 2^16 outcomes and
post-selecting to the sector is identical in law to K ~ Binomial(n, q(sector)) followed by a multinomial
draw of K over the sector with probabilities q|sector / q(sector). So all draws are over the 3920 sector
configurations only.

Bit order: the Statevector index i and the qiskit counts bitstring format(i, '016b') (rightmost char =
qubit 0) map to the same Fock integer through the job's bitstring_to_config = int(bs, 2) (the CORRECT
map; the reversed companion map int(bs[::-1], 2) is the documented bug). Checked below three ways:
(a) the job's postselect() over all 2^16 bitstrings returns exactly the sector si; (b) the stored Aer
dry-run counts of the same circuits (heron_counts_matched_L8_DRYRUN.json) put 100% of their shots on
configurations with P_k > 0 and reproduce the stored dry-run Delta_0; (c) under the reversed map the
noiseless in-sector mass is reported (the bug rejects every noiseless shot).

rel_L1 (the B4 comparator: relative L1 error of the reconstructed A(omega), eta = 0.15, 600-pt grid of
the job) is computed by Lanczos with full reorthogonalisation on H restricted to S (Gauss quadrature of
the Lorentzian-broadened spectral function); the result is validated in the JSON against the stored
exact-eigh values of the four device runs and of the Aer dry run.

Comparator of the "sd / analytic sigma" ratio: the analytic sigma is sd_delta0 of the SUPERSEDED model
data/_superseded/analytic_coverage_leak.json (draws from |phi|^2 at the device's pooled kept N, Poisson
independence). Both Monte Carlo comparators are reported: sd(i) raw-matched and sd(ii) retained-matched.

Output: data/2026-09-27_delta0_reference_mc.json  (new file; no existing data file is touched).
Figure data: src/export_device_dat.py reads that JSON and writes paper/figs/device_*.dat.

USAGE (0 QPU):  python src/delta0_reference_mc.py        (about 10 min on 2 BLAS threads)
Env overrides: SEED, R_RAW, R_KEPT, R_MIX, R_UNI, R_REL_RAW, R_REL_KEPT, R_REL_UNI, LANCZOS_M.
"""
import os, sys, json, time, zlib, hashlib, platform, importlib.util
import numpy as np
from scipy.linalg import eigh_tridiagonal

T_START = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "data"))
OUT = os.environ.get("OUT_JSON", os.path.join(RES, "2026-09-27_delta0_reference_mc.json"))
ANALYTIC = os.path.join(RES, "_superseded", "analytic_coverage_leak.json")
DRYRUN = os.path.join(RES, "heron_counts_matched_L8_DRYRUN.json")
DEVICE_FILES = {50000: "heron_counts_matched_L8_da983qse74ec73ajfr2g.json",
                30000: "heron_counts_matched_L8_da98k16rbfbs73chckvg.json",
                16000: "heron_counts_matched_L8_da98icsjbipc73fej78g.json",
                4000: "heron_counts_matched_L8_da98cksjbipc73fej1tg.json"}
BUDGETS = (50000, 30000, 16000, 4000)
P_MIX = (0.0, 0.1, 0.25, 0.5, 0.75, 1.0)

SEED = int(os.environ.get("SEED", 20260927))
R_RAW = int(os.environ.get("R_RAW", 2000))
R_KEPT = int(os.environ.get("R_KEPT", 2000))
R_MIX = int(os.environ.get("R_MIX", 1000))
R_UNI = int(os.environ.get("R_UNI", 2000))
R_REL_RAW = int(os.environ.get("R_REL_RAW", 80))
R_REL_KEPT = int(os.environ.get("R_REL_KEPT", 40))
R_REL_UNI = int(os.environ.get("R_REL_UNI", 10))
LANCZOS_M = int(os.environ.get("LANCZOS_M", 300))

_trapz = getattr(np, "trapezoid", None) or np.trapz
log = lambda *a: print(f"[{time.time() - T_START:7.1f}s]", *a, flush=True)
timing = {}


def rng_for(tag):
    """Independent, reproducible stream per block (does not depend on the replica counts of other blocks)."""
    return np.random.default_rng(np.random.SeedSequence([SEED, zlib.crc32(tag.encode())]))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


# ---------------------------------------------------------------------------------------------
# 0. Job machinery (H, exact phi, circuits, post-selection) -- imported, never re-derived
# ---------------------------------------------------------------------------------------------
log("importing hardware_matched_job_L8.py (builds H, diagonalises both sectors; no QPU) ...")
_spec = importlib.util.spec_from_file_location("hwjob_L8", os.path.join(HERE, "hardware_matched_job_L8.py"))
mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mod)
timing["import_job_module_s"] = time.time() - T_START

si = np.asarray(mod.si)
NSEC = len(si)
NQ = 2 * mod.L
DIM = 1 << NQ
pos = mod.pos
phi = np.asarray(mod.phi)[si]
assert np.abs(phi.imag).max() < 1e-12
phi = phi.real.copy()
w = phi ** 2
m0 = float(mod.m0_op)
m1 = float(mod.m1_op)
E0 = float(mod.E0)
H_add = mod.H[si][:, si].tocsr()
grid = np.asarray(mod.grid)
A_exact = np.asarray(mod.A_exact)
eta = float(mod.eta)
A_exact_norm = float(_trapz(A_exact, grid))

# ---------------------------------------------------------------------------------------------
# 1. Exact noiseless output distributions of the 7 circuits
# ---------------------------------------------------------------------------------------------
t0 = time.time()
from qiskit.quantum_info import Statevector
P = np.array([np.abs(Statevector(qc.remove_final_measurements(inplace=False)).data) ** 2
              for qc in mod.circuits])                      # K x 2^16, index i <-> bitstring format(i,'016b')
K = P.shape[0]
Q0 = P[:, si]                                               # in-sector probabilities per circuit
mass0 = Q0.sum(axis=1)
timing["statevectors_s"] = time.time() - t0
_ws = np.sort(w)[::-1]
union_inf = (Q0 > 1e-12).any(axis=0)
phi_weight_stats = dict(
    n_configs=NSEC, n_nonzero_weight=int((w > 1e-14).sum()), max_weight=float(_ws[0]),
    mean_weight=float(w.mean()), median_weight=float(np.median(w)),
    share_of_m0_in_top_1pct_configs=float(_ws[:NSEC // 100].sum() / m0),
    share_of_m0_in_top_10pct_configs=float(_ws[:NSEC // 10].sum() / m0),
    note="w_x = |phi_x|^2 is heavy-tailed, so Monte Carlo means of rare-miss Delta_0 are dominated by rare "
         "misses of heavy configurations; the analytic E[Delta_0] values in this file are exact")
noiseless_infinite_shot = dict(
    union_support_size=int(union_inf.sum()), delta0_limit=float(m0 - w[union_inf].sum()),
    support_size_per_circuit=[int((Q0[k] > 1e-12).sum()) for k in range(Q0.shape[0])],
    note="support of the exact output distributions (P_k > 1e-12): the shots -> infinity limit of (i)")
log(f"statevectors: K={K}, in-sector mass per circuit = {np.round(mass0, 14).tolist()}")

# ---------------------------------------------------------------------------------------------
# 2. Bit-order verification (the job's CORRECT int(bs,2) map)
# ---------------------------------------------------------------------------------------------
t0 = time.time()
all_bs = [format(i, f"0{NQ}b") for i in range(DIM)]
seen_all, _, _ = mod.postselect({bs: 1 for bs in all_bs})
check_a = (seen_all == set(int(x) for x in si))
check_a_identity = all(mod.bitstring_to_config(all_bs[i]) == i for i in range(DIM))
rev_idx = np.array([int(bs[::-1], 2) for bs in all_bs])     # reversed companion map
insec = np.zeros(DIM, bool); insec[si] = True
rev_mass = [float(P[k][insec[rev_idx]].sum()) for k in range(K)]
dry = json.load(open(DRYRUN, encoding="utf-8"))
dry_on_support = []; dry_tv = []
pooled = {}
for k, cnt in enumerate(dry["per_time_counts"]):
    n = sum(cnt.values()); on = 0
    emp = np.zeros(DIM)
    for bs, c in cnt.items():
        x = mod.bitstring_to_config(bs); emp[x] += c
        if P[k][x] > 1e-12:
            on += c
        pooled[bs] = pooled.get(bs, 0) + c
    dry_on_support.append(on / n)
    dry_tv.append(0.5 * float(np.abs(emp / n - P[k]).sum()))
_rng_tv = rng_for("tv_check")
dry_tv_expected_one_draw = [0.5 * float(np.abs(_rng_tv.multinomial(sum(c_.values()), P[k] / P[k].sum())
                                               / sum(c_.values()) - P[k]).sum())
                            for k, c_ in enumerate(dry["per_time_counts"])]
dry_seen, dry_kept, dry_total = mod.postselect(pooled)
dry_d0 = abs(m0 - float(sum(w[pos[x]] for x in dry_seen)))
bitorder = dict(
    map_used="hardware_matched_job_L8.bitstring_to_config = int(bs, 2) (rightmost char = qubit 0)",
    postselect_over_all_65536_bitstrings_equals_sector=bool(check_a),
    statevector_index_equals_config_int=bool(check_a_identity),
    noiseless_insector_mass_correct_map=mass0.tolist(),
    noiseless_insector_mass_reversed_map=rev_mass,
    aer_dryrun_fraction_of_shots_on_Pk_support=dry_on_support,
    aer_dryrun_total_variation_vs_Pk=dry_tv,
    total_variation_of_one_simulated_draw_from_Pk_same_n=dry_tv_expected_one_draw,
    noiseless_support_size_per_circuit=[int((P[k] > 1e-12).sum()) for k in range(K)],
    aer_dryrun_delta0_recomputed=dry_d0,
    aer_dryrun_delta0_stored=float(dry["falsifier"]["delta0"]),
)
assert check_a and check_a_identity and min(dry_on_support) == 1.0
assert abs(dry_d0 - dry["falsifier"]["delta0"]) < 1e-12
timing["bitorder_checks_s"] = time.time() - t0
log(f"bit order OK: reversed-map noiseless in-sector mass = {max(rev_mass):.2e}; dry-run Delta0 {dry_d0:.6f}")


# ---------------------------------------------------------------------------------------------
# 3. Falsifier arithmetic on a sector mask (identical to the job's, vectorised)
# ---------------------------------------------------------------------------------------------
def residuals(seen):
    v = np.where(seen, phi, 0.0)
    m0b = float(v @ v)
    m1b = float(v @ (H_add @ v)) - E0 * m0b
    return abs(m0 - m0b), int(seen.sum()), abs(m1 - m1b)


def rel_l1(seen, m=LANCZOS_M):
    """rel L1 of A_S vs A_exact on the job grid; Lanczos (full DGKS reorth.) Gauss quadrature on H_S."""
    Sidx = np.flatnonzero(seen)
    HS = H_add[Sidx][:, Sidx].tocsr()
    v0 = phi[Sidx]; nrm = float(np.linalg.norm(v0))
    n = len(Sidx); m = min(m, n)
    V = np.zeros((m, n)); a = np.zeros(m); b = np.zeros(m)
    V[0] = v0 / nrm; kk = m
    for j in range(m):
        x = HS @ V[j]; a[j] = V[j] @ x
        nx = np.linalg.norm(x)
        x -= V[:j + 1].T @ (V[:j + 1] @ x)
        if np.linalg.norm(x) < 0.7071 * nx:                   # DGKS: second pass only when cancellation
            x -= V[:j + 1].T @ (V[:j + 1] @ x)
        if j + 1 < m:
            b[j] = np.linalg.norm(x)
            if b[j] < 1e-12:
                kk = j + 1; break
            V[j + 1] = x / b[j]
    th, Y = eigh_tridiagonal(a[:kk], b[:kk - 1])
    wt = nrm ** 2 * Y[0] ** 2
    A = ((wt[None, :] * (eta / np.pi)) / ((grid[:, None] - (th - E0)[None, :]) ** 2 + eta ** 2)).sum(1)
    return float(_trapz(np.abs(A - A_exact), grid) / A_exact_norm)


def draw_support(rng, nk, Q, mass):
    """One replica: per circuit n_k shots from q_k over 2^16, post-selected to the sector (exact in law)."""
    seen = np.zeros(NSEC, bool); kept = 0
    for k in range(len(nk)):
        mk = min(float(mass[k]), 1.0)
        Kk = int(nk[k]) if mk >= 1.0 else int(rng.binomial(int(nk[k]), mk))
        if Kk:
            seen |= rng.multinomial(Kk, Q[k] / Q[k].sum()) > 0
        kept += Kk
    return seen, kept


def summarize(x):
    x = np.asarray(x, float)
    if len(x) == 0:
        return None
    q = np.percentile(x, [2.5, 16, 50, 84, 97.5])
    return dict(n=int(len(x)), mean=float(x.mean()), sd=float(x.std(ddof=1)) if len(x) > 1 else 0.0,
                min=float(x.min()), max=float(x.max()), p2_5=float(q[0]), p16=float(q[1]), p50=float(q[2]),
                p84=float(q[3]), p97_5=float(q[4]))


def rank(x, v):
    x = np.asarray(x, float)
    if len(x) == 0:
        return None
    return dict(frac_le=float((x <= v).mean()), frac_ge=float((x >= v).mean()), n=int(len(x)))


def exact_E(nk, Q):
    """Exact E[Delta_0] and E[|S|] for independent per-circuit multinomials (miss prob prod_k (1-q_kx)^n_k)."""
    with np.errstate(divide="ignore"):
        lm = np.log1p(-np.minimum(Q, 1 - 1e-300))
    miss = np.exp((np.asarray(nk, float)[:, None] * lm).sum(0))
    return float((w * miss).sum()), float((1 - miss).sum())


def run_block(tag, nk, Q, mass, R, R_rel):
    rng = rng_for(tag)
    d0 = np.empty(R); S = np.empty(R); d1 = np.empty(R); kept = np.empty(R); rels = []
    for r in range(R):
        seen, kp = draw_support(rng, nk, Q, mass)
        d0[r], S[r], d1[r] = residuals(seen); kept[r] = kp
        if r < R_rel:
            rels.append(rel_l1(seen))
    return d0, S, d1, kept, np.array(rels)


# ---------------------------------------------------------------------------------------------
# 4. Device rows, recomputed from the retained raw counts (per circuit raw and kept shots)
# ---------------------------------------------------------------------------------------------
t0 = time.time()
device = {}; device_seen = {}
for NS in BUDGETS:
    fp = os.path.join(RES, DEVICE_FILES[NS])
    d = json.load(open(fp, encoding="utf-8"))
    assert not d.get("_dry_run") and d["shots"] == NS
    raw_k, kept_k = [], []; pooled = {}
    for cnt in d["per_time_counts"]:
        s_k, kp_k, tot_k = mod.postselect(cnt)
        raw_k.append(int(tot_k)); kept_k.append(int(kp_k))
        for bs, c in cnt.items():
            pooled[bs] = pooled.get(bs, 0) + c
    seen_set, kept, total = mod.postselect(pooled)
    seen = np.zeros(NSEC, bool); seen[[pos[x] for x in seen_set]] = True
    device_seen[NS] = seen
    per_k = []
    masks_k = []
    for cnt in d["per_time_counts"]:
        s_k = mod.postselect(cnt)[0]
        mk_ = np.zeros(NSEC, bool); mk_[[pos[x] for x in s_k]] = True; masks_k.append(mk_)
    for k in range(len(masks_k)):
        others = np.zeros(NSEC, bool)
        for j in range(len(masks_k)):
            if j != k:
                others |= masks_k[j]
        per_k.append(dict(k=k, S_k=int(masks_k[k].sum()), delta0_this_circuit_alone=residuals(masks_k[k])[0],
                          delta0_union_without_this_circuit=residuals(others)[0],
                          unique_configs_only_from_this_circuit=int((masks_k[k] & ~others).sum())))
    dd0, dS, dd1 = residuals(seen)
    drel = rel_l1(seen, m=400)
    f = d["falsifier"]
    assert dS == f["S_size"] and abs(dd0 - f["delta0"]) < 1e-12 and kept == f["kept_shots"]
    device[NS] = dict(file=DEVICE_FILES[NS], job_id=d["job_id"], backend=d["backend"],
                      raw_per_circuit=raw_k, kept_per_circuit=kept_k, kept_total=int(kept), raw_total=int(total),
                      retention=kept / total, S=dS, coverage=dS / NSEC, delta0=dd0, delta1=dd1,
                      rel_L1_lanczos400=drel, rel_L1_stored_exact_eigh=float(f["rel_L1"]),
                      delta0_stored=float(f["delta0"]), delta1_stored=float(f["delta1"]),
                      per_circuit_support=per_k)
    log(f"device {NS}: |S|={dS} kept={kept}/{total} Delta0={dd0:.6f} Delta1={dd1:.6f} "
        f"rel_L1={drel:.6f} (stored {f['rel_L1']:.6f})")
dry_seen_mask = np.zeros(NSEC, bool); dry_seen_mask[[pos[x] for x in dry_seen]] = True
lanczos_validation = dict(
    method=f"Lanczos, full DGKS reorthogonalisation, m={LANCZOS_M} (device rows m=400), vs stored exact eigh",
    dryrun_50k_seed1=dict(stored=float(dry["falsifier"]["rel_L1"]),
                          lanczos_m=rel_l1(dry_seen_mask), lanczos_m400=rel_l1(dry_seen_mask, m=400)),
    device_m400_minus_stored={str(NS): device[NS]["rel_L1_lanczos400"] - device[NS]["rel_L1_stored_exact_eigh"]
                              for NS in BUDGETS},
    device_m_minus_stored={},
)
for NS in BUDGETS:
    lanczos_validation["device_m_minus_stored"][str(NS)] = (rel_l1(device_seen[NS])
                                                            - device[NS]["rel_L1_stored_exact_eigh"])
timing["device_rows_s"] = time.time() - t0

analytic = {int(r["shots"]): r for r in json.load(open(ANALYTIC, encoding="utf-8"))["rows"]}
RMC = json.load(open(os.path.join(RES, "2026-09-05_retention_matched_control.json"), encoding="utf-8"))
AER1 = json.load(open(os.path.join(RES, "aer_noiseless_baseline.json"), encoding="utf-8"))

# ---------------------------------------------------------------------------------------------
# 5. Reference distributions per budget
# ---------------------------------------------------------------------------------------------
UNI_SECTOR_FRAC = NSEC / DIM
budgets_out = {}
for NS in BUDGETS:
    tb = time.time()
    D = device[NS]
    raw_k = np.array(D["raw_per_circuit"]); kept_k = np.array(D["kept_per_circuit"])
    row = dict(device=dict(delta0=D["delta0"], S=D["S"], delta1=D["delta1"], rel_L1=D["rel_L1_stored_exact_eigh"],
                           kept_total=D["kept_total"], retention=D["retention"]))

    # (i) noiseless raw-matched
    t0 = time.time()
    d0, S, d1, kp, rels = run_block(f"raw{NS}", raw_k, Q0, mass0, R_RAW, R_REL_RAW)
    E_raw, ES_raw = exact_E(raw_k, Q0)
    row["noiseless_raw_matched"] = dict(
        shots_per_circuit=raw_k.tolist(), replicas=R_RAW,
        delta0=summarize(d0), S=summarize(S), delta1=summarize(d1), kept=summarize(kp),
        rel_L1=summarize(rels), rel_L1_replicas=len(rels),
        exact_E_delta0=E_raw, exact_E_S=ES_raw,
        device_rank_delta0=rank(d0, D["delta0"]), device_rank_S=rank(S, D["S"]),
        device_rank_delta1=rank(d1, D["delta1"]), device_rank_rel_L1=rank(rels, D["rel_L1_stored_exact_eigh"]),
        runtime_s=time.time() - t0)
    if NS == 50000:
        row["noiseless_raw_matched"]["aer_dryrun_seed1_rank_delta0"] = rank(d0, dry_d0)
        row["noiseless_raw_matched"]["aer_dryrun_seed1_rank_rel_L1"] = rank(rels, float(dry["falsifier"]["rel_L1"]))
    log(f"{NS} raw: Delta0 {d0.mean():.4f}+-{d0.std(ddof=1):.4f} (exact E {E_raw:.4f}); "
        f"device frac_le={rank(d0, D['delta0'])['frac_le']:.4f}; rel_L1 {rels.mean():.3f}+-{rels.std(ddof=1):.3f}")

    # (ii) noiseless retained-matched (per-circuit device kept counts)
    t0 = time.time()
    d0k, Sk, d1k, kpk, relk = run_block(f"kept{NS}", kept_k, Q0, mass0, R_KEPT, R_REL_KEPT)
    E_kept, ES_kept = exact_E(kept_k, Q0)
    row["noiseless_retained_matched"] = dict(
        shots_per_circuit=kept_k.tolist(), replicas=R_KEPT,
        delta0=summarize(d0k), S=summarize(Sk), delta1=summarize(d1k),
        rel_L1=summarize(relk), rel_L1_replicas=len(relk),
        exact_E_delta0=E_kept, exact_E_S=ES_kept,
        device_rank_delta0=rank(d0k, D["delta0"]), device_rank_S=rank(Sk, D["S"]),
        device_rank_delta1=rank(d1k, D["delta1"]), device_rank_rel_L1=rank(relk, D["rel_L1_stored_exact_eigh"]),
        runtime_s=time.time() - t0)
    log(f"{NS} kept: Delta0 {d0k.mean():.4f}+-{d0k.std(ddof=1):.4f} (exact E {E_kept:.4f}); "
        f"device frac_le={rank(d0k, D['delta0'])['frac_le']:.4f}")

    # (ii-b) noiseless retained-matched, EQUAL split round(kept_total/7) per circuit: the allocation used by
    #        retention_matched_control.py (single Aer seed 1) and by the analytic kept-count model
    t0 = time.time()
    ns_eq = int(round(D["kept_total"] / K))
    eq_k = np.full(K, ns_eq)
    d0e, Se, d1e, kpe, rele = run_block(f"keptEq{NS}", eq_k, Q0, mass0, R_KEPT, R_REL_KEPT // 2)
    E_eq, ES_eq = exact_E(eq_k, Q0)
    ctrl = {int(r_["raw_per_circuit"]): r_ for r_ in RMC["rows"]}[NS]
    row["noiseless_retained_matched_equal_split"] = dict(
        shots_per_circuit=eq_k.tolist(), replicas=R_KEPT,
        note="same total usable shots as the device, but spread equally over the 7 circuits; the device "
             "instead keeps ~92% of circuit k=0 (a single determinant when noiseless) and ~12% of k=1..6",
        delta0=summarize(d0e), S=summarize(Se), delta1=summarize(d1e),
        rel_L1=summarize(rele), rel_L1_replicas=len(rele),
        exact_E_delta0=E_eq, exact_E_S=ES_eq,
        device_rank_delta0=rank(d0e, D["delta0"]), device_rank_S=rank(Se, D["S"]),
        retention_matched_control_seed1_delta0=float(ctrl["ctrl_d0"]),
        retention_matched_control_seed1_rank=rank(d0e, float(ctrl["ctrl_d0"])),
        runtime_s=time.time() - t0)
    log(f"{NS} kept-equal-split: Delta0 {d0e.mean():.4f}+-{d0e.std(ddof=1):.4f} (exact E {E_eq:.4f}); "
        f"device frac_le={rank(d0e, D['delta0'])['frac_le']:.4f}")
    aer1 = {int(r_["shots"]): r_ for r_ in AER1["rows"]}[NS]
    row["noiseless_raw_matched"]["published_aer_seed1_delta0"] = float(aer1["d0_noiseless"])
    row["noiseless_raw_matched"]["published_aer_seed1_rank"] = rank(d0, float(aer1["d0_noiseless"]))

    # sd / analytic sigma (superseded |phi|^2 model at kept N)
    A_ = analytic[NS]
    row["ratio_reference_mean_over_device_delta0"] = dict(
        noiseless_raw_matched=float(d0.mean() / D["delta0"]),
        noiseless_retained_matched=float(d0k.mean() / D["delta0"]),
        noiseless_retained_matched_equal_split=float(d0e.mean() / D["delta0"]),
        retention_matched_control_seed1_single_draw=float(ctrl["ctrl_d0"] / D["delta0"]))
    row["sd_over_analytic_sigma"] = dict(
        analytic_source="data/_superseded/analytic_coverage_leak.json rows[shots].sd_delta0 "
                        "(draws from |phi|^2 at the device's pooled kept N, Poisson independence)",
        analytic_sigma=float(A_["sd_delta0"]), analytic_E=float(A_["E_delta0"]), analytic_kept=int(A_["kept"]),
        ratio_sd_rawmatched_MC=float(d0.std(ddof=1) / A_["sd_delta0"]),
        ratio_sd_retainedmatched_MC=float(d0k.std(ddof=1) / A_["sd_delta0"]),
        ratio_sd_retainedmatched_equal_split_MC=float(d0e.std(ddof=1) / A_["sd_delta0"]))

    # (iii) uniform-noise mixture, raw-matched
    t0 = time.time()
    p_ret = float(np.clip((1 - D["retention"]) / (1 - UNI_SECTOR_FRAC), 0, 1))
    mix = []
    for p in list(P_MIX) + [p_ret]:
        Qp = (1 - p) * Q0 + p / DIM
        massp = Qp.sum(axis=1)
        d0m, Sm, d1m, kpm, _ = run_block(f"mix{NS}_{p:.6f}", raw_k, Qp, massp, R_MIX, 0)
        Em, ESm = exact_E(raw_k, Qp)
        mix.append(dict(p=p, retention_matched_extra=(p == p_ret and p not in P_MIX), replicas=R_MIX,
                        expected_retention=float((raw_k * np.minimum(massp, 1)).sum() / raw_k.sum()),
                        delta0=summarize(d0m), S=summarize(Sm), delta1=summarize(d1m),
                        exact_E_delta0=Em, exact_E_S=ESm,
                        device_rank_delta0=rank(d0m, D["delta0"]), device_rank_S=rank(Sm, D["S"])))
        log(f"{NS} mix p={p:.3f}: Delta0 {d0m.mean():.5f}+-{d0m.std(ddof=1):.5f} (exact E {Em:.5f}) "
            f"ret={mix[-1]['expected_retention']:.4f}")
    e_seq = [m_["exact_E_delta0"] for m_ in mix[:len(P_MIX)]]
    row["uniform_noise_mixture_raw_matched"] = dict(
        definition="per circuit (1-p) P_k + p/2^16, device raw shots per circuit, same post-selection",
        p_retention_matched=p_ret, rows=mix,
        exact_E_monotone_decreasing_in_p=bool(all(np.diff(e_seq) < 0)),
        mc_mean_monotone_decreasing_in_p=bool(all(np.diff([m_["delta0"]["mean"] for m_ in mix[:len(P_MIX)]]) < 0)),
        runtime_s=time.time() - t0)

    # (iv) uniform in-sector sampler at the device's pooled kept count
    t0 = time.time()
    N = D["kept_total"]
    rng = rng_for(f"uni{NS}")
    d0u = np.empty(R_UNI); Su = np.empty(R_UNI); d1u = np.empty(R_UNI); relu = []
    pu = np.full(NSEC, 1.0 / NSEC)
    for r in range(R_UNI):
        seen = rng.multinomial(N, pu) > 0
        d0u[r], Su[r], d1u[r] = residuals(seen)
        if r < R_REL_UNI:
            relu.append(rel_l1(seen))
    E_uni = m0 * (1 - 1 / NSEC) ** N
    row["uniform_in_sector_at_device_kept"] = dict(
        N_kept=int(N), replicas=R_UNI, analytic_E_delta0=float(E_uni),
        analytic_E_S=float(NSEC * (1 - (1 - 1 / NSEC) ** N)),
        delta0=summarize(d0u), S=summarize(Su), delta1=summarize(d1u),
        rel_L1=summarize(relu), rel_L1_replicas=len(relu),
        mc_check_missed_configs=dict(
            mc_mean=float((NSEC - Su).mean()), analytic=float(NSEC * (1 - 1 / NSEC) ** N),
            z=float(((NSEC - Su).mean() - NSEC * (1 - 1 / NSEC) ** N)
                    / max((NSEC - Su).std(ddof=1) / np.sqrt(R_UNI), 1e-300)) if (NSEC - Su).std(ddof=1) > 0 else None,
            note="the MC check uses the missed-configuration count (light-tailed); the MC mean of Delta_0 itself is "
                 "dominated by rare misses of heavy configurations (see phi_weight_stats), the analytic E is exact"),
        device_minus_uniform_E=float(D["delta0"] - E_uni),
        runtime_s=time.time() - t0)
    log(f"{NS} uniform in-sector N={N}: E={E_uni:.3e}, MC {d0u.mean():.3e}+-{d0u.std(ddof=1):.3e}; "
        f"rel_L1 {np.mean(relu):.3f}")
    row["runtime_s"] = time.time() - tb
    budgets_out[str(NS)] = row

# ---------------------------------------------------------------------------------------------
# 6. Write
# ---------------------------------------------------------------------------------------------
import scipy, qiskit
runtime = time.time() - T_START
out = dict(
    provenance=dict(
        script="src/delta0_reference_mc.py", date="2026-09-27", seed=SEED,
        seed_scheme="np.random.default_rng(SeedSequence([seed, crc32(block_tag)])) per block",
        runtime_s=runtime, timing=timing,
        replicas=dict(R_RAW=R_RAW, R_KEPT=R_KEPT, R_MIX=R_MIX, R_UNI=R_UNI, R_REL_RAW=R_REL_RAW,
                      R_REL_KEPT=R_REL_KEPT, R_REL_UNI=R_REL_UNI, LANCZOS_M=LANCZOS_M),
        inputs={os.path.relpath(p, os.path.join(RES, "..")).replace("\\", "/"): sha256(p)
                for p in [os.path.join(RES, f) for f in DEVICE_FILES.values()] + [DRYRUN, ANALYTIC,
                          os.path.join(RES, "2026-09-05_retention_matched_control.json"),
                          os.path.join(RES, "aer_noiseless_baseline.json"),
                          os.path.join(HERE, "hardware_matched_job_L8.py")]},
        versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                      qiskit=qiskit.__version__),
        qpu=False),
    system=dict(L=int(mod.L), U=float(mod.U), eta=eta, K=K, n_trot=int(mod.n_trot), dt_total=float(mod.dt_total),
                sector_dim=NSEC, fock_dim=DIM, m0=m0, m1=m1, E0=E0, uniform_sector_fraction_of_2to16=UNI_SECTOR_FRAC),
    definitions=dict(
        delta0="|m0 - ||P_S phi||^2|, exact missed weight of phi by the post-selected support S",
        delta1="|m1 - (<P_S phi|H|P_S phi> - E0 ||P_S phi||^2)|",
        rel_L1="int |A_S - A_exact| / int A_exact on the job's 600-pt grid, eta=0.15 (the B4 metric)",
        frac_le="fraction of reference replicas with value <= the device value (empirical percentile rank)",
        noiseless_raw_matched="(i) multinomial of the device's raw shots per circuit from the exact P_k",
        noiseless_retained_matched="(ii) multinomial of the device's retained shots per circuit from the exact P_k",
        noiseless_retained_matched_equal_split="(ii-b) as (ii) with the device's total retained shots split "
                                                "equally over the 7 circuits (retention_matched_control.py allocation)"),
    phi_weight_stats=phi_weight_stats,
    noiseless_infinite_shot=noiseless_infinite_shot,
    bit_order_checks=bitorder,
    lanczos_validation=lanczos_validation,
    device=device,
    budgets=budgets_out,
)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1)
log(f"wrote {OUT}  (runtime {runtime:.1f} s)")
