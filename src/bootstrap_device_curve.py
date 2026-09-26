r"""
Bootstrap confidence intervals on the device-side falsifier (Delta_0, Delta_1)
from the RETAINED raw ibm_fez counts, for every shot level of the L=8 curve.

Purpose (referee-facing): SEPARATE the two contributions to Delta_0 that the
one-number summary conflated:
  (a) the DETERMINISTIC residual coverage leak  = m0_op * (missed weight fraction),
      set by which determinants the support S happens to miss;
  (b) the STOCHASTIC shot-noise dispersion       = spread of Delta_0 when the same
      finite shot budget is re-drawn (nonparametric bootstrap of the raw counts).

Method: for each of the K=7 Krylov-time circuits, resample its raw shots WITH
REPLACEMENT (multinomial over the empirical per-bitstring frequencies, INCLUDING
out-of-sector strings so the ~23.5% post-selection survival is itself resampled),
re-pool, re-post-select -> a bootstrap support S*, recompute Delta_0*, Delta_1*.
B repetitions give the sampling distribution.

Fast falsifier (no eigendecomposition, exact): since U_m is unitary,
  m0_bar = || a_S ||^2 = || phi_S ||^2                       (sum of |amp|^2 on S)
  m1_bar = phi_S^dag (H_S - E0) phi_S = v^dag H_add v - E0 ||v||^2
with v = phi[si] masked to the sampled support. Validated against the stored
falsifier delta0/delta1 to <1e-9 before trusting the bootstrap.

Reuses the exact machinery (H, phi, si, postselect) by importing the job module
(import does NOT submit to QPU: the RUN=1 path is under __main__/main()).
"""
import os, sys, json, glob, importlib.util
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
JOB  = os.path.join(HERE, "hardware_matched_job_L8.py")
RES  = os.path.normpath(os.path.join(HERE, "..", "data"))
B    = int(os.environ.get("NBOOT", 1000))
SEED = int(os.environ.get("SEED", 1))

# --- import the job module (builds H, si, phi, postselect; no QPU) -------------
spec_ = importlib.util.spec_from_file_location("hwjob_L8", JOB)
mod = importlib.util.module_from_spec(spec_)
print(f"[boot] importing job machinery (builds L={os.environ.get('L','8')} H, diagonalizes sectors) ...")
spec_.loader.exec_module(mod)

si     = np.asarray(mod.si)
phi_si = np.asarray(mod.phi)[si]                 # exact classical amplitudes on the addition sector
H_add  = mod.H[si][:, si].tocsr()               # sparse addition-sector Hamiltonian
pos    = mod.pos
E0     = float(mod.E0)
m0_op  = float(mod.m0_op)
m1_op  = float(mod.m1_op)
postselect = mod.postselect
Nsi = len(si)
print(f"[boot] machinery ready: |addition sector|={Nsi}, m0_op={m0_op:.6f}, m1_op={m1_op:.6f}, E0={E0:.6f}")


def falsifier_fast(seen):
    """Exact (m0_bar,m1_bar) via norms/quadratic form -- no eigh. seen = set of Fock ints."""
    Sidx = np.fromiter((pos[x] for x in seen if x in pos), dtype=int)
    v = np.zeros(Nsi, dtype=phi_si.dtype)
    v[Sidx] = phi_si[Sidx]
    m0_bar = float(np.vdot(v, v).real)
    m1_bar = float(np.vdot(v, H_add @ v).real) - E0 * m0_bar
    return len(Sidx), m0_bar, abs(m0_op - m0_bar), abs(m1_op - m1_bar)


def pooled_from(per_time_counts):
    pooled = {}
    for cnt in per_time_counts:
        for bs, ct in cnt.items():
            pooled[bs] = pooled.get(bs, 0) + ct
    return pooled


def main():
    rng = np.random.default_rng(SEED)
    files = sorted(glob.glob(os.path.join(RES, "heron_counts_matched_L8_*.json")),
                   key=lambda p: -json.load(open(p, encoding="utf-8")).get("shots", 0))
    out = []
    for fp in files:
        d = json.load(open(fp, encoding="utf-8"))
        if d.get("_dry_run"):
            continue
        ptc = d["per_time_counts"]; shots = d["shots"]

        # ---- point estimate + self-check vs stored JSON ----
        seen, kept, total = postselect(pooled_from(ptc))
        Ssz, m0b, d0, d1 = falsifier_fast(seen)
        d0_json = d["falsifier"]["delta0"]; d1_json = d["falsifier"]["delta1"]
        assert abs(d0 - d0_json) < 1e-9 and abs(d1 - d1_json) < 1e-9, \
            f"fast falsifier mismatch vs JSON at {shots}: d0 {d0} vs {d0_json}, d1 {d1} vs {d1_json}"

        # ---- pre-extract per-circuit empirical distributions (all strings) ----
        circ = []
        for cnt in ptc:
            keys = list(cnt.keys())
            c = np.array([cnt[k] for k in keys], dtype=np.float64)
            n = int(c.sum()); p = c / c.sum()
            p = p / p.sum()                       # guard against fp round-off > 1
            circ.append((keys, p, n))

        # ---- bootstrap ----
        d0s = np.empty(B); d1s = np.empty(B); Ss = np.empty(B)
        for b in range(B):
            pooled = {}
            for keys, p, n in circ:
                nc = rng.multinomial(n, p)
                for j in np.nonzero(nc)[0]:
                    pooled[keys[j]] = pooled.get(keys[j], 0) + int(nc[j])
            seenb, _, _ = postselect(pooled)
            Sb, _, d0b, d1b = falsifier_fast(seenb)
            d0s[b] = d0b; d1s[b] = d1b; Ss[b] = Sb

        row = dict(
            shots=shots, job_id=d["job_id"], backend=d["backend"],
            S_point=int(Ssz), coverage=Ssz / Nsi,
            delta0_point=d0, delta1_point=d1,
            delta0_boot_mean=float(d0s.mean()), delta0_boot_std=float(d0s.std(ddof=1)),
            delta0_ci95=[float(np.percentile(d0s, 2.5)), float(np.percentile(d0s, 97.5))],
            delta1_boot_mean=float(d1s.mean()), delta1_boot_std=float(d1s.std(ddof=1)),
            delta1_ci95=[float(np.percentile(d1s, 2.5)), float(np.percentile(d1s, 97.5))],
            S_boot_mean=float(Ss.mean()), S_boot_std=float(Ss.std(ddof=1)),
        )
        out.append(row)
        print(f"\nshots={shots:>6}  |S|={Ssz:>4} ({100*Ssz/Nsi:.1f}%)  [self-check OK vs JSON]")
        print(f"   Delta0 = {d0:.4f}  boot mean {d0s.mean():.4f} +/- {d0s.std(ddof=1):.4f}  "
              f"95%CI [{np.percentile(d0s,2.5):.4f}, {np.percentile(d0s,97.5):.4f}]")
        print(f"   Delta1 = {d1:.4f}  boot mean {d1s.mean():.4f} +/- {d1s.std(ddof=1):.4f}  "
              f"95%CI [{np.percentile(d1s,2.5):.4f}, {np.percentile(d1s,97.5):.4f}]")
        print(f"   shot-noise dispersion / signal = {d0s.std(ddof=1)/max(d0,1e-12):.3f}  "
              f"(small => Delta0 dominated by deterministic coverage leak, not shot noise)")

    outpath = os.path.join(RES, "bootstrap_device_curve.json")
    json.dump({"B": B, "seed": SEED, "m0_op": m0_op, "m1_op": m1_op, "rows": out},
              open(outpath, "w", encoding="utf-8"), indent=2)
    print(f"\n[boot] wrote {outpath}")


if __name__ == "__main__":
    main()
