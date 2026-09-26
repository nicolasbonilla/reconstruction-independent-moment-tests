r"""
Noiseless Aer baseline for the device-side coverage-leak falsifier.

Runs the SAME K=7 Krylov-Trotter circuits used on ibm_fez, on a noiseless Aer
simulator, at the matched shot budgets (50k/30k/16k/4k per circuit). This isolates
the finite-sampling + Trotter contribution to Delta_0 from the device-noise
contribution, giving the honest decomposition

    Delta_0(device) = Delta_0(noiseless Aer)   [sampling + Trotter]
                    + [device-noise excess]     [hardware infidelity -> support collapse]

The measured device-noise excess calibrates the bias budget b0 of the pre-registered
threshold tau0 = z*sqrt(Var_shot/N + b0^2), against which 50k passes and the
undersampled points fire -- a clean, calibrated pass rather than a marginal one.

Reuses the exact machinery + the O(1) fast falsifier (no eigh) from the bootstrap.
"""
import os, json, glob, importlib.util
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
JOB  = os.path.join(HERE, "hardware_matched_job_L8.py")
RES  = os.path.normpath(os.path.join(HERE, "..", "data"))
SEED = int(os.environ.get("SEED", 1))
LEVELS = [50000, 30000, 16000, 4000]

spec_ = importlib.util.spec_from_file_location("hwjob_L8", JOB)
mod = importlib.util.module_from_spec(spec_)
print("[aer] importing job machinery ...")
spec_.loader.exec_module(mod)

from qiskit import transpile
from qiskit_aer import AerSimulator

si     = np.asarray(mod.si)
phi_si = np.asarray(mod.phi)[si]
H_add  = mod.H[si][:, si].tocsr()
pos    = mod.pos
E0     = float(mod.E0); m0_op = float(mod.m0_op); m1_op = float(mod.m1_op)
postselect = mod.postselect
circuits = mod.circuits
Nsi = len(si)


def falsifier_fast(seen):
    Sidx = np.fromiter((pos[x] for x in seen if x in pos), dtype=int)
    v = np.zeros(Nsi, dtype=phi_si.dtype); v[Sidx] = phi_si[Sidx]
    m0_bar = float(np.vdot(v, v).real)
    m1_bar = float(np.vdot(v, H_add @ v).real) - E0 * m0_bar
    return len(Sidx), abs(m0_op - m0_bar), abs(m1_op - m1_bar)


def main():
    sim = AerSimulator()
    tcirc = [transpile(qc, sim) for qc in circuits]
    # device points keyed by shots
    dev = {}
    for fp in glob.glob(os.path.join(RES, "heron_counts_matched_L8_*.json")):
        d = json.load(open(fp, encoding="utf-8"))
        if not d.get("_dry_run"):
            dev[d["shots"]] = d["falsifier"]["delta0"]

    rows = []
    print(f"\n{'shots':>7} {'kept':>7} {'|S|':>5} {'d0_noiseless':>13} {'d0_device':>10} {'excess(dev-noise)':>18}")
    for NS in LEVELS:
        pooled = {}
        for k, tqc in enumerate(tcirc):
            cnt = sim.run(tqc, shots=NS, seed_simulator=SEED + k).result().get_counts()
            for bs, c in cnt.items():
                pooled[bs] = pooled.get(bs, 0) + c
        seen, kept, total = postselect(pooled)
        Ssz, d0n, d1n = falsifier_fast(seen)
        d0d = dev.get(NS, float("nan"))
        excess = d0d - d0n
        rows.append(dict(shots=NS, kept=int(kept), S=int(Ssz),
                         d0_noiseless=d0n, d1_noiseless=d1n, d0_device=d0d,
                         device_noise_excess=excess))
        print(f"{NS:>7} {kept:>7} {Ssz:>5} {d0n:>13.4f} {d0d:>10.4f} {excess:>18.4f}")

    b0 = float(np.median([r["device_noise_excess"] for r in rows if np.isfinite(r["device_noise_excess"])]))
    print(f"\n[aer] median device-noise excess (bias budget b0 proxy) = {b0:.4f}")
    json.dump({"seed": SEED, "levels": LEVELS, "rows": rows, "b0_proxy": b0},
              open(os.path.join(RES, "aer_noiseless_baseline.json"), "w", encoding="utf-8"), indent=2)
    print(f"[aer] wrote {os.path.join(RES,'aer_noiseless_baseline.json')}")


if __name__ == "__main__":
    main()
