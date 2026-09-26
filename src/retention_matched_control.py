r"""D2 control: the device/noiseless comparison of Fig.(device) is at matched RAW shots.
The noiseless Aer run post-selects to 100% (particle-number-conserving Trotter circuits),
the device to ~23.5%. This script (a) records both retention fractions, and (b) re-runs the
noiseless simulation at RETENTION-MATCHED budgets so |S| is compared at equal usable samples.
"""
import os, json, glob, importlib.util
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = HERE
RES = os.path.normpath(os.path.join(HERE, "..", "data"))
JOB = os.path.join(HERE, "hardware_matched_job_L8.py")
SEED = 1

spec_ = importlib.util.spec_from_file_location("hwjob_L8", JOB)
mod = importlib.util.module_from_spec(spec_)
print("[ctrl] importing job machinery ...")
spec_.loader.exec_module(mod)

from qiskit import transpile
from qiskit_aer import AerSimulator

si = np.asarray(mod.si)
phi_si = np.asarray(mod.phi)[si]
H_add = mod.H[si][:, si].tocsr()
pos = mod.pos
E0 = float(mod.E0); m0_op = float(mod.m0_op)
postselect = mod.postselect
circuits = mod.circuits
Nsi = len(si)


def falsifier_fast(seen):
    Sidx = np.fromiter((pos[x] for x in seen if x in pos), dtype=int)
    v = np.zeros(Nsi, dtype=phi_si.dtype); v[Sidx] = phi_si[Sidx]
    m0_bar = float(np.vdot(v, v).real)
    return len(Sidx), abs(m0_op - m0_bar)


# --- device rows ---
dev = {}
for fp in glob.glob(os.path.join(RES, "heron_counts_matched_L8_*.json")):
    d = json.load(open(fp, encoding="utf-8"))
    if d.get("_dry_run"):
        continue
    f = d["falsifier"]
    dev[d["shots"]] = dict(S=f["S_size"], kept=f["kept_shots"], total=f["total_shots"],
                           ret=f["kept_shots"] / f["total_shots"], d0=f["delta0"])

base = json.load(open(os.path.join(RES, "aer_noiseless_baseline.json"), encoding="utf-8"))
nl = {r["shots"]: r for r in base["rows"]}

print("\n=== retention, as published (matched RAW shots) ===")
print(f"{'raw/circ':>9} {'dev kept':>9} {'dev ret':>8} {'nl kept':>9} {'nl ret':>7} "
      f"{'dev |S|':>8} {'nl |S|':>7} {'dev>nl?':>8}")
for NS in (50000, 30000, 16000, 4000):
    D, N = dev[NS], nl[NS]
    print(f"{NS:>9} {D['kept']:>9} {D['ret']:>8.4f} {N['kept']:>9} "
          f"{N['kept']/(NS*7):>7.4f} {D['S']:>8} {N['S']:>7} {str(D['S'] > N['S']):>8}")

# --- retention-matched noiseless re-run ---
sim = AerSimulator()
tcirc = [transpile(qc, sim) for qc in circuits]
print("\n=== retention-matched control: noiseless at the DEVICE's kept-shot count ===")
out = []
for NS in (50000, 30000, 16000, 4000):
    D = dev[NS]
    ns_eq = int(round(D["kept"] / len(tcirc)))   # noiseless keeps 100% -> raw == kept
    pooled = {}
    for k, tqc in enumerate(tcirc):
        cnt = sim.run(tqc, shots=ns_eq, seed_simulator=SEED + k).result().get_counts()
        for bs, c in cnt.items():
            pooled[bs] = pooled.get(bs, 0) + c
    seen, kept, total = postselect(pooled)
    S, d0 = falsifier_fast(seen)
    out.append(dict(raw_per_circuit=NS, device_kept=D["kept"], device_ret=D["ret"],
                    device_S=D["S"], device_d0=D["d0"],
                    noiseless_raw_matched=nl[NS]["S"], noiseless_d0_rawmatched=nl[NS]["d0_noiseless"],
                    ctrl_shots_per_circuit=ns_eq, ctrl_kept=int(kept), ctrl_S=int(S), ctrl_d0=d0))
    print(f"raw {NS:>6}/circ | device kept {D['kept']:>6} |S|={D['S']:>5} d0={D['d0']:.4f}  ||  "
          f"noiseless @ kept {kept:>6} |S|={S:>5} d0={d0:.4f}   (raw-matched noiseless |S|="
          f"{nl[NS]['S']}, d0={nl[NS]['d0_noiseless']:.4f})")

json.dump({"_provenance": {"script": "retention_matched_control.py", "seed": SEED,
                           "note": "noiseless retains 100% (particle-number-conserving Trotter "
                                   "circuits); device ~23.5% after post-selection"},
           "rows": out},
          open(os.path.join(RES, "2026-09-05_retention_matched_control.json"), "w",
               encoding="utf-8"), indent=2)
print("\n[ctrl] wrote 2026-09-05_retention_matched_control.json")
