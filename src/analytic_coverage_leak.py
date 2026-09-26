r"""
Analytic shot-noise model of the coverage-leak falsifier Delta_0(N).

SUPERSEDED MODEL (fp_rate = 1.0); kept as a record, not evidence. Its committed output is
data/_superseded/analytic_coverage_leak.json, which is where a re-run writes.

Replaces the (upward-biased) nonparametric bootstrap of the support functional with
the CORRECT decomposition from the exact amplitudes. A determinant x of the addition
sector carries weight w_x = |phi_x|^2 (sum_x w_x = m0_op) and ideal sampling
probability q_x = w_x / m0_op (sum_x q_x = 1). After N post-selected shots it is
absent from the support S with probability (1-q_x)^N, so the missed weight is

    Delta_0(N) = sum_x w_x * 1[x not in S],     E[.] deterministic, Var[.] = shot noise:
    E[Delta_0](N)   = sum_x w_x (1-q_x)^N                     (deterministic coverage leak)
    Var[Delta_0](N) = sum_x w_x^2 (1-q_x)^N (1-(1-q_x)^N)     (Poisson/independent approx)

This separates cleanly the two contributions the bootstrap conflated. A calibrated
z-test fires when Delta_0 > tau0(N) = z * sd[Delta_0](N); tau0 is fixed by the shot
budget and the exact amplitudes, NOT by the measured outcome -> pre-registrable.
A "clean pass" at the operating budget means the measured Delta_0 sits below tau0
with a modeled false-positive rate at the nominal alpha, and the firing points sit
far above tau0.

Also inverts the model: the number of post-selected shots N* to reach a target
expected leak (=> a target coverage / clean-pass margin), which tells us whether a
higher-shot QPU run (budget renews each cycle) is worth spending, before spending it.
"""
import os, json, glob, importlib.util
import numpy as np
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
JOB  = os.path.join(HERE, "hardware_matched_job_L8.py")
RES  = os.path.normpath(os.path.join(HERE, "..", "data"))
Z    = float(os.environ.get("ZLEVEL", norm.ppf(0.975)))      # 1.96 -> two-sided 5%
ALPHA = 2 * (1 - norm.cdf(Z))

spec_ = importlib.util.spec_from_file_location("hwjob_L8", JOB)
mod = importlib.util.module_from_spec(spec_)
print("[cov] importing job machinery (builds H, diagonalizes sectors) ...")
spec_.loader.exec_module(mod)

si     = np.asarray(mod.si)
w      = np.abs(np.asarray(mod.phi)[si]) ** 2                 # weights, sum = m0_op
m0     = float(w.sum())
q      = w / m0                                               # ideal sampling probs, sum=1
print(f"[cov] sector={len(si)}  m0_op={m0:.6f}  max q={q.max():.3e}  min q>0={q[q>0].min():.3e}")


def log1m_q():
    with np.errstate(divide="ignore"):
        return np.log1p(-q)                                    # log(1-q), stable


LQ = log1m_q()


def E_delta0(N):
    return float(np.sum(w * np.exp(N * LQ)))                   # sum w_x (1-q_x)^N


def sd_delta0(N):
    miss = np.exp(N * LQ)
    var = np.sum(w ** 2 * miss * (1.0 - miss))
    return float(np.sqrt(var))


def tau0(N):
    return Z * sd_delta0(N)                                    # calibrated pre-registered threshold


def fp_rate(N):
    """Modeled false-positive prob that a re-drawn Delta_0 exceeds tau0(N) when the
    reconstruction is exact-by-coverage (Gaussian approx around E[Delta_0])."""
    mu, sd, t = E_delta0(N), sd_delta0(N), tau0(N)
    if sd == 0:
        return 0.0
    return float(1.0 - norm.cdf((t - mu) / sd))


def N_for_target_leak(target):
    lo, hi = 1.0, 5e7
    for _ in range(80):
        mid = np.sqrt(lo * hi)
        if E_delta0(mid) > target:
            lo = mid
        else:
            hi = mid
    return mid


def main():
    files = sorted(glob.glob(os.path.join(RES, "heron_counts_matched_L8_*.json")),
                   key=lambda p: -json.load(open(p, encoding="utf-8")).get("shots", 0))
    print(f"\n=== analytic coverage-leak model  (z={Z:.3f}, nominal alpha={100*ALPHA:.1f}%) ===")
    print(f"{'shots':>7} {'kept N':>8} {'meas d0':>9} {'E[d0]':>9} {'sd[d0]':>9} "
          f"{'tau0':>8} {'z-score':>8} {'verdict':>9} {'modeled FP@this N':>18}")
    rows = []
    for fp in files:
        d = json.load(open(fp, encoding="utf-8"))
        if d.get("_dry_run"):
            continue
        N = int(d["falsifier"]["kept_shots"]); meas = float(d["falsifier"]["delta0"])
        Ed, sd, t = E_delta0(N), sd_delta0(N), tau0(N)
        z = meas / sd if sd else float("inf")
        verdict = "FIRE" if meas > t else "pass"
        fpr = fp_rate(N)
        rows.append(dict(shots=d["shots"], kept=N, meas_delta0=meas, E_delta0=Ed,
                         sd_delta0=sd, tau0=t, zscore=z, verdict=verdict, fp_rate=fpr))
        print(f"{d['shots']:>7} {N:>8} {meas:>9.4f} {Ed:>9.4f} {sd:>9.4f} "
              f"{t:>8.4f} {z:>8.2f} {verdict:>9} {100*fpr:>16.2f}%")

    # inversion: shots to reach clean-pass margins
    print("\n=== shots needed (post-selected N) for target expected leak E[Delta_0] ===")
    inv = {}
    for tgt in (0.010, 0.005, 0.002, 0.001):
        Nc = N_for_target_leak(tgt)
        raw = Nc / 0.235 / 7.0                                 # ~23.5% survival, 7 circuits
        inv[tgt] = dict(N_kept=Nc, shots_per_circuit=raw)
        print(f"  E[Delta_0] <= {tgt:.3f}  ->  N_kept ~ {Nc:>9.0f}   (~{raw:>8.0f} shots/circuit x7)")

    json.dump({"z": Z, "alpha": ALPHA, "rows": rows, "inversion": inv,
               "note": "E=deterministic coverage leak, sd=shot noise, tau0=z*sd pre-registered"},
              open(os.path.join(RES, "_superseded", "analytic_coverage_leak.json"), "w", encoding="utf-8"), indent=2)
    print(f"\n[cov] wrote {os.path.join(RES,'_superseded','analytic_coverage_leak.json')}")


if __name__ == "__main__":
    main()
