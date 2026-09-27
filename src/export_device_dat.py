# -*- coding: utf-8 -*-
"""Write the fig_device coordinates (plan R11/R12; items B2, B3) from committed JSON only.

Reads data/2026-09-27_delta0_reference_mc.json (written by src/delta0_reference_mc.py, which recomputes
the device rows from the retained ibm_fez counts) and writes, in paper/figs/:

  device_points.dat   device rows:  shots_k delta0 S cov_pct kept retention rel_L1
  device_band.dat     noiseless sampling reference at the device's RAW shots per circuit
                      (multinomial replicas of the exact output distributions of the same 7 circuits):
                      shots_k mean lo hi sd p2_5 p50 p97_5 n     with lo/hi = mean -/+ 1 sd
  device_uniform.dat  uniform in-sector sampler at the device's pooled RETAINED count N:
                      shots_k E mc_mean N                          with E = m0 (1 - 1/3920)^N exactly

No computation beyond formatting; the numbers are the JSON's. Run after delta0_reference_mc.py.
"""
import os, json

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_JSON = os.path.normpath(os.path.join(HERE, "..", "data", "2026-09-27_delta0_reference_mc.json"))
OUT = os.path.normpath(os.path.join(HERE, "..", "paper", "figs"))
BUDGETS = (4000, 16000, 30000, 50000)          # ascending x for pgfplots


def main():
    d = json.load(open(SRC_JSON, encoding="utf-8"))
    dev, bud, nsec = d["device"], d["budgets"], d["system"]["sector_dim"]
    with open(os.path.join(OUT, "device_points.dat"), "w", encoding="utf-8") as f:
        f.write("shots_k delta0 S cov_pct kept retention rel_L1\n")
        for NS in BUDGETS:
            r = dev[str(NS)]
            f.write(f"{NS/1000:g} {r['delta0']:.6f} {r['S']} {100*r['S']/nsec:.2f} {r['kept_total']} "
                    f"{r['retention']:.4f} {r['rel_L1_stored_exact_eigh']:.4f}\n")
    with open(os.path.join(OUT, "device_band.dat"), "w", encoding="utf-8") as f:
        f.write("shots_k mean lo hi sd p2_5 p50 p97_5 n\n")
        for NS in BUDGETS:
            s = bud[str(NS)]["noiseless_raw_matched"]["delta0"]
            f.write(f"{NS/1000:g} {s['mean']:.6f} {s['mean']-s['sd']:.6f} {s['mean']+s['sd']:.6f} {s['sd']:.6f} "
                    f"{s['p2_5']:.6f} {s['p50']:.6f} {s['p97_5']:.6f} {s['n']}\n")
    with open(os.path.join(OUT, "device_uniform.dat"), "w", encoding="utf-8") as f:
        f.write("shots_k E mc_mean N\n")
        for NS in BUDGETS:
            u = bud[str(NS)]["uniform_in_sector_at_device_kept"]
            f.write(f"{NS/1000:g} {u['analytic_E_delta0']:.6e} {u['delta0']['mean']:.6e} {u['N_kept']}\n")
    for NS in BUDGETS:
        s = bud[str(NS)]["noiseless_raw_matched"]["delta0"]; u = bud[str(NS)]["uniform_in_sector_at_device_kept"]
        print(f"{NS:>6}: device {dev[str(NS)]['delta0']:.4f} | noiseless raw {s['mean']:.4f} +- {s['sd']:.4f} "
              f"| uniform in-sector E {u['analytic_E_delta0']:.2e}")
    print("wrote device_points.dat, device_band.dat, device_uniform.dat ->", OUT)


if __name__ == "__main__":
    main()
