# -*- coding: utf-8 -*-
"""Recount Fig. `separating` from the SEALED blinded record and export its plotted data.

Inputs (sealed, read-only):
    data/blind_verdicts.json       the frozen verdicts (g, delta per moment; carries the prereg SHA-256)
    data/blind_labels_sealed.json  the class labels (unsealed after scoring)
    data/prereg.sha256             checked against the SHA-256 the verdicts carry
Every one of the 300 plotted points is computed here from those two files; the frozen rule's own
fire flags are re-derived (r_k > 1) and checked against the sealed `moment_fire` for all 300.

THE NORMALISATION, taken from blind_classify.py:23-25:
    g[k]     = |m_hat[k] - m_bar[k]|                   the raw residual
    delta[k] = z * sqrt(var/Ns + (bias*|m_hat[k]|)^2)  the frozen per-moment threshold
                                                       (z = 1.96, Ns = 5e4, bias = 2%)
    r[k]     = g[k] / delta[k]                         so r[k] > 1 is exactly a fire
    x = r[0]              the diagonal total-weight channel m_0
    y = max(r[1], r[2])   the larger of the two off-diagonal residuals
The shaded region of the figure is x < 1 (m_0 does not fire) and y > 1 (an off-diagonal moment fires).
No value is clipped: the figure's axes (x in [5e-4, 60], y in [0.02, 60]) contain every point.

THE RECOUNT (sealed record, prereg SHA-256 0e4c3368...):
    26 instances lie in the shaded region =
        20 weight-preserving spurious-feature artifacts (ac, preserve_m0 = true)
      +  4 total-weight-changing spurious-feature artifacts (ac, preserve_m0 = false)
      +  1 determinant truncation (id 162; the only truncation with r_0 < 1)
      +  1 clean false positive (id 184; fires on m_2 alone).
    m_1 fires on 16 of the 26; m_2 alone on the other 10.
The earlier record data/separating_demonstration.json (24 / 15 / 9) counted the spurious-feature
instances only; it is kept unchanged. Its exemplar id 222 is weight-PRESERVING (preserve_m0 = true).

Outputs:
    paper/figs/sep_clean.dat           112 faithful reconstructions (includes id 184)
    paper/figs/sep_krylov.dat           71 under-converged Krylov
    paper/figs/sep_trunc.dat            56 determinant truncations (includes id 162)
    paper/figs/sep_ac_other.dat         37 spurious-feature instances outside the region
    paper/figs/sep_ac_pres_region.dat   20 weight-preserving spurious features in the region
    paper/figs/sep_ac_chg_region.dat     4 weight-changing spurious features in the region
    paper/figs/sep_region_nonac.dat      2 non-spurious instances in the region (ids 162, 184),
                                           drawn as an extra ring on top of their class marker
    data/2026-09-27_separating_counts.json   the counts, ids and composition (a new file)
Columns of every .dat: id x y r0 r1 r2 U (x, y as above; r_k per moment).

Deterministic (no random numbers). Run from src/:  python separating_counts.py
"""
import hashlib
import json
import os
import platform
import time
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "data"))
FIGS = os.path.normpath(os.path.join(HERE, "..", "paper", "figs"))
OUT_JSON = os.path.join(RES, "2026-09-27_separating_counts.json")
DAT_FILES = ("sep_clean", "sep_krylov", "sep_trunc", "sep_ac_other",
             "sep_ac_pres_region", "sep_ac_chg_region", "sep_region_nonac")


def sha256(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def load():
    vj = json.load(open(os.path.join(RES, "blind_verdicts.json")))
    V = {v["id"]: v for v in vj["verdicts"]}
    L = json.load(open(os.path.join(RES, "blind_labels_sealed.json")))["labels"]
    sha = vj.get("_prereg_sha256", "")
    if sha != open(os.path.join(RES, "prereg.sha256")).read().strip():
        raise SystemExit("blind_verdicts.json was not scored under the sealed prereg")
    return V, L, sha


def norm(v):
    """(x, y, r) with r the per-moment threshold-normalised residuals."""
    r = [v["g"][k] / v["delta"][k] for k in range(3)]
    return r[0], max(r[1], r[2]), r


def in_region(v):
    x, y, _ = norm(v)
    return x < 1.0 and y > 1.0


def classify(lab):
    m = lab["mode"]
    if m == "ac":
        return "ac_preserve_m0" if lab["params"].get("preserve_m0") else "ac_weight_changing"
    return m


def plot_class(lab, v):
    """which .dat file the instance's class marker is drawn from."""
    m = lab["mode"]
    if m == "ac":
        if not in_region(v):
            return "sep_ac_other"
        return "sep_ac_pres_region" if lab["params"].get("preserve_m0") else "sep_ac_chg_region"
    return {"clean": "sep_clean", "krylov": "sep_krylov", "trunc": "sep_trunc"}[m]


def write_dat(name, rows):
    with open(os.path.join(FIGS, name + ".dat"), "w", newline="\n") as f:
        f.write("id x y r0 r1 r2 U\n")
        for i, v, lab in rows:
            x, y, r = norm(v)
            f.write(f"{i} {x:.6e} {y:.6e} {r[0]:.6e} {r[1]:.6e} {r[2]:.6e} {lab['params']['U']:g}\n")


def main():
    t0 = time.time()
    V, L, sha = load()
    Lb = {lab["id"]: lab for lab in L}

    # integrity: r_k > 1 must reproduce the sealed moment_fire flag on every instance
    agree = sum((max(norm(V[i])[2]) > 1.0) == V[i]["moment_fire"] for i in V)
    if agree != len(V):
        raise SystemExit(f"normalisation does not reproduce moment_fire ({agree}/{len(V)})")

    cls = defaultdict(list)
    for lab in L:
        cls[classify(lab)].append(V[lab["id"]])

    print("=" * 70)
    print("Fig. `separating` -- recount from the sealed record")
    print(f"prereg SHA-256: {sha[:16]}...   r_k>1 reproduces moment_fire on {agree}/{len(V)}")
    print("=" * 70)
    print(f"{'class':<20}{'n':>6}{'in shaded region':>20}")
    region = []
    for k in sorted(cls):
        inr = [v for v in cls[k] if in_region(v)]
        region += [(k, v) for v in inr]
        print(f"{k:<20}{len(cls[k]):>6}{len(inr):>20}")
    print(f"{'TOTAL':<20}{sum(len(x) for x in cls.values()):>6}{len(region):>20}")

    m1 = [v["id"] for _, v in region if norm(v)[2][1] > 1.0]
    m2only = [v["id"] for _, v in region if norm(v)[2][2] > 1.0 and norm(v)[2][1] <= 1.0]
    print(f"\nof the {len(region)} in the region:  m_1 fires on {len(m1)}, "
          f"m_2 alone on {len(m2only)}  (sum {len(m1) + len(m2only)})")
    comp = defaultdict(list)
    for k, v in region:
        comp[k].append(v["id"])
    print(f"composition: {len(comp['ac_preserve_m0'])} weight-preserving + {len(comp['ac_weight_changing'])} weight-changing "
          f"spurious-feature, {len(comp['trunc'])} truncation {comp['trunc']}, "
          f"{len(comp['clean'])} clean false positive {comp['clean']}")

    tr = cls["trunc"]
    above = sum(1 for v in tr if norm(v)[0] > 1.0)
    print(f"\ndeterminant truncation: {len(tr)} instances, {above} with m_0 residual > 1, "
          f"{len(tr) - above} caught only off-diagonally")

    # ---- export the figure data ------------------------------------------------------------
    files = {n: [] for n in DAT_FILES}
    for i in sorted(V):
        lab, v = Lb[i], V[i]
        files[plot_class(lab, v)].append((i, v, lab))
        if in_region(v) and lab["mode"] != "ac":
            files["sep_region_nonac"].append((i, v, lab))
    for n in DAT_FILES:
        write_dat(n, files[n])
    nplotted = sum(len(files[n]) for n in DAT_FILES if n != "sep_region_nonac")
    assert nplotted == len(V), "every instance must be drawn exactly once as a class marker"
    print("\nwrote " + ", ".join(f"{n}.dat ({len(files[n])})" for n in DAT_FILES))

    xs = [norm(V[i])[0] for i in V]
    ys = [norm(V[i])[1] for i in V]
    ex = V[222]
    out = {
        "provenance": {
            "script": "src/separating_counts.py", "date": "2026-09-27", "seed": None,
            "note": "deterministic recount; no random numbers",
            "runtime_s": None, "python": platform.python_version(),
            "inputs_sha256": {f: sha256(os.path.join(RES, f)) for f in
                              ("blind_verdicts.json", "blind_labels_sealed.json", "prereg.sha256")},
            "code_sha256": sha256(os.path.abspath(__file__)),
            "prereg_sha256": sha},
        "normalisation": "r_k = |m_hat_k - m_bar_k| / delta_k with the frozen delta_k of blind_classify.py; "
                         "x = r_0, y = max(r_1, r_2); region = x < 1 and y > 1",
        "r_k_gt_1_reproduces_sealed_moment_fire": f"{agree}/{len(V)}",
        "n_in_region": len(region),
        "region_composition": {
            "ac_weight_preserving": len(comp["ac_preserve_m0"]), "ac_weight_changing": len(comp["ac_weight_changing"]),
            "trunc": len(comp["trunc"]), "clean_false_positive": len(comp["clean"]),
            "ids": {k: sorted(v) for k, v in comp.items()}},
        "m1_fires": len(m1), "m2_alone": len(m2only),
        "m1_fires_ids": sorted(m1), "m2_alone_ids": sorted(m2only),
        "per_class_n_and_in_region": {k: {"n": len(cls[k]), "in_region": sum(in_region(v) for v in cls[k])}
                                      for k in sorted(cls)},
        "trunc_with_m0_residual_gt_1": above, "trunc_n": len(tr),
        "range_x": [min(xs), max(xs)], "range_y": [min(ys), max(ys)],
        "exemplar_222": {"mode": Lb[222]["mode"], "params": Lb[222]["params"],
                         "r_k": norm(ex)[2], "in_region": in_region(ex)},
        "supersedes": "data/separating_demonstration.json (24/15/9: spurious-feature instances only; kept "
                      "unchanged as the earlier record)",
        "dat_files": {n + ".dat": len(files[n]) for n in DAT_FILES},
    }
    out["provenance"]["runtime_s"] = round(time.time() - t0, 3)
    json.dump(out, open(OUT_JSON, "w"), indent=1)
    print(f"wrote {OUT_JSON}")
    return len(region), len(m1), len(m2only)


if __name__ == "__main__":
    main()
