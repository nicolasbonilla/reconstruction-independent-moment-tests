# -*- coding: utf-8 -*-
"""Recount Fig. `separating` from the SEALED record, and make the figure reproducible.

WHY THIS EXISTS. The figure's data files had no committed generator. This script derives
the figure's quantities from the sealed artifacts alone (blind_verdicts.json +
blind_labels_sealed.json); all 300 instances reproduce the sealed record. The export of
paper/figs/sep_*.dat from this recount is pending.

THE NORMALISATION, recovered from blind_classify.py:23-25 rather than guessed:
    g[k]     = |m_hat[k] - m_bar[k]|                 the RAW residual
    delta[k] = z * sqrt(var/Ns + (bias*|m_hat[k]|)^2) the per-moment threshold
    a moment fires iff g[k] > delta[k]
so the normalised residual the caption describes is r[k] = g[k]/delta[k], and "> 1 is
exactly a fire" holds by construction.

    x = r[0]                      the diagonal m_0 channel (hardware-executable)
    y = max(r[1], r[2])           the larger of the two off-diagonal channels

The shaded region is x < 1 (m_0 blind) and y > 1 (an off-diagonal moment fires).

Inputs are the two sealed files only:
    data/blind_verdicts.json       (carries the prereg SHA-256 it was scored under)
    data/blind_labels_sealed.json  (the class labels, unsealed after scoring)

Run: python separating_counts.py
"""
import json
import os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "data"))


def load():
    vj = json.load(open(os.path.join(RES, "blind_verdicts.json")))
    V = {v["id"]: v for v in vj["verdicts"]}
    L = json.load(open(os.path.join(RES, "blind_labels_sealed.json")))["labels"]
    return V, L, vj.get("_prereg_sha256", "")


def norm(v):
    """(x, y, r) with r the per-moment threshold-normalised residuals."""
    r = [v["g"][k] / v["delta"][k] for k in range(3)]
    return r[0], max(r[1], r[2]), r


def classify(lab):
    m = lab["mode"]
    if m == "ac":
        return "ac_sep" if lab["params"].get("preserve_m0") else "ac_other"
    return m


def main():
    V, L, sha = load()
    cls = defaultdict(list)
    for lab in L:
        cls[classify(lab)].append(V[lab["id"]])

    print("=" * 70)
    print("Fig. `separating` -- recount from the sealed record")
    print(f"prereg SHA-256: {sha[:16]}...")
    print("=" * 70)
    print(f"{'class':<10}{'n':>6}{'in shaded region':>20}")
    region = []
    for k in sorted(cls):
        inr = [v for v in cls[k] if norm(v)[0] < 1.0 and norm(v)[1] > 1.0]
        region += [(k, v) for v in inr]
        print(f"{k:<10}{len(cls[k]):>6}{len(inr):>20}")
    print(f"{'TOTAL':<10}{sum(len(x) for x in cls.values()):>6}{len(region):>20}")

    m1 = sum(1 for _, v in region if norm(v)[2][1] > 1.0)
    m2only = sum(1 for _, v in region if norm(v)[2][2] > 1.0 and norm(v)[2][1] <= 1.0)
    print(f"\nof the {len(region)} in the region:  m_1 fires on {m1}, "
          f"m_2 alone on {m2only}  (sum {m1 + m2only})")

    ac = sum(1 for k, _ in region if k.startswith("ac"))
    print(f"composition: {ac} spurious-feature (ac), "
          f"{sum(1 for k, _ in region if k == 'clean')} clean false positive, "
          f"{sum(1 for k, _ in region if k == 'trunc')} determinant-truncation")

    tr = [v for v in cls["trunc"]]
    above = sum(1 for v in tr if norm(v)[0] > 1.0)
    print(f"\ndeterminant truncation: {len(tr)} instances, {above} with m_0 residual > 1, "
          f"{len(tr) - above} caught only off-diagonally")
    return len(region), m1, m2only


if __name__ == "__main__":
    main()
