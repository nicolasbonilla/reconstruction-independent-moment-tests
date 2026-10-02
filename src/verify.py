# -*- coding: utf-8 -*-
"""Fast reproducibility smoke test (numpy/scipy only, seconds).

Recomputes, from the exact sector-Lanczos engine, the current-probe spectral moments of the doped
L=6, U/t=4 Hubbard ring and the shot-budget thresholds tau_k = z*delta_k at N_s = 5x10^4 shots (the
per-circuit budget of the companion's ibm_fez run) with an ASSUMED 2% bias budget. It then scans
determinant truncations of the exact ground state for one that the lone first moment passes and the
second moment rejects, and checks that the scan lands on the case quoted in the manuscript (d=97).
That is one catch by m2 of a lone-m1 miss, not a guarantee of power: the every-d scan in
interval_battery.py finds 12 truncations that the whole (m0, m1, m2)+Hankel battery passes at 2% bias,
and within_sector_lp.py has redistributions that pass it by construction.

Truncation order (2026-09-27): the same DETERMINISTIC key as the interval scripts,
np.lexsort((index, -round(|psi0|^2, 12))): descending |psi0|^2, ties broken by ascending sector index.
|psi0|^2 comes in exactly degenerate 6/12/24-fold shells, and the cut at d=97 falls inside a 24-fold
shell (7 of 24 kept), so the gaps depend on which shell members are kept. The earlier np.argsort order
broke those ties by floating-point noise and printed |g1| = 0.155, |g2| = 2.560 on the authors' machine;
the deterministic order gives the manuscript's |g1| = 0.186, |g2| = 2.729. The m1-pass / m2-reject
verdict at d=97 holds for every tie choice tested (300 random shell choices, interval_battery.py ->
data/2026-09-27_interval_battery.json, tie_break_sensitivity).

Record check (record_check below, JSON reads only): every number paper/main.tex prints from the records
data/2026-08-28_inverse_moment_falsifier.json, data/heron_spectral.json, data/2026-09-27_theory_numerics.json
(keys R8_inverse_moment_sensitivity, R9_christoffel_max_Wn, m2_rescaled_frame), data/2026-08-31_offdiag_gsurface.json
(forecast percentages at 100 and 295 CZ) and paper/figs/bracket_maxW.dat is matched against its key, at the precision the paper prints it. Each check is anchored on the sentence that
prints the number; a sentence absent from main.tex is reported as ABSENT (not a failure), a key missing from its
record or a printed number that differs from its key is a failure. All other numbers: docs/REPRODUCE.md, Key
numbers.

Run:  python src/verify.py    (or `make verify`)
"""
import os
import re
import sys
import json
import numpy as np
import scipy.sparse as sp
import spectral_lanczos as sl

TOL = 1e-3
ROUND = 12            # decimals that identify the exactly degenerate |psi0|^2 shells (as in interval_battery.py)
FAIL = []


def approx(name, got, want, tol=TOL):
    ok = abs(got - want) <= tol * max(1.0, abs(want))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<40} = {got:.5f}  (expected {want:.5f})")
    if not ok:
        FAIL.append(name)
    return ok


def exact(name, got, want):
    ok = got == want
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<40} = {got}  (expected {want})")
    if not ok:
        FAIL.append(name)
    return ok


def main():
    L, U, Ns, z, bias = 6, 4.0, 50_000, 1.96, 0.02
    print(f"Reconstruction-independent moment tests -- reproducibility check (L={L}, U/t={U:.0f})\n")

    # --- exact ground state of the doped 2/3-filled sector + current probe ---
    nup = nd = L // 3
    Tu, Su, _ = sl.hop(L, nup); Td, Sd, _ = sl.hop(L, nd)
    Du, Dd = len(Su), len(Sd)
    upocc, dnocc = sl.occ_matrix(Su, L), sl.occ_matrix(Sd, L)
    Hs = (sp.kron(Tu, sp.identity(Dd)) + sp.kron(sp.identity(Du), Td)
          + sp.diags(U * (upocc @ dnocc.T).ravel())).tocsr()
    Ju, Jd = sl.current_string(L, nup), sl.current_string(L, nd)
    Js = (sp.kron(Ju, sp.identity(Dd)) + sp.kron(sp.identity(Du), Jd)).tocsr()
    E0, psi0 = sl._sub_gs(Hs)
    Hc = (Hs - E0 * sp.identity(Hs.shape[0])).tocsr()
    Jp = Js @ psi0

    # --- reconstruction-independent moments m_k = <0| J (H-E0)^k J |0> (exact state, no sampling) ---
    Hphi = Hc @ Jp
    m0 = float(np.vdot(Jp, Jp).real)
    m1 = float(np.vdot(Jp, Hphi).real)
    m2 = float(np.vdot(Jp, Hc @ Hphi).real)
    print("Reconstruction-independent current moments (operator algebra, never read off A(w)):")
    approx("m0 = <J^2>", m0, 0.5542)
    approx("m1 = 1/2<[J,[H,J]]>", m1, 4.2271)
    approx("m2 = <J (H-E0)^2 J>", m2, 33.3887)

    # --- shot-budget thresholds (local-estimator variance = a LOWER bound on the per-shot variance;
    #     2% bias budget ASSUMED, not measured) ---
    p = np.abs(psi0) ** 2
    mask = np.abs(psi0) > 1e-12

    def locvar(Opsi):     # local estimator O_k^loc(x) = (O_k psi0)_x / psi0_x; variance over p_x = |psi0_x|^2
        loc = np.zeros(len(psi0), dtype=np.result_type(Opsi, psi0)); loc[mask] = Opsi[mask] / psi0[mask]
        loc = np.real(loc)    # the imaginary part is identically zero here
        return float(np.sum(p * loc ** 2) - np.sum(p * loc) ** 2)

    var0, var1, var2 = locvar(Js @ Jp), locvar(Js @ Hphi), locvar(Js @ (Hc @ Hphi))
    tau0 = z * np.sqrt(var0 / Ns + (bias * abs(m0)) ** 2)
    tau1 = z * np.sqrt(var1 / Ns + (bias * abs(m1)) ** 2)
    tau2 = z * np.sqrt(var2 / Ns + (bias * abs(m2)) ** 2)
    print(f"\nShot-budget thresholds tau_k = z*delta_k (Ns={Ns}, z={z}, assumed {100*bias:.0f}% bias budget):")
    approx("tau0", tau0, 0.02868)
    approx("tau1", tau1, 0.22894)
    approx("tau2", tau2, 1.84760)
    print(f"  (tau1 = {100*tau1/m1:.1f}% of m1)")

    # --- deterministic truncation order (same key as interval_battery.py / interval_moment.py) ---
    key = np.round(p, ROUND)
    order = np.lexsort((np.arange(len(p)), -key))
    nsup = int((p > 1e-14).sum())

    def trunc_moments(d):
        idx = np.sort(order[:d])
        Hsub, Jsub = Hs[idx][:, idx], Js[idx][:, idx]
        e0, g = sl._sub_gs(Hsub); Hcs = (Hsub - e0 * sp.identity(d)).tocsr(); Jg = Jsub @ g
        return (float(np.vdot(Jg, Jg).real), float(np.vdot(Jg, Hcs @ Jg).real),
                float(np.vdot(Jg, Hcs @ (Hcs @ Jg)).real))

    def hankel_ok(a, b, c):   # point check on the reconstruction's own moments (cannot fire for a positive measure)
        return bool(np.linalg.eigvalsh(np.array([[a, b], [b, c]]))[0] >= -1e-9 and b >= -1e-9)

    # --- scan (auto-calibrated, no hardcoded d): from d = nsup-1 downward, the first truncation that the
    #     lone first moment passes (|g1| <= tau1) and the second moment rejects (|g2| > tau2) ---
    demo = None
    for d in range(nsup - 1, max(2, nsup // 3), -1):
        b0, b1, b2 = trunc_moments(d)
        g0, g1, g2 = abs(b0 - m0), abs(b1 - m1), abs(b2 - m2)
        if g1 <= tau1 and g2 > tau2:
            demo = (d, b0, b1, b2, g0, g1, g2)
            break

    if demo is None:
        print(f"\n  [FAIL] no truncation in d = {nsup - 1}..{max(2, nsup // 3) + 1} is passed by m1 alone and "
              f"rejected by m2 (expected d=97; see interval_battery.py)")
        FAIL.append("scan")
    else:
        d, b0, b1, b2, g0, g1, g2 = demo
        pc = key[order[d - 1]]
        shell, above = int(np.sum(key == pc)), int(np.sum(key > pc))
        battery_reject = (g0 > tau0) or (g1 > tau1) or (g2 > tau2) or (not hankel_ok(b0, b1, b2))
        print(f"\nTruncation found by the scan (deterministic lexsort order), n_support = {nsup} determinants:")
        exact("truncation depth d", d, 97)
        approx("|m1_trunc - m1| (passes: <= tau1)", g1, 0.18600)
        approx("|m2_trunc - m2| (rejects: > tau2)", g2, 2.72896)
        print(f"  |m0_trunc - m0| = {g0:.4f} vs tau0 = {tau0:.4f} -> {'reject' if g0 > tau0 else 'pass'}; "
              f"Hankel point check on the truncated moments: {'pass' if hankel_ok(b0, b1, b2) else 'reject'}")
        print(f"  the cut keeps {d - above} of a {shell}-fold degenerate |psi0|^2 shell ({above} determinants above it);"
              f" tie-robustness over 300 shell choices: interval_battery.py")
        if not battery_reject:
            FAIL.append("battery")
        print(f"\n  [{'PASS' if battery_reject else 'FAIL'}] the joint (m0,m1,m2)+Hankel battery rejects, via m2, a "
              f"truncation the lone first moment passes (one catch, not a guarantee of power)")

    record_check()

    print("\n" + ("=" * 62))
    if FAIL:
        print(f"RESULT: {len(FAIL)} CHECK(S) FAILED: {', '.join(FAIL)}")
        sys.exit(1)
    print("RESULT: ALL CHECKS PASS -- the L=6 moments, thresholds and the d=97 catch reproduce from source, and the\n"
          "numbers main.tex prints from the records listed in record_check match their JSON keys (other numbers:\n"
          "docs/REPRODUCE.md, Key numbers).")


# =====================================================================================================
# record_check: numbers printed in paper/main.tex against their JSON keys (no recomputation, seconds)
# =====================================================================================================
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))


def _load(rel):
    with open(os.path.join(ROOT, rel), encoding="utf-8") as f:
        return json.load(f) if rel.endswith(".json") else f.read()


def _same_as_printed(printed, value):
    """True if `value` rounds to the printed decimal string at the printed number of decimals."""
    nd = len(printed.split(".")[1]) if "." in printed else 0
    return f"{value:.{nd}f}" == printed


def _sci_as_printed(mant, expo, value):
    """True if `value` written to one significant figure is mant x 10^expo (printed as $m\\times10^{e}$)."""
    m, e = f"{value:.0e}".split("e")
    return int(m) == int(mant) and int(e) == int(expo)


def _maxW_dat():
    rows = [r.split() for r in _load("paper/figs/bracket_maxW.dat").strip().splitlines()[1:]]
    return {int(n): float(w) for n, w in rows}


def _records():
    """(id, record/key, regex on main.tex, function(match groups) -> list of (printed, value, kind))."""
    inv = _load("data/2026-08-28_inverse_moment_falsifier.json")
    dep, summ = inv["deployment"], inv["summary"]
    tn = _load("data/2026-09-27_theory_numerics.json")
    r8, r9, m2f = tn["R8_inverse_moment_sensitivity"], tn["R9_christoffel_max_Wn"], tn["m2_rescaled_frame"]
    ring, opn = r8["ring"], r8["open_chain"]
    her = _load("data/heron_spectral.json")
    last = dep["sweep"][-1]
    exmax = {int(n): v["max_W"] for n, v in r9["exact_max"].items()}
    obc = [c for c in inv["identity_checks"] if not c["pbc"]]
    obc_rel = max(c["resid_identity_halfnegT"] / abs(c["m_minus1"]) for c in obc)

    def first_poles(geo):   # spurious-pole positions at which m_-1 crosses its threshold first (full rule)
        return [s["om_spur"] for s in geo["om_spur_scan"]
                if s["f_star_full_rule"]["m_-1"] < min(s["f_star_full_rule"]["m1"], s["f_star_full_rule"]["m2"])]

    P, X = "p", "sci"
    INV = "2026-08-28_inverse_moment_falsifier.json"
    R8 = "2026-09-27_theory_numerics.json:R8_inverse_moment_sensitivity"
    R9 = "2026-09-27_theory_numerics.json:R9_christoffel_max_Wn"
    HER = "heron_spectral.json"
    gs = {r["G"]: r for r in _load("data/2026-08-31_offdiag_gsurface.json")["gsurface"]}
    GS = "2026-08-31_offdiag_gsurface.json:gsurface[G]"
    return [
        ("inv-relshift-text", INV + ":summary.max_Delta_{-1,1}_over_m",
         r"moves most in relative terms \(\$(\d+)\\%\$ against \$(\d+)\\%\$ for \$m_1\$",
         lambda g: [(g[0], 100 * summ["max_Delta_minus1_over_m"], P), (g[1], 100 * summ["max_Delta_1_over_m"], P)]),
        ("inv-relshift-caption", INV + ":deployment.sweep[-1].{d1,d2,dm1}",
         r"move by up to \$(\d+)\\%\$, \$(\d+)\\%\$ and \$(\d+)\\%\$",
         lambda g: [(g[0], 100 * last["d1"], P), (g[1], 100 * last["d2"], P), (g[2], 100 * last["dm1"], P)]),
        ("inv-m-1-caption", INV + ":deployment.m_minus1_poles",
         r"m_\{-1\}=\\tfrac12\\langle-\\hat T\\rangle-\\mathcal D=([\d.]+)",
         lambda g: [(g[0], dep["m_minus1_poles"], P)]),
        ("inv-m0-blind", INV + ":deployment.sweep[*].d0",
         r"\\Delta_0/m_0<10\^\{(-?\d+)\}",
         lambda g: [(g[0], max(abs(r["d0"]) for r in dep["sweep"]) < 10.0 ** int(g[0]), "bool")]),
        ("inv-cancellation", INV + ":deployment.{m_minus1_poles,half_negT,kohn_stiffness_D}",
         r"\$m_\{-1\}=([\d.]+)\$ is the small difference of \$\\tfrac12\\langle-\\hat T\\rangle=([\d.]+)\$ and "
         r"\$\\mathcal D=([\d.]+)\$: the stiffness cancels \$(\d+)\\%\$",
         lambda g: [(g[0], dep["m_minus1_poles"], P), (g[1], dep["half_negT"], P), (g[2], dep["kohn_stiffness_D"], P),
                    (g[3], 100 * dep["kohn_stiffness_D"] / dep["half_negT"], P)]),
        ("inv-cancellation-caption", INV + ":deployment.{kohn_stiffness_D,half_negT}",
         r"cancels \$(\d+)\\%\$ of \$\\tfrac12\\langle-\\hat T\\rangle\$",
         lambda g: [(g[0], 100 * dep["kohn_stiffness_D"] / dep["half_negT"], P)]),
        ("inv-ring-identity", INV + ":deployment.resid_fsum_vs_poles",
         r"verified at \$L=12\$ to a relative \$(\d)\\times10\^\{(-\d+)\}\$",
         lambda g: [((g[0], g[1]), dep["resid_fsum_vs_poles"], X)]),
        ("inv-obc-identity", INV + ":identity_checks[OBC] + " + R8 + ".open_chain.resid_poles_vs_fsum_rel",
         r"to a relative \$(\d)\\times10\^\{(-\d+)\}\$ at \$L=6\$ and \$8\$ and \$(\d)\\times10\^\{(-\d+)\}\$ at \$L=12\$",
         lambda g: [((g[0], g[1]), obc_rel, X), ((g[2], g[3]), opn["resid_poles_vs_fsum_rel"], X)]),
        ("r8-thresholds", R8 + ".ring.{tau_full_rule,tau_full_rule_over_m}",
         r"\\tau_\{-1\}=([\d.]+)\$, which is \$(\d+)\\%\$ of \$m_\{-1\}\$, whereas \$\\tau_1\$ and \$\\tau_2\$ are "
         r"\$([\d.]+)\\%\$",
         lambda g: [(g[0], ring["tau_full_rule"]["m_-1"], P), (g[1], 100 * ring["tau_full_rule_over_m"]["m_-1"], P),
                    (g[2], 100 * ring["tau_full_rule_over_m"]["m1"], P),
                    (g[2], 100 * ring["tau_full_rule_over_m"]["m2"], P)]),
        ("r8-thresholds-caption", R8 + ".ring.tau_full_rule_over_m",
         r"(?:Deployed|Default) thresholds: \$([\d.]+)\\%\$ of \$m_1\$ and \$m_2\$ but \$(\d+)\\%\$ of \$m_\{-1\}\$",
         lambda g: [(g[0], 100 * ring["tau_full_rule_over_m"]["m1"], P),
                    (g[0], 100 * ring["tau_full_rule_over_m"]["m2"], P),
                    (g[1], 100 * ring["tau_full_rule_over_m"]["m_-1"], P)]),
        ("r8-ring-order", R8 + ".ring.f_star_full_rule",
         r"So \$m_2\$ rejects first \(\$f\\approx([\d.]+)\$\), then \$m_1\$ \(\$([\d.]+)\$\), and \$m_\{-1\}\$ only at "
         r"\$f\\approx([\d.]+)\$",
         lambda g: [(g[0], ring["f_star_full_rule"]["m2"], P), (g[1], ring["f_star_full_rule"]["m1"], P),
                    (g[2], ring["f_star_full_rule"]["m_-1"], P)]),
        ("r8-ring-order-caption", R8 + ".ring.f_star_full_rule",
         r"rejection at \$f\\approx([\d.]+)\$ \(\$m_2\$\), \$([\d.]+)\$ \(\$m_1\$\), \$([\d.]+)\$ \(\$m_\{-1\}\$\)",
         lambda g: [(g[0], ring["f_star_full_rule"]["m2"], P), (g[1], ring["f_star_full_rule"]["m1"], P),
                    (g[2], ring["f_star_full_rule"]["m_-1"], P)]),
        ("r8-open-chain", R8 + ".open_chain.{frac_below_split,m_minus1_poles}",
         r"\$(\d+)\\%\$ of its regular weight already lies below \$6\\,t\$ and \$m_\{-1\}=([\d.]+)\$",
         lambda g: [(g[0], 100 * opn["frac_below_split"], P), (g[1], opn["m_minus1_poles"], P)]),
        ("r8-open-order", R8 + ".open_chain.f_star_full_rule",
         r"crossed at \$f\\approx([\d.]+)\$ for \$m_\{-1\}\$, \$([\d.]+)\$ for \$m_1\$ and \$([\d.]+)\$ for \$m_2\$",
         lambda g: [(g[0], opn["f_star_full_rule"]["m_-1"], P), (g[1], opn["f_star_full_rule"]["m1"], P),
                    (g[2], opn["f_star_full_rule"]["m2"], P)]),
        ("r8-low-pole", R8 + ".ring.om_spur_scan[om_spur=0.25].f_star_full_rule",
         r"at \$0\.25\\,t\$ on the ring, \$f\\approx([\d.]+)\$ against \$([\d.]+)\$ for \$m_1\$ and \$([\d.]+)\$ for \$m_2\$",
         lambda g: [(g[0], ring["om_spur_scan"][0]["f_star_full_rule"]["m_-1"], P),
                    (g[1], ring["om_spur_scan"][0]["f_star_full_rule"]["m1"], P),
                    (g[2], ring["om_spur_scan"][0]["f_star_full_rule"]["m2"], P),
                    ("0.25", ring["om_spur_scan"][0]["om_spur"], P)]),
        ("r8-first-poles", R8 + ".{ring,open_chain}.om_spur_scan",
         r"rejects first only for poles at \$([\d.]+)\\,t\$ and \$([\d.]+)\\,t\$, on both geometries",
         lambda g: [(f"{float(g[0])}/{float(g[1])}", first_poles(ring) == first_poles(opn) == [float(g[0]), float(g[1])],
                     "bool")]),
        ("r8-open-slopes", R8 + ".open_chain.slope_relative_per_unit_f",
         r"relative slopes? \$([\d.]+)\$ for \$m_\{-1\}\$ (?:vs|against) \$([\d.]+)\$ for \$m_1\$",
         lambda g: [(g[0], opn["slope_relative_per_unit_f"]["m_-1"], P), (g[1], opn["slope_relative_per_unit_f"]["m1"], P)]),
        ("r9-maxW-caption", R9 + ".exact_max + paper/figs/bracket_maxW.dat",
         r"\\max_xW_n=([\d.]+)\$, \$([\d.]+)\$, \$([\d.]+)\$, \$([\d.]+)\$ for \$n=1,\\dots,4\$",
         lambda g: [(g[i], exmax[i + 1], P) for i in range(4)]),
        ("r9-maxW-appendix", R9 + ".exact_max",
         r"\\max_xW_n=([\d.]+),\\ ([\d.]+),\\ ([\d.]+),\\ ([\d.]+)\$ for \$n=1,\\dots,4\$",
         lambda g: [(g[i], exmax[i + 1], P) for i in range(4)]),
        ("r9-maxW4", R9 + ".exact_max.4.max_W",
         r"\\max_xW_4=([\d.]+)",
         lambda g: [(g[0], exmax[4], P)]),
        ("r9-dominant-atom", R9 + ".{dominant_atom_omega,dominant_residue}",
         r"dominant atom \(\$\\omega=([\d.]+)\\,t\$, residue \$([\d.]+)\$\)",
         lambda g: [(g[0], r9["dominant_atom_omega"], P), (g[1], r9["dominant_residue"], P)]),
        ("r9-residue-floor", R9 + ".dominant_residue",
         r"equals its residue \$([\d.]+)\$, below which",
         lambda g: [(g[0], r9["dominant_residue"], P)]),
        ("r9-one-pole", R9 + ".dominant_residue",
         r"one pole carries \$(\d+)\\%\$ of its weight",
         lambda g: [(g[0], 100 * r9["dominant_residue"], P)]),
        ("r9-mean-variance", R9 + ".omega_mean + m2_rescaled_frame.variance",
         r"mean \$\\bar\\omega=([\d.]+)\\,t\$ and variance \$\\sigma\^2=([\d.]+)\\,t\^2\$",
         lambda g: [(g[0], r9["omega_mean"], P), (g[1], m2f["variance"], P)]),
        ("r9-mean-offset", R9 + ".{omega_mean,dominant_atom_omega}",
         r"lies at the mean, \$([\d.]+)\\,t\$ from the dominant atom \(\$([\d.]+)\\,t\$\)",
         lambda g: [(g[0], abs(r9["omega_mean"] - r9["dominant_atom_omega"]), P), (g[1], r9["dominant_atom_omega"], P)]),
        ("r9-mean-appendix", R9 + ".omega_mean",
         r"mean frequency \$\\bar\\omega=([\d.]+)\\,t\$",
         lambda g: [(g[0], r9["omega_mean"], P)]),
        ("offdiag-forecast-app", GS + ".{pct_exact,pct_folklore,pct_coherent} at G=100, 295",
         r"bias at \$(\d+)\\%\$ of the signal at \$100\$ controlled-\$Z\$ \(CZ\) gates and \$(\d+)\\%\$ at \$295\$ CZ; "
         r"global depolarizing noise at the same depths gives \$(\d+)\\%\$ and \$(\d+)\\%\$, and a coherent-\$ZZ\$ model "
         r"\$(\d+)\\%\$ at \$295\$ CZ",
         lambda g: [(g[0], gs[100]["pct_exact"], P), (g[1], gs[295]["pct_exact"], P),
                    (g[2], gs[100]["pct_folklore"], P), (g[3], gs[295]["pct_folklore"], P),
                    (g[4], gs[295]["pct_coherent"], P)]),
        ("offdiag-forecast-limits", GS + ".{pct_exact,pct_folklore} at G=295",
         r"at \$295\$ CZ and two-qubit error \$2\.5\\times10\^\{-3\}\$ is \$(\d+)\\%\$ of the signal under a static "
         r"local noise model and \$(\d+)\\%\$ under global",
         lambda g: [(g[0], gs[295]["pct_exact"], P), (g[1], gs[295]["pct_folklore"], P)]),
        ("offdiag-forecast-range", GS + ".{pct_exact,pct_folklore} at G=295 (Table scope, Conclusions)",
         r"bias (?:for a surrogate probe is )?\$(\d+)\$--\$(\d+)\\%\$",
         lambda g: [(g[0], gs[295]["pct_exact"], P), (g[1], gs[295]["pct_folklore"], P)]),
        ("heron-coverage", HER + ":{hw_S,nsector}",
         r"\|\\mathcal S\|=(\d+)/(\d+)",
         lambda g: [(g[0], her["hw_S"], P), (g[1], her["nsector"], P)]),
        ("heron-sector", HER + ":{nsector,eta}",
         r"\$4\\!\\uparrow,3\\!\\downarrow\$ \(\$(\d+)\$ configurations\), \$\\eta=([\d.]+)\\,t\$",
         lambda g: [(g[0], her["nsector"], P), (g[1], her["eta"], P)]),
        ("heron-model", HER + ":{L,U}",
         r"\$L=(\d+)\$ Hubbard model, \$U/t=([\d.]+)\$, addition sector",
         lambda g: [(g[0], her["L"], P), (g[1], her["U"], P)]),
        ("heron-shots-job", HER + ":{shots,job_id,backend} (7 circuits: the companion's K=7)",
         r"\\texttt\{ibm\\_fez\}, \$(\d+)\$ qubits, \$([\d.]+)\\times10\^\{5\}\$ shots \(\$(\d+)\\times10\^\{4\}\$ on each "
         r"of seven circuits; job \\texttt\{([a-z0-9]+)\}",
         lambda g: [(g[0], 2 * her["L"], P), (g[1], 7 * her["shots"] / 1e5, P), (g[2], her["shots"] / 1e4, P),
                    (g[3], g[3] == her["job_id"] and her["backend"] == "ibm_fez", "bool")]),
        ("heron-retained", HER + ":hw_S",
         r"all \$(\d+)\$ retained determinants",
         lambda g: [(g[0], her["hw_S"], P)]),
    ], exmax


def record_check():
    print("\nRecord check: numbers printed in paper/main.tex against their JSON keys:")
    tex_path = os.path.join(ROOT, "paper", "main.tex")
    if not os.path.exists(tex_path):
        print("  [SKIP] paper/main.tex not found")
        return
    tex = _load("paper/main.tex")
    try:
        records, exmax = _records()
    except (KeyError, IndexError, OSError, ValueError) as e:
        print(f"  [FAIL] a record or key is missing: {e!r}")
        FAIL.append("record_check:keys")
        return
    # data-file consistency: bracket_maxW.dat (export_bracketing_dat.py) = R9 exact maxima (export_christoffel_dat.py)
    dat = _maxW_dat()
    ok = all(abs(dat[n] - exmax[n]) < 5e-7 for n in (1, 2, 3, 4))
    print(f"  [{'PASS' if ok else 'FAIL'}] {'bracket_maxW.dat = R9 exact max_x W_n':<40} = "
          f"{[dat[n] for n in (1, 2, 3, 4)]}")
    if not ok:
        FAIL.append("record_check:bracket_maxW")
    n_pass = n_abs = 0
    for rid, key, pat, fn in records:
        hits = re.findall(pat, tex)
        if not hits:
            n_abs += 1
            print(f"  [ABSENT] {rid:<26} (sentence not in main.tex; key {key})")
            continue
        for h in hits:
            g = h if isinstance(h, tuple) else (h,)
            for printed, value, kind in fn(g):
                if kind == "bool":
                    good, shown = bool(value), f"{printed}: {bool(value)}"
                elif kind == "sci":
                    good, shown = _sci_as_printed(printed[0], printed[1], value), \
                        f"{printed[0]}e{printed[1]} vs {value:.3e}"
                else:
                    good, shown = _same_as_printed(printed, value), f"{printed} vs {value:.6g}"
                n_pass += good
                print(f"  [{'PASS' if good else 'FAIL'}] {rid:<26} {shown}   ({key})")
                if not good:
                    FAIL.append(f"record_check:{rid}")
    print(f"  {n_pass} printed numbers match their keys; {n_abs} anchored sentences absent from main.tex")


if __name__ == "__main__":
    main()
