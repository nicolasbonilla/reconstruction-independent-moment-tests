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

Run:  python src/verify.py    (or `make verify`)
"""
import sys
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

    print("\n" + ("=" * 62))
    if FAIL:
        print(f"RESULT: {len(FAIL)} CHECK(S) FAILED: {', '.join(FAIL)}")
        sys.exit(1)
    print("RESULT: ALL CHECKS PASS -- the screen's headline numbers reproduce from source.")


if __name__ == "__main__":
    main()
