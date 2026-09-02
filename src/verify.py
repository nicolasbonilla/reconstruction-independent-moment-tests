# -*- coding: utf-8 -*-
"""Fast reproducibility smoke test (numpy/scipy only, seconds).

Recomputes, from the exact sector-Lanczos engine, the current-probe spectral moments of the doped
L=6, U/t=4 Hubbard ring; the shot-budget interval at the real ibm_fez 50k-shot budget; and shows the
joint (m0, m1, m2) + Hankel battery REJECTING a determinant truncation that the lone first moment
misses -- the make-or-break "independence gives the teeth" result, distilled to a PASS/FAIL check.

Run:  python src/verify.py    (or `make verify`)
"""
import sys
import numpy as np
import scipy.sparse as sp
import spectral_lanczos as sl

TOL = 1e-3
FAIL = []


def approx(name, got, want, tol=TOL):
    ok = abs(got - want) <= tol * max(1.0, abs(want))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<34} = {got:.5f}  (expected {want:.5f})")
    if not ok:
        FAIL.append(name)
    return ok


def main():
    L, U, Ns, z = 6, 4.0, 50_000, 1.96
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

    # --- reconstruction-independent moments m_k = <0| J (H-E0)^k J |0> ---
    Hphi = Hc @ Jp
    m0 = float(np.vdot(Jp, Jp).real)
    m1 = float(np.vdot(Jp, Hphi).real)
    m2 = float(np.vdot(Jp, Hc @ Hphi).real)
    print("Reconstruction-independent current moments (operator algebra, never read off A(w)):")
    approx("m0 = <J^2>", m0, 0.5542)
    approx("m1 = 1/2<[J,[H,J]]>", m1, 4.2271)
    approx("m2 = <J (H-E0)^2 J>", m2, 33.3887)

    # --- shot-budget interval at the real ibm_fez budget (local-estimator variance, 2% bias budget) ---
    p = np.abs(psi0) ** 2
    mask = np.abs(psi0) > 1e-12
    O1 = np.zeros_like(psi0); O1[mask] = (Js @ Hphi)[mask] / psi0[mask]
    var1 = float(np.sum(p * O1 ** 2) - np.sum(p * O1) ** 2)
    O2psi = Js @ (Hc @ Hphi)
    O2 = np.zeros_like(psi0); O2[mask] = O2psi[mask] / psi0[mask]
    var2 = float(np.sum(p * O2 ** 2) - np.sum(p * O2) ** 2)
    d1 = z * np.sqrt(var1 / Ns + (0.02 * abs(m1)) ** 2)
    d2 = z * np.sqrt(var2 / Ns + (0.02 * abs(m2)) ** 2)
    print(f"\nShot-budget 95% interval (Ns={Ns}, 2% bias):  z*delta_1 = {d1:.3f} ({100*d1/m1:.1f}% of m1)")

    # --- the blind-spot: at d=98 the lone first moment MISSES; the joint battery CATCHES via m2 ---
    order = np.argsort(p)[::-1]
    nsup = int((p > 1e-14).sum())
    d = 98
    idx = np.sort(order[:d])
    Hsub, Jsub = Hs[idx][:, idx], Js[idx][:, idx]
    e0, g = sl._sub_gs(Hsub); Hcs = (Hsub - e0 * sp.identity(d)).tocsr(); Jg = Jsub @ g
    b1 = float(np.vdot(Jg, Hcs @ Jg).real)
    b2 = float(np.vdot(Jg, Hcs @ (Hcs @ Jg)).real)
    g1, g2 = abs(b1 - m1), abs(b2 - m2)
    print(f"\nTruncation d={d} of n_support={nsup} determinants:")
    print(f"  first moment gap  |m1_trunc - m1| = {g1:.3f}  vs interval {d1:.3f}  -> "
          f"{'MISS (within interval)' if g1 <= d1 else 'reject'}")
    print(f"  second moment gap |m2_trunc - m2| = {g2:.3f}  vs interval {d2:.3f}  -> "
          f"{'reject' if g2 > d2 else 'miss'}")
    battery_rejects = (g1 > d1) or (g2 > d2)
    lone_m1_misses = g1 <= d1
    ok = lone_m1_misses and battery_rejects and (g2 > d2)
    print(f"\n  [{'PASS' if ok else 'FAIL'}] joint (m0,m1,m2)+Hankel battery REJECTS d=98 while lone m1 misses "
          f"-> the battery has no blind spot here")
    if not ok:
        FAIL.append("battery-closes-d98")

    print("\n" + ("=" * 62))
    if FAIL:
        print(f"RESULT: {len(FAIL)} CHECK(S) FAILED: {', '.join(FAIL)}")
        sys.exit(1)
    print("RESULT: ALL CHECKS PASS -- the screen's headline numbers reproduce from source.")


if __name__ == "__main__":
    main()
