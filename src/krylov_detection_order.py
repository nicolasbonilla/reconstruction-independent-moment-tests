# -*- coding: utf-8 -*-
"""Post hoc (NOT sealed) analysis of 2026-11-08: from which moment order does an
under-converged Krylov (n-node Lanczos) reconstruction become detectable, and at what shot cost?

System: the blinded catalog's calibration system -- current probe J of the doped L=6 Hubbard ring,
(N_up, N_dn) = (2, 2), sector dimension 225, U/t in {3, 4, 5, 6} (blind_generate.build_U).
Reconstruction: the n-node Lanczos/Haydock pole representation of J|0> (spectral_lanczos.haydock_poles),
n = 1..4, exactly as the sealed catalog builds its 'krylov' instances.
Oracle estimator: exact full-operator moments m_k = <J0|(H-E0)^k|J0> (no noise); thresholds of the
default rule extended to k <= 2*nmax + 2 = 10: tau_k = z sqrt(Var(O_k)/N_s + (b m_k)^2), z = 1.96, N_s = 5e4,
b = 0.02, with Var(O_k) the exact-state local-estimator variance (within_sector_lp.deployed_taus).

Checks the classical Gauss-quadrature remainder: for an n-node Lanczos rule started from phi = J|0>,
m_k - mbar_k = 0 for k <= 2n-1 and m_{2n} - mbar_{2n} = ||phi||^2 prod_{j=1}^{n} beta_j^2 > 0.

Reuses the repository code (src/) read-only. Runs in seconds on one process.
"""
import os, sys, json, time, platform
import numpy as np
import scipy, scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
# Works in two places: inside the repository's src/ (next to spectral_lanczos.py, output to ../data) or outside
# it, where it imports the repository's src/ from the path held in SRC (read-only) and writes next to this script.
IN_REPO = os.path.exists(os.path.join(HERE, 'spectral_lanczos.py'))
SRC = HERE if IN_REPO else os.path.normpath(os.path.join(HERE, '..', '..', '..', 'release', 'src'))
OUTDIR = os.path.normpath(os.path.join(HERE, '..', 'data')) if IN_REPO else HERE
sys.path.insert(0, SRC)
import spectral_lanczos as sl            # noqa: E402  (repository code, unchanged)
from blind_generate import build_U       # noqa: E402  (same system builder as the sealed catalog)

NS, Z, BIAS = 50000, 1.96, 0.02
Z_POWER80 = 1.96 + 0.8416                # z_{1-alpha/2} + z_{1-beta}, power 0.80 at alpha = 0.05
UGRID = [3.0, 4.0, 5.0, 6.0]
NMAX = 4
KMAX = 2 * NMAX + 2                      # moments 0..10


ETA = 0.15
GRID_WIDE = np.linspace(-40.0, 80.0, 120001)     # as full_support_rr.py (d omega = 1e-3 t)
_trapz = getattr(np, 'trapezoid', None) or np.trapz


def rel_l1(om_a, w_a, om_b, w_b, grid=GRID_WIDE):
    """int |A_b - A_a| / int A_a, Lorentzian eta = 0.15 t (definition of full_support_rr.rel_l1; A_a exact)."""
    def br(om, w):
        return (w[None, :] * (ETA / np.pi) / ((grid[:, None] - om[None, :]) ** 2 + ETA ** 2)).sum(1)
    A = br(om_a, w_a); B = br(om_b, w_b)
    return float(_trapz(np.abs(B - A), grid) / _trapz(A, grid))


def lanczos_betas(Hc, phi, n):
    """Lanczos coefficients (alpha_0..alpha_{n-1}, beta_1..beta_n) from phi, same recurrence as
    spectral_lanczos.haydock_poles; beta_n is the residual norm after n steps."""
    v = phi / np.linalg.norm(phi); vp = np.zeros_like(v); b = 0.0; a, bs = [], []
    for _ in range(n):
        w = Hc @ v; an = float(np.vdot(v, w).real); a.append(an)
        w = w - an * v - b * vp
        bn = float(np.linalg.norm(w)); bs.append(bn)
        if bn < 1e-12:
            break
        vp = v; v = w / bn; b = bn
    return np.array(a), np.array(bs)


def exact_moments_and_var(psi0, Js, Hc, kmax):
    """m_k and the local-estimator variance of O_k = J (H-E0)^k J over |psi0(x)|^2 (as deployed_taus)."""
    p = np.abs(psi0) ** 2; mask = np.abs(psi0) > 1e-12
    phi = Js @ psi0; v = phi.copy(); m, var = [], []
    for k in range(kmax + 1):
        Okpsi = Js @ v
        loc = np.zeros(len(psi0), dtype=complex); loc[mask] = Okpsi[mask] / psi0[mask]
        mean = float(np.sum(p * loc).real)
        var.append(max(float(np.sum(p * np.abs(loc) ** 2) - mean ** 2), 0.0))
        m.append(float(np.vdot(phi, v).real))
        v = Hc @ v
    return np.array(m), np.array(var), phi


def main():
    t0 = time.time()
    rows = []
    for U in UGRID:
        Hs, Js, Hc, psi0 = build_U(6, U)
        m, var, phi = exact_moments_and_var(psi0, Js, Hc, KMAX)
        Ev = np.linalg.eigh(Hs.toarray())[0]
        tau = Z * np.sqrt(var / NS + (BIAS * np.abs(m)) ** 2)
        tau_shot = Z * np.sqrt(var / NS)
        # number of distinct weighted poles, from the full eigendecomposition (as within_sector_lp)
        Ev_, Vv_ = np.linalg.eigh(Hs.toarray())
        wts = np.abs(Vv_.T @ phi) ** 2; om_ = Ev_ - Ev_[0]
        keep = wts > 1e-10
        n_poles = int(len(np.unique(np.round(om_[keep], 6))))
        om_ex, w_ex = om_[keep], wts[keep]
        for n in range(1, NMAX + 1):
            theta, wj = sl.haydock_poles(lambda x: Hc @ x, phi, n)
            mbar = np.array([float(np.sum(wj * theta ** k)) for k in range(KMAX + 1)])
            dlt = m - mbar
            a, b = lanczos_betas(Hc, phi, n)
            gauss_rem = float(np.vdot(phi, phi).real * np.prod(b[:n] ** 2))
            k_star = 2 * n
            # Chebyshev-Jensen bound (Corollary 2 of the manuscript): r_2n <= 4 [D/(4 wbar)]^(2n), with D the
            # diameter of the support of the spectral measure of phi and wbar = m1/m0 its mean.
            D_supp = float(om_ex.max() - om_ex.min()); wbar = float(m[1] / m[0])
            bound_rel = 4.0 * (D_supp / (4.0 * wbar)) ** (2 * n)
            bound_abs = 4.0 * float(m[0]) * (D_supp / 4.0) ** (2 * n)
            first_exposed = next((k for k in range(KMAX + 1) if abs(dlt[k]) > tau[k]), None)
            first_exposed_shot = next((k for k in range(KMAX + 1) if abs(dlt[k]) > tau_shot[k]), None)
            rel = abs(dlt[k_star]) / abs(m[k_star])
            # shots for exposure (noise-free shift above threshold) and for 80% power, no bias / 2% bias
            def shots(zz, bias):
                den = dlt[k_star] ** 2 - (zz * bias * m[k_star]) ** 2
                return float(zz ** 2 * var[k_star] / den) if den > 0 else None
            rows.append(dict(
                U=U, n_nodes=n, n_weighted_poles=n_poles,
                rel_L1_eta0p15_wide_grid=rel_l1(om_ex, w_ex, theta, wj),
                ritz_nodes=[float(x) for x in theta], ritz_weights=[float(x) for x in wj],
                weighted_poles=[float(x) for x in np.unique(np.round(om_ex, 6))],
                support_diameter=D_supp, mean_frequency=wbar, D_over_4wbar=D_supp / (4.0 * wbar),
                cheb_jensen_bound_rel=bound_rel, cheb_bound_abs=bound_abs,
                bound_holds=bool(abs(dlt[k_star]) <= bound_abs * (1 + 1e-12) and
                                 abs(dlt[k_star]) / abs(m[k_star]) <= bound_rel * (1 + 1e-12)),
                bound_over_actual=float(bound_rel / (abs(dlt[k_star]) / abs(m[k_star]))),
                Ns_lower_bound_from_cheb=float(Z ** 2 * var[k_star] / bound_abs ** 2),
                bias_floor_excludes_order_2n=bool(bound_rel <= Z * BIAS),
                max_rel_residual_below_2n=float(max(abs(dlt[k]) / abs(m[k]) for k in range(k_star))),
                k_first_nonzero=k_star,
                delta_2n=float(dlt[k_star]), gauss_remainder=gauss_rem,
                gauss_identity_relerr=float(abs(dlt[k_star] - gauss_rem) / gauss_rem),
                m_2n=float(m[k_star]), rel_residual_2n=float(rel),
                tau_2n_default=float(tau[k_star]), ratio_default=float(abs(dlt[k_star]) / tau[k_star]),
                tau_2n_shot=float(tau_shot[k_star]), ratio_shot=float(abs(dlt[k_star]) / tau_shot[k_star]),
                rel_residual_2n_over_bias_floor=float(rel / (Z * BIAS)),
                Ns_expose_nobias=shots(Z, 0.0), Ns_expose_bias2pct=shots(Z, BIAS),
                Ns_power80_nobias=shots(Z_POWER80, 0.0), Ns_power80_bias2pct=shots(Z_POWER80, BIAS),
                bias_frac_max_for_exposure=float(rel / Z),
                first_order_exposed_default=first_exposed, first_order_exposed_shot_only=first_exposed_shot,
                ratio_default_by_k=[float(abs(dlt[k]) / tau[k]) for k in range(KMAX + 1)],
                rel_residual_by_k=[float(abs(dlt[k]) / abs(m[k])) for k in range(KMAX + 1)],
                tau_rel_default_by_k=[float(tau[k] / abs(m[k])) for k in range(KMAX + 1)],
                tau_rel_shot_by_k=[float(tau_shot[k] / abs(m[k])) for k in range(KMAX + 1)],
                var_local_by_k=[float(v) for v in var], m_exact_by_k=[float(x) for x in m],
            ))
        rows_U = [r for r in rows if r['U'] == U]
        print(f"U={U}: poles={n_poles}  m0={m[0]:.4f}  E0={Ev[0]:.6f}")
        for r in rows_U:
            print(f"  n={r['n_nodes']}: below-2n max rel {r['max_rel_residual_below_2n']:.1e} | "
                  f"rel D_2n={r['rel_residual_2n']:.3e} (floor zb={Z*BIAS:.4f}) | D/tau_def={r['ratio_default']:.3f} "
                  f"D/tau_shot={r['ratio_shot']:.2f} | Ns_exp(b=0)={r['Ns_expose_nobias']} | "
                  f"first exposed (default, shot-only) = {r['first_order_exposed_default']}, "
                  f"{r['first_order_exposed_shot_only']} | gauss relerr {r['gauss_identity_relerr']:.1e}")
    out = dict(
        _what="Post hoc, unsealed: detection order and shot cost of n-node Lanczos reconstructions, "
              "current probe of the doped L=6 ring (blinded catalog system), oracle (exact) moments.",
        _date="2026-11-08",
        _script="src/krylov_detection_order.py" if IN_REPO else "work/PRB/sim/krylov_detection_order.py",
        _uses=["release/src/spectral_lanczos.py (haydock_poles)", "release/src/blind_generate.py (build_U)"],
        rule=dict(Ns=NS, z=Z, bias_frac=BIAS, z_power80=Z_POWER80, variance="exact-state local estimator"),
        python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
        runtime_s=round(time.time() - t0, 2), rows=rows)
    path = os.path.join(OUTDIR, '2026-11-08_krylov_detection_order.json')
    json.dump(out, open(path, 'w'), indent=1)
    print(f"wrote {path}  ({out['runtime_s']} s)")


if __name__ == '__main__':
    main()
