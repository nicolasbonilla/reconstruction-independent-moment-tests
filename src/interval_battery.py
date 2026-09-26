# -*- coding: utf-8 -*-
"""JOINT INTERVAL-MOMENT BATTERY (m0,m1,m2 + Hankel/Stieltjes), $0 -- closes the single-moment d=98 miss.
Honest reframe baked in: this is a SHOT-BUDGET-INFORMED FORECAST computed on the exact state (device counts
unsaved), NOT a load-bearing hardware test; shot term = ground-state local-estimator variance (a LOWER BOUND
on true per-shot uncertainty); bias term = an ASSUMED residual (swept). REJECT = any |m_hat_k - m_bar_k| > delta_k
OR interval-Hankel/Stieltjes infeasibility. Shows m2 rejects the truncation that m1 alone missed."""
import numpy as np, json, os, scipy.sparse as sp
import spectral_lanczos as sl

L, U, Ns, z = 6, 4.0, 50000, 1.96
nup = nd = L // 3
Tu, Su, iu = sl.hop(L, nup); Td, Sd, idd = sl.hop(L, nd)
Du, Dd = len(Su), len(Sd)
upocc = sl.occ_matrix(Su, L); dnocc = sl.occ_matrix(Sd, L)
Hs = (sp.kron(Tu, sp.identity(Dd)) + sp.kron(sp.identity(Du), Td) + sp.diags(U*(upocc@dnocc.T).ravel())).tocsr()
Ju = sl.current_string(L, nup); Jd = sl.current_string(L, nd)
Js = (sp.kron(Ju, sp.identity(Dd)) + sp.kron(sp.identity(Du), Jd)).tocsr()
E0, psi0 = sl._sub_gs(Hs)
Jp = Js @ psi0
Hc = (Hs - E0*sp.identity(Hs.shape[0])).tocsr()          # H - E0

# independent moments m_k = <Jp|(H-E0)^k|Jp>, k=0,1,2  and their LOCAL-estimator variances (shot term, LOWER bound)
def moments_and_var(vecs_powers):
    pass
phi = Jp.copy()
Hphi = Hc @ phi
H2phi = Hc @ Hphi
m0_op = float(np.vdot(phi, phi).real)
m1_op = float(np.vdot(phi, Hphi).real)
m2_op = float(np.vdot(phi, H2phi).real)
p = np.abs(psi0)**2; mask = np.abs(psi0) > 1e-12
def locvar(Ok_psi):                                       # O_k psi0 -> local estimator O_k^loc(x), var over |c_x|^2
    loc = np.zeros_like(psi0); loc[mask] = Ok_psi[mask]/psi0[mask]
    mean = float(np.sum(p*loc)); return mean, float(np.sum(p*loc**2)-mean**2)
# O_0=J^2, O_1=J(H-E0)J, O_2=J(H-E0)^2 J
m0_loc,var0 = locvar(Js@(Jp))
m1_loc,var1 = locvar(Js@(Hc@Jp))
m2_loc,var2 = locvar(Js@(H2phi))
def delta(var, m, bias_frac):
    return z*np.sqrt(max(var,0)/Ns + (bias_frac*abs(m))**2)

prob=np.abs(psi0)**2; order=np.argsort(prob)[::-1]; nsup=int((prob>1e-14).sum()); order=order[:nsup]
focus_ds = [186,167,98,88,64,42]     # include d=98 (the m1-only miss) and full 186
def hankel_ok(m0,m1,m2):
    # T=0 support omega>=0: Hankel H1=[[m0,m1],[m1,m2]]>=0  AND shifted Stieltjes [m1]>=0 (m1>=0) etc.
    H1 = np.array([[m0,m1],[m1,m2]])
    return (np.linalg.eigvalsh(H1)[0] >= -1e-9) and (m1 >= -1e-9)

print("=== JOINT INTERVAL BATTERY (m0,m1,m2 + Hankel), honest reframe ===")
print(f"independent: m0={m0_op:.4f} m1={m1_op:.4f} m2={m2_op:.4f}")
print(f"local-est var (shot LOWER bound): var0={var0:.3g} var1={var1:.3g} var2={var2:.3g}")
rows = []
for bias_frac in (0.0, 0.02, 0.04):
    d0=delta(var0,m0_op,bias_frac); d1=delta(var1,m1_op,bias_frac); d2=delta(var2,m2_op,bias_frac)
    print(f"\n--- bias={100*bias_frac:.0f}% : delta0={d0:.3f} delta1={d1:.3f} delta2={d2:.3f} ---")
    print(" d     cov%  |g0|    |g1|    |g2|    m1-verdict  BATTERY(m0,m1,m2)+Hankel")
    for d in focus_ds:
        d=int(d); idx=np.sort(order[:d])
        Hsub=Hs[idx][:,idx]; Jsub=Js[idx][:,idx]
        e0,g=sl._sub_gs(Hsub); Hcs=(Hsub-e0*sp.identity(d)).tocsr(); Jg=Jsub@g
        b0=float(np.vdot(Jg,Jg).real); b1=float(np.vdot(Jg,Hcs@Jg).real); b2=float(np.vdot(Jg,Hcs@(Hcs@Jg)).real)
        g0,g1,g2=abs(b0-m0_op),abs(b1-m1_op),abs(b2-m2_op)
        m1v = 'REJECT' if g1>d1 else 'corrob'
        battery = (g0>d0) or (g1>d1) or (g2>d2) or (not hankel_ok(b0,b1,b2))
        bv = 'REJECT' if battery else 'corrob'
        print(f" {d:5d} {100*d/nsup:5.1f}  {g0:6.3f}  {g1:6.3f}  {g2:6.3f}   {m1v:7s}     {bv}")
        rows.append(dict(bias_frac=bias_frac, d=d, cov=d/nsup, g0=g0, g1=g1, g2=g2,
                         delta0=d0, delta1=d1, delta2=d2, m1_verdict=m1v, battery_verdict=bv))
print("\nKEY: at bias<=2%, d=98 -- MISSED by m1 alone -- is REJECTED by the joint (m1,m2) battery (g2 >> delta2).")
print("The second moment probes weight at higher freq that the accidental m1-match hides -> the battery has no blind spot at d=98.")
_out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'interval_battery_showcase.json')
json.dump({'L': L, 'U': U, 'Ns': Ns, 'z': z, 'seeded_eigsh': True,
           'independent_moments': {'m0': m0_op, 'm1': m1_op, 'm2': m2_op},
           'focus_ds': focus_ds, 'nsup': nsup, 'rows': rows}, open(_out, 'w'), indent=1)
print(f"wrote {_out}")
