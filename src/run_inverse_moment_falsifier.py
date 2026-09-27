# -*- coding: utf-8 -*-
r"""
run_inverse_moment_falsifier.py  --  the NEGATIVE-ORDER (inverse) f-sum moment m_{-1} deployed as an
independent LOW-FREQUENCY falsifier, alongside the positive moments m_0,m_1,m_2 of the current
spectral function A_J(omega) = sum_n |<n|J|0>|^2 delta(omega - omega_n),  omega_n = E_n - E_0 > 0.

WHY m_{-1}.  The flagship error (Fig. fig_falsifier) is a SPURIOUS LOW-frequency (Drude-like) peak that
misplaces mid-infrared weight to small omega while preserving the total weight m_0. The positive power
moments m_k = int omega^k A_J domega weight HIGH frequency (kernel omega^k), so they are WEAKEST exactly
where this error lives; m_0 is blind by construction. The inverse moment
    m_{-1} = int A_J(omega)/omega domega = sum_n |<n|J|0>|^2 / omega_n
carries the OPPOSITE (1/omega) kernel: it is maximally sensitive to low-frequency weight and its baseline
is smallest (all true weight sits at high omega), so a low-omega error is amplified twice over.

THE EXACT OPERATOR IDENTITY (verified here to ~1e-14).  m_{-1} is the classic optical f-sum rule and is a
GROUND-STATE KINETIC ENERGY -- a ground-state expectation value, evaluated here by exact diagonalization (the
hopping operator is off-diagonal in the occupation basis, so on a device it would need rotated-basis
measurements, not the computational-basis samples alone; no moment is estimated on hardware in this work):

  * OPEN chain (polarization P = sum_j j n_j single-valued, current J = i[H,P] exact):
        m_{-1} = int A_J/omega domega = (1/2) <-That>  = (1/2) <[P,[H,P]]>   (EXACT, residual ~1e-15)
    with That = -t sum_<ij>,s (c^dag_{i,s} c_{j,s} + h.c.) the hopping (kinetic) operator.
    The optical-conductivity f-sum is the same statement with the Kubo prefactor:
        int_0^inf Re sigma_reg(omega) domega = pi * m_{-1} = (pi/2) <-That>   [Kohn1964, Maldague1977].

  * PERIODIC ring (P multivalued -> J != i[H,P] at the boundary bond): the identity acquires the exact
    charge-stiffness (Drude) correction Kohn established,
        m_{-1} = (1/2) <-That> - D ,     D = (1/2) d^2 E_0 / d theta^2   (Kohn flux-curvature stiffness),
    verified here against the finite-difference flux curvature. On the open chain D=0 and the clean
    kinetic-energy form is recovered. On ANY geometry m_{-1} = <0| J (H-E0)^{-1} Q J |0> (reduced
    resolvent, Q = 1-|0><0|) is exact and reconstruction-independent (verified ~1e-15).

We therefore DEPLOY the falsifier on the doped OPEN Hubbard chain, where m_{-1}=(1/2)<-That> is exact and
is the cheapest of all the moment estimators (one-body operator: JW weight-2 XX/YY on adjacent qubits;
two qubit-wise-commuting settings, a subset of what m_0/m_1 already require -- Sec. 'measurement cost').

Outputs (no matplotlib; native-pgfplots .dat + JSON only):
  data/<date>_inverse_moment_falsifier.json                (all verified numbers)
  paper/figs/inverse_falsifier_spectrum.dat                (A_true, A_wrong, 1/omega kernel; panel a)
  paper/figs/inverse_falsifier_sweep.dat                   (Delta_k/m_k vs misplaced fraction; panel b)
  paper/figs/inverse_falsifier_table.tex                   (\input-able booktabs-free tabular)
"""
import os, sys, json, datetime
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh, LinearOperator
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)   # hubbard_ed.py is vendored in src/
import spectral_lanczos as SL   # strings(), occ_matrix(), haydock_poles() reused verbatim

DATE = datetime.date.today().isoformat()
RES  = os.path.normpath(os.path.join(HERE, '..', 'data'))
# write native-pgfplots inputs to the paper's figs dir
FIGS_DIRS = [os.path.normpath(os.path.join(HERE, '..', 'paper', 'figs'))]
FIGS = FIGS_DIRS[0]

def _write_figs(name, text):
    for d in FIGS_DIRS:
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, name), 'w') as fh:
            fh.write(text)

# ----------------------------------------------------------------------------------------------------
# sector-basis operators with an OPEN/PERIODIC boundary flag (SL.hop/current_string are PBC-only)
# ----------------------------------------------------------------------------------------------------
def _bonds(L, pbc):
    return [(i, (i + 1) % L) for i in range(L if (pbc and L > 2) else L - 1)]

def hop_string(L, n, pbc, t=1.0):
    """single-spin hopping T_s = -t sum_bonds (c^dag_i c_j + h.c.) in the n-string basis."""
    S, idx = SL.strings(L, n); D = len(S); rows = []; cols = []; val = []
    for a, m in enumerate(S):
        for (i, j) in _bonds(L, pbc):
            for (p, q) in ((i, j), (j, i)):
                if (m >> q) & 1 and not (m >> p) & 1:
                    m2 = (m & ~(1 << q)) | (1 << p)
                    lo, hi = min(p, q), max(p, q)
                    mask = m & (((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1))
                    sign = -1.0 if (bin(mask).count('1') & 1) else 1.0
                    rows.append(idx[m2]); cols.append(a); val.append(-t * sign)
    return sp.csr_matrix((val, (rows, cols)), shape=(D, D))

def current_string(L, n, pbc, t=1.0):
    """single-spin current J_s = -i t sum_bonds (c^dag_{i+1} c_i - c^dag_i c_{i+1}) in the n-string basis."""
    S, idx = SL.strings(L, n); D = len(S); rows = []; cols = []; val = []
    for a, m in enumerate(S):
        for (i, j) in _bonds(L, pbc):
            for (p, q, coeff) in ((j, i, -1j * t), (i, j, +1j * t)):   # c^dag_p c_q
                if (m >> q) & 1 and not (m >> p) & 1:
                    m2 = (m & ~(1 << q)) | (1 << p)
                    lo, hi = min(p, q), max(p, q)
                    mask = m & (((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1))
                    sign = -1.0 if (bin(mask).count('1') & 1) else 1.0
                    rows.append(idx[m2]); cols.append(a); val.append(coeff * sign)
    return sp.csr_matrix((val, (rows, cols)), shape=(D, D))

def sector_matrices(L, U, nup, ndn, pbc, t=1.0):
    """dense sector H, hopping T, current J, polarization P for the (nup,ndn) product basis."""
    Tu = hop_string(L, nup, pbc, t); Td = hop_string(L, ndn, pbc, t)
    Su, _ = SL.strings(L, nup); Sd, _ = SL.strings(L, ndn)
    Du, Dd = len(Su), len(Sd)
    Iu, Id = sp.identity(Du), sp.identity(Dd)
    upocc = SL.occ_matrix(Su, L); dnocc = SL.occ_matrix(Sd, L)
    diagU = U * (upocc @ dnocc.T).ravel()
    Tfull = (sp.kron(Tu, Id) + sp.kron(Iu, Td))
    H = (Tfull + sp.diags(diagU)).tocsr()
    Ju = current_string(L, nup, pbc, t); Jd = current_string(L, ndn, pbc, t)
    J = (sp.kron(Ju, Id) + sp.kron(Iu, Jd)).tocsr()
    # polarization P = sum_j j (n^up_j + n^dn_j), diagonal in the product basis
    xup = (upocc @ np.arange(L));  xdn = (dnocc @ np.arange(L))
    Pdiag = (xup[:, None] + xdn[None, :]).ravel()
    P = sp.diags(Pdiag).tocsr()
    return H.toarray(), Tfull.toarray(), J.toarray(), P.toarray(), Du * Dd

# ----------------------------------------------------------------------------------------------------
# Part 1 -- EXACT identity verification (dense sector diagonalization)
# ----------------------------------------------------------------------------------------------------
def kohn_stiffness(L, U, nup, ndn, t=1.0, h=1e-4):
    """D = (1/2) d^2 E0/d theta^2 via Peierls per-bond flux theta (periodic ring), 2nd-order stencil."""
    def E0(theta):
        Su, iu = SL.strings(L, nup); Sd, idd = SL.strings(L, ndn)
        Du, Dd = len(Su), len(Sd)
        def hop_theta(n):
            S, idx = SL.strings(L, n); D = len(S); r = []; c = []; v = []
            for a, m in enumerate(S):
                for (i, j) in _bonds(L, True):
                    for (p, q, ph) in ((i, j, +1j * theta), (j, i, -1j * theta)):  # c^dag_p c_q, forward=+theta
                        if (m >> q) & 1 and not (m >> p) & 1:
                            m2 = (m & ~(1 << q)) | (1 << p)
                            lo, hi = min(p, q), max(p, q)
                            mask = m & (((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1))
                            sign = -1.0 if (bin(mask).count('1') & 1) else 1.0
                            r.append(idx[m2]); c.append(a); v.append(-t * np.exp(ph) * sign)
            return sp.csr_matrix((v, (r, c)), shape=(D, D))
        Tu = hop_theta(nup); Td = hop_theta(ndn)
        upocc = SL.occ_matrix(Su, L); dnocc = SL.occ_matrix(Sd, L)
        diagU = U * (upocc @ dnocc.T).ravel()
        H = (sp.kron(Tu, sp.identity(Dd)) + sp.kron(sp.identity(Du), Td) + sp.diags(diagU)).tocsr()
        return float(eigsh(H, k=1, which='SA')[0][0])
    return 0.5 * (E0(h) - 2 * E0(0.0) + E0(-h)) / h ** 2

def verify_identity(L, U, nup, ndn, pbc, t=1.0):
    H, T, J, P, dim = sector_matrices(L, U, nup, ndn, pbc, t)
    Ev, Vv = np.linalg.eigh(H)
    E0 = Ev[0]; psi0 = Vv[:, 0]
    negT = -float(np.real(psi0.conj() @ (T @ psi0)))
    Jpsi = J @ psi0
    ov = Vv.conj().T @ Jpsi; w = np.abs(ov) ** 2; om = Ev - E0
    reg = om > 1e-9; om_r, w_r = om[reg], w[reg]
    m0 = float(w_r.sum()); m1 = float((om_r * w_r).sum())
    m2 = float((om_r ** 2 * w_r).sum()); mm1 = float((w_r / om_r).sum())
    # reduced-resolvent form (exact on any geometry): m_{-1} = <Jpsi| (H-E0)^+ Q |Jpsi>
    Q = np.eye(dim) - np.outer(psi0, psi0.conj())
    x = np.linalg.pinv(H - E0 * np.eye(dim)) @ (Q @ Jpsi)
    mm1_res = float(np.real(np.vdot(Jpsi, x)))
    # polarization double commutator (1/2)<[P,[H,P]]>
    HP = H @ P - P @ H; PHP = P @ HP - HP @ P
    poldc = 0.5 * float(np.real(psi0.conj() @ (PHP @ psi0)))
    out = dict(L=L, U=U, nup=nup, ndn=ndn, pbc=pbc, dim=dim,
               half_negT=0.5 * negT, m0=m0, m1=m1, m2=m2, m_minus1=mm1,
               m_minus1_resolvent=mm1_res, half_pol_doublecomm=poldc,
               resid_resolvent=abs(mm1 - mm1_res), resid_poldc=abs(mm1 - poldc))
    if pbc:
        D = kohn_stiffness(L, U, nup, ndn, t)
        out['kohn_stiffness_D'] = D
        out['resid_decomp'] = abs((0.5 * negT) - (D + mm1))   # (1/2)<-T> = D + m_{-1}
    else:
        out['resid_identity_halfnegT'] = abs(mm1 - 0.5 * negT)  # m_{-1} = (1/2)<-T> EXACT
    return out

# ----------------------------------------------------------------------------------------------------
# Part 2 -- deployment: sensitivity of m_{-1} vs m_0,m_1,m_2 to the Drude/low-omega error (flagship scale)
# ----------------------------------------------------------------------------------------------------
def sector_H_linop(L, U, nup, ndn, pbc, t=1.0):
    Tu = hop_string(L, nup, pbc, t); Td = hop_string(L, ndn, pbc, t)
    Su, _ = SL.strings(L, nup); Sd, _ = SL.strings(L, ndn)
    Du, Dd = len(Su), len(Sd)
    upocc = SL.occ_matrix(Su, L); dnocc = SL.occ_matrix(Sd, L)
    diagU = U * (upocc @ dnocc.T).ravel()
    def matvec(x):
        X = x.reshape(Du, Dd)
        Y = Tu @ X + (Td @ X.T).T
        return Y.ravel() + diagU * x
    return LinearOperator((Du * Dd, Du * Dd), matvec=matvec, dtype=float), Du, Dd, Tu, Td, diagU

def _kohn_D_large(L, U, nup, ndn, t=1.0, h=1e-3):
    """Kohn charge stiffness D=(1/2)d^2E0/dtheta^2 at flagship scale via matrix-free sparse GS (Peierls flux)."""
    def E0(theta):
        def hop_theta(n):
            S, idx = SL.strings(L, n); D = len(S); r = []; c = []; v = []
            for a, m in enumerate(S):
                for (i, j) in _bonds(L, True):
                    for (p, q, ph) in ((i, j, +1j * theta), (j, i, -1j * theta)):
                        if (m >> q) & 1 and not (m >> p) & 1:
                            m2 = (m & ~(1 << q)) | (1 << p)
                            lo, hi = min(p, q), max(p, q)
                            mask = m & (((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1))
                            sign = -1.0 if (bin(mask).count('1') & 1) else 1.0
                            r.append(idx[m2]); c.append(a); v.append(-t * np.exp(ph) * sign)
            return sp.csr_matrix((v, (r, c)), shape=(D, D))
        Su, _ = SL.strings(L, nup); Sd, _ = SL.strings(L, ndn)
        Du, Dd = len(Su), len(Sd)
        Tu = hop_theta(nup); Td = hop_theta(ndn)
        upocc = SL.occ_matrix(Su, L); dnocc = SL.occ_matrix(Sd, L)
        diagU = U * (upocc @ dnocc.T).ravel()
        Iu, Id = sp.identity(Du), sp.identity(Dd)
        H = (sp.kron(Tu, Id) + sp.kron(Iu, Td) + sp.diags(diagU)).tocsr()
        return float(eigsh(H, k=1, which='SA', ncv=20, maxiter=5000)[0][0])
    return 0.5 * (E0(h) - 2 * E0(0.0) + E0(-h)) / h ** 2

def deploy(L=12, U=8.0, t=1.0, nl=260, om_split=6.0, om_spur=2.5, fmax=0.5, nf=11, eta=0.35):
    """Doped PERIODIC Hubbard ring (2/3 filling) -- the EXACT flagship system of fig_falsifier. On the ring
    the coherent low-omega current weight is the separate Drude stiffness D, so the REGULAR spectrum is
    depleted at low frequency and its inverse moment m_{-1} is small; a spurious low-omega peak therefore
    produces a LARGE relative discrepancy Delta_{-1}/m_{-1}, while the high-frequency-weighted positive
    moments respond only mildly. m_{-1} is the optical f-sum, tied to the ground state by the exact
    decomposition  m_{-1} = (1/2)<-That> - D  (open-chain / insulating limit: pure kinetic energy)."""
    nup = ndn = L // 3                                   # 2/3 filling (matches run_current / fig_falsifier)
    Hlin, Du, Dd, Tu, Td, diagU = sector_H_linop(L, U, nup, ndn, True, t)   # PERIODIC (flagship)
    # ground state (sparse Lanczos)
    Eg, Vg = eigsh(Hlin, k=1, which='SA', ncv=20, maxiter=5000)
    E0 = float(Eg[0]); psi0 = Vg[:, 0]
    # exact kinetic energy  <-That>  (one-body, sample-friendly) and Kohn stiffness D
    Psi = psi0.reshape(Du, Dd)
    TPsi = Tu @ Psi + (Td @ Psi.T).T
    negT = -float(np.real(np.vdot(psi0, TPsi.ravel())))
    D_kohn = _kohn_D_large(L, U, nup, ndn, t)
    m_minus1_fsum = 0.5 * negT - D_kohn                  # exact f-sum decomposition value on the ring
    # current Ritz poles
    Ju = current_string(L, nup, True, t); Jd = current_string(L, ndn, True, t)
    JPsi = (Ju @ Psi) + (Psi @ Jd.T)
    v0 = JPsi.ravel()
    om, w = SL.haydock_poles(lambda x: Hlin.matvec(x) - E0 * x, v0, nl)
    keep = (w > 1e-12) & (om > 1e-9)
    om, w = om[keep], w[keep]
    m0 = float(w.sum()); m1 = float((om * w).sum()); m2 = float((om ** 2 * w).sum())
    m_minus1 = float((w / om).sum())
    resid_kinetic = abs(m_minus1 - m_minus1_fsum) / abs(m_minus1)   # poles vs f-sum decomposition check

    hi = om >= om_split
    W_hi = float(w[hi].sum())

    def wrong(f):
        # steal fraction f of the mid-IR band into a spurious pole at om_spur (m_0 preserved exactly)
        oms = np.append(om, om_spur)
        wts = np.append(w.copy(), f * W_hi); wts[:-1][hi] *= (1.0 - f)
        return oms, wts

    # sensitivity sweep
    fs = np.linspace(0.0, fmax, nf)
    rows = []
    for f in fs:
        oms, wts = wrong(f)
        M0 = float(wts.sum()); M1 = float((oms * wts).sum())
        M2 = float((oms ** 2 * wts).sum()); Mm1 = float((wts / oms).sum())
        rows.append(dict(f=float(f),
                         d0=abs(M0 - m0) / m0,        # blind (=0 by construction)
                         d1=abs(M1 - m1) / m1,        # modest
                         d2=abs(M2 - m2) / m2,        # weakest at low freq (grows slowly)
                         dm1=abs(Mm1 - m_minus1) / m_minus1))   # LARGE (1/omega amplification)

    # ---- write native-pgfplots .dat (panel a: spectrum + 1/omega kernel; panel b: sensitivity) ----
    wgrid = np.linspace(0.0, 15.0, 600)
    def broaden(oms, wts):
        A = np.zeros_like(wgrid)
        for e, wt in zip(oms, wts):
            if wt > 1e-10:
                A += wt * (eta / np.pi) / ((wgrid - e) ** 2 + eta ** 2)
        return A
    A_true = broaden(om, w)
    ow, ww = wrong(0.40); A_wrong = broaden(ow, ww)
    ker = np.where(wgrid > 0.25, 1.0 / np.maximum(wgrid, 0.25), 1.0 / 0.25)  # 1/omega falsifier kernel
    kmax = float(A_true.max())
    spec = "omega Atrue Awrong kernel\n" + "".join(
        f"{x:.4f} {a:.6f} {b:.6f} {kmax*k/ker.max():.6f}\n"
        for x, a, b, k in zip(wgrid, A_true, A_wrong, ker))
    _write_figs('inverse_falsifier_spectrum.dat', spec)
    swp = "f dm1 d1 d2 d0\n" + "".join(
        f"{r['f']*100:.3f} {r['dm1']*100:.5f} {r['d1']*100:.5f} {r['d2']*100:.5f} {r['d0']*100:.6e}\n"
        for r in rows)
    _write_figs('inverse_falsifier_sweep.dat', swp)

    # ---- \input-able table (relative discrepancies at representative fractions) ----
    pick = [rows[i] for i in (2, 4, 6, 8, 10)]   # f = 10,20,30,40,50 %
    tab = ("% auto-generated by run_inverse_moment_falsifier.py -- do not edit by hand\n"
           "\\begin{tabular}{lccccc}\n\\hline\\hline\n"
           "misplaced fraction $f$ & 10\\% & 20\\% & 30\\% & 40\\% & 50\\% \\\\\n\\hline\n"
           "$\\Delta_0/m_0$ (total weight) & " + " & ".join("$<\\!10^{-13}$" for _ in pick) + " \\\\\n"
           "$\\Delta_1/m_1$ (first moment) & " + " & ".join(f"${r['d1']*100:.1f}\\%$" for r in pick) + " \\\\\n"
           "$\\Delta_2/m_2$ (second moment) & " + " & ".join(f"${r['d2']*100:.1f}\\%$" for r in pick) + " \\\\\n"
           "$\\Delta_{-1}/m_{-1}$ (\\textbf{inverse}) & "
           + " & ".join(f"$\\mathbf{{{r['dm1']*100:.0f}\\%}}$" for r in pick) + " \\\\\n"
           "\\hline\\hline\n\\end{tabular}\n")
    _write_figs('inverse_falsifier_table.tex', tab)

    return dict(L=L, U=U, filling='2/3', nup=nup, ndn=ndn, pbc=True, E0=E0, npoles=int(len(om)),
                negT=negT, half_negT=0.5 * negT, kohn_stiffness_D=D_kohn,
                m_minus1_fsum_decomp=m_minus1_fsum, m_minus1_poles=m_minus1,
                resid_fsum_vs_poles=resid_kinetic, m0=m0, m1=m1, m2=m2,
                W_low_below_split=float(w[~hi].sum()), frac_low=float(w[~hi].sum() / m0),
                om_split=om_split, om_spur=om_spur, W_midIR=W_hi, sweep=rows)

# ----------------------------------------------------------------------------------------------------
def main():
    os.makedirs(RES, exist_ok=True)
    print("=== Part 1: EXACT identity verification (dense sector ED) ===")
    checks = []
    for (L, U, nu, nd, pbc, tag) in [
        (6, 8.0, 3, 3, False, 'OBC half-filled (Mott)'),
        (6, 8.0, 2, 2, False, 'OBC doped 2/3'),
        (8, 8.0, 3, 3, False, 'OBC doped 3/8'),
        (6, 8.0, 3, 3, True,  'PBC half-filled (Mott)'),
        (6, 8.0, 2, 2, True,  'PBC doped 2/3'),
    ]:
        r = verify_identity(L, U, nu, nd, pbc); r['tag'] = tag; checks.append(r)
        if pbc:
            print(f"[{tag}] m_-1={r['m_minus1']:.10f}  (1/2)<-T>={r['half_negT']:.10f}  "
                  f"D_Kohn={r['kohn_stiffness_D']:.8f}  |(1/2)<-T>-(D+m_-1)|={r['resid_decomp']:.2e}  "
                  f"|m_-1-resolvent|={r['resid_resolvent']:.1e}")
        else:
            print(f"[{tag}] m_-1={r['m_minus1']:.10f}  (1/2)<-T>={r['half_negT']:.10f}  "
                  f"|m_-1-(1/2)<-T>|={r['resid_identity_halfnegT']:.2e}  "
                  f"|m_-1-(1/2)<[P,[H,P]]>|={r['resid_poldc']:.2e}  |m_-1-resolvent|={r['resid_resolvent']:.1e}")

    print("\n=== Part 2: deploy m_{-1} on the flagship doped PERIODIC ring (L=12, 2/3 filling) ===")
    dep = deploy(L=12, U=8.0, nl=260)
    print(f"L={dep['L']} U={dep['U']} filling={dep['filling']} ring  poles={dep['npoles']}  E0={dep['E0']:.6f}")
    print(f"<-That>={dep['negT']:.6f}  (1/2)<-That>={dep['half_negT']:.6f}  D_Kohn={dep['kohn_stiffness_D']:.6f}")
    print(f"m_-1=(1/2)<-That>-D={dep['m_minus1_fsum_decomp']:.6f}  m_-1(from poles)={dep['m_minus1_poles']:.6f}"
          f"  rel.resid={dep['resid_fsum_vs_poles']:.2e}")
    print(f"regular low-omega weight below split = {dep['frac_low']*100:.1f}% of m_0 (depleted -> m_-1 small)")
    print(f"m0=<J^2>={dep['m0']:.4f}  m1={dep['m1']:.4f}  m2={dep['m2']:.4f}  m_-1={dep['m_minus1_poles']:.4f}")
    print("sensitivity of each moment to the spurious low-omega (Drude) peak (relative discrepancy):")
    print(f"  {'f':>5} | {'D-1/m-1':>10} {'D1/m1':>8} {'D2/m2':>8} {'D0/m0':>10}")
    for r in dep['sweep'][::2]:
        print(f"  {r['f']*100:5.0f}% | {r['dm1']*100:9.1f}% {r['d1']*100:7.1f}% {r['d2']*100:7.1f}% {r['d0']*100:9.1e}%")

    verdict = (dep['sweep'][-1]['dm1'] > 2 * dep['sweep'][-1]['d1']) and (dep['sweep'][-1]['d0'] < 1e-9)
    out = {'_provenance': {'script': 'A_payload_echoes/03_src/run_inverse_moment_falsifier.py',
                           'date': DATE,
                           'claim': ('m_{-1}=int A_J/omega = (1/2)<-That> (optical f-sum, ground-state kinetic '
                                     'energy) on the open chain; = (1/2)<-That> - D_Kohn on the ring. Deployed '
                                     'as a LOW-frequency falsifier: Delta_{-1}/m_{-1} >> Delta_1/m_1 >> Delta_0/m_0=0 '
                                     'for the spurious-low-omega (Drude) error the positive moments are weakest on.')},
           'identity_checks': checks,
           'deployment': dep,
           'summary': {'inverse_moment_more_sensitive': bool(verdict),
                       'max_Delta_minus1_over_m': dep['sweep'][-1]['dm1'],
                       'max_Delta_1_over_m': dep['sweep'][-1]['d1'],
                       'sensitivity_ratio_minus1_to_1': dep['sweep'][-1]['dm1'] / max(dep['sweep'][-1]['d1'], 1e-12)}}
    jpath = os.path.join(RES, f'{DATE}_inverse_moment_falsifier.json')
    with open(jpath, 'w') as fh:
        json.dump(out, fh, indent=2)
    print(f"\nVERDICT: inverse moment is the more sensitive low-frequency falsifier = {verdict}")
    print("saved:", jpath)
    print("saved:", os.path.join(FIGS, 'inverse_falsifier_spectrum.dat'))
    print("saved:", os.path.join(FIGS, 'inverse_falsifier_sweep.dat'))
    print("saved:", os.path.join(FIGS, 'inverse_falsifier_table.tex'))


# ----------------------------------------------------------------------------------------------------
# R8 (added 2026-09-27): noise-normalized sensitivity of m_{-1} vs m_1, m_2, ring vs open chain.
# `python run_inverse_moment_falsifier.py --r8` writes ONLY the key 'R8_inverse_moment_sensitivity' of
# data/2026-09-27_theory_numerics.json; it writes no figure file and no *_inverse_moment_falsifier.json.
# main() (no flag) is unchanged.
# ----------------------------------------------------------------------------------------------------
R8_NS, R8_Z, R8_BIAS = 50000, 1.96, 0.02


def _r8_measure(L, U, t, nl, pbc, seed):
    """current Ritz poles, <-T>, and local-estimator variances of the moment operators on the sector GS."""
    nup = ndn = L // 3
    Hlin, Du, Dd, Tu, Td, diagU = sector_H_linop(L, U, nup, ndn, pbc, t)
    v0 = np.random.default_rng(seed).standard_normal(Du * Dd)
    Eg, Vg = eigsh(Hlin, k=2, which='SA', ncv=30, maxiter=5000, v0=v0, tol=1e-12)
    o = np.argsort(Eg); E0 = float(Eg[o[0]]); psi0 = Vg[:, o[0]]; gap = float(Eg[o[1]] - Eg[o[0]])
    Psi = psi0.reshape(Du, Dd)
    TPsi = (Tu @ Psi + (Td @ Psi.T).T).ravel()
    negT = -float(np.real(np.vdot(psi0, TPsi)))
    Ju = current_string(L, nup, pbc, t); Jd = current_string(L, ndn, pbc, t)
    Jop = lambda x: ((Ju @ x.reshape(Du, Dd)) + (x.reshape(Du, Dd) @ Jd.T)).ravel()
    Hm = lambda x: Hlin.matvec(x) - E0 * x
    v1 = Jop(psi0.astype(complex))
    om, w = SL.haydock_poles(Hm, v1, nl)
    keep = (w > 1e-12) & (om > 1e-9); om, w = om[keep], w[keep]
    # local estimators over |psi0(x)|^2 (shot term of the deployed rule; a LOWER bound on per-shot variance)
    p = np.abs(psi0) ** 2; mask = np.abs(psi0) > 1e-12

    def locvar(Opsi):
        loc = np.zeros(len(psi0), dtype=complex); loc[mask] = Opsi[mask] / psi0[mask]
        mean = float(np.real(np.sum(p * loc))); return mean, float(np.sum(p * np.abs(loc) ** 2) - mean ** 2)
    mhalfT, var_halfT = locvar(-0.5 * TPsi)                 # (1/2)<-T>
    Hv1 = Hm(v1); H2v1 = Hm(Hv1)
    m0l, var0 = locvar(Jop(v1)); m1l, var1 = locvar(Jop(Hv1)); m2l, var2 = locvar(Jop(H2v1))
    return dict(E0=E0, gap_to_next=gap, negT=negT, half_negT=0.5 * negT, om=om, w=w,
                local=dict(half_negT=(mhalfT, var_halfT), m0=(m0l, var0), m1=(m1l, var1), m2=(m2l, var2)))


def _r8_analysis(meas, D_kohn, om_split=6.0, om_spur=2.5, fs=(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)):
    om, w = meas['om'], meas['w']
    m0 = float(w.sum()); m1 = float((om * w).sum()); m2 = float((om ** 2 * w).sum())
    mm1 = float((w / om).sum())
    hi = om >= om_split; W_hi = float(w[hi].sum())
    # exact linear shifts per unit f (fraction f of the weight at om>=om_split moved to om_spur; m0 exact)
    dM = {k: W_hi * om_spur ** k - float((w[hi] * om[hi] ** k).sum()) for k in (-1, 1, 2)}
    mk = {-1: mm1, 1: m1, 2: m2}
    slope = {k: abs(dM[k]) / mk[k] for k in dM}
    est = {-1: meas['half_negT'], 1: m1, 2: m2}             # the estimated operator expectation (bias base)
    var = {-1: meas['local']['half_negT'][1], 1: meas['local']['m1'][1], 2: meas['local']['m2'][1]}
    tau_b = {k: R8_Z * R8_BIAS * abs(est[k]) for k in dM}
    tau_f = {k: R8_Z * np.sqrt(var[k] / R8_NS + (R8_BIAS * est[k]) ** 2) for k in dM}
    fstar = lambda tau: {('m_-1' if k == -1 else f'm{k}'): (float(tau[k] / abs(dM[k])) if abs(dM[k]) > 0 else None)
                         for k in dM}
    name = lambda k: 'm_-1' if k == -1 else f'm{k}'
    return dict(
        m0=m0, m1=m1, m2=m2, m_minus1_poles=mm1, half_negT=meas['half_negT'], D_kohn=D_kohn,
        m_minus1_fsum=meas['half_negT'] - D_kohn, resid_poles_vs_fsum_rel=abs(mm1 - (meas['half_negT'] - D_kohn)) / mm1,
        cancellation_D_over_half_negT=D_kohn / meas['half_negT'],
        E0=meas['E0'], gap_to_next=meas['gap_to_next'], npoles=int(len(om)),
        W_above_split=W_hi, frac_below_split=float(w[~hi].sum() / m0),
        shift_per_unit_f_abs={name(k): dM[k] for k in dM},
        slope_relative_per_unit_f={name(k): slope[k] for k in slope},
        slope_ratio_minus1_over_1=slope[-1] / slope[1], slope_ratio_minus1_over_2=slope[-1] / slope[2],
        rel_shift_at_f=[dict(f=f, dm1=slope[-1] * f, d1=slope[1] * f, d2=slope[2] * f) for f in fs],
        tau_bias_only={name(k): tau_b[k] for k in tau_b},
        tau_bias_only_over_m={name(k): tau_b[k] / mk[k] for k in tau_b},
        tau_full_rule={name(k): tau_f[k] for k in tau_f},
        tau_full_rule_over_m={name(k): tau_f[k] / mk[k] for k in tau_f},
        local_variance={name(k): var[k] for k in var},
        f_star_bias_only=fstar(tau_b), f_star_full_rule=fstar(tau_f))


def r8_main(L=12, U=8.0, t=1.0, nl=260, seed=0):
    import time as _time
    from small_checks import record_theory_numerics
    t0 = _time.time()
    ring_meas = _r8_measure(L, U, t, nl, True, seed)
    D = _kohn_D_large(L, U, L // 3, L // 3, t)
    ring = _r8_analysis(ring_meas, D)
    obc_meas = _r8_measure(L, U, t, nl, False, seed)
    obc = _r8_analysis(obc_meas, 0.0)
    # secondary: does the conclusion depend on where the spurious pole sits? (same measures, exact arithmetic)
    for tag, meas, Dk, rec in (('ring', ring_meas, D, ring), ('open chain', obc_meas, 0.0, obc)):
        scan = []
        for osp in (0.25, 0.5, 1.0, 1.5, 2.5, 4.0):
            r = _r8_analysis(meas, Dk, om_spur=osp)
            scan.append(dict(om_spur=osp, f_star_bias_only=r['f_star_bias_only'],
                             f_star_full_rule=r['f_star_full_rule'],
                             slope_relative_per_unit_f=r['slope_relative_per_unit_f']))
        rec['om_spur_scan'] = scan
        print(f"[{tag}] om_spur scan (f* bias-only m_-1 / m1 / m2): " + "; ".join(
            f"{s['om_spur']}: {s['f_star_bias_only']['m_-1']:.3f}/{s['f_star_bias_only']['m1']:.3f}/"
            f"{s['f_star_bias_only']['m2']:.3f}" for s in scan))
    for tag, r in (('ring', ring), ('open chain', obc)):
        print(f"[{tag}] m_-1(poles)={r['m_minus1_poles']:.5f} (1/2)<-T>={r['half_negT']:.5f} D={r['D_kohn']:.5f} "
              f"m1={r['m1']:.4f} m2={r['m2']:.3f}; slopes rel/f: m_-1 {r['slope_relative_per_unit_f']['m_-1']:.3f} "
              f"m1 {r['slope_relative_per_unit_f']['m1']:.3f} m2 {r['slope_relative_per_unit_f']['m2']:.3f}")
        print(f"      tau/m (bias only) m_-1 {r['tau_bias_only_over_m']['m_-1']:.3f} m1 {r['tau_bias_only_over_m']['m1']:.4f};"
              f" f* bias-only {r['f_star_bias_only']}  f* full rule {r['f_star_full_rule']}")
    rt = _time.time() - t0
    R8 = dict(
        system=f'L={L}, U/t={U}, 2/3 filling (N_up=N_dn={L//3}); ring = deployment of this script; open chain = '
               'same sector with open boundary (D = 0)',
        corruption='fraction f of the regular weight at omega >= 6 t moved into a spurious pole at 2.5 t '
                   '(m0 exact) -- the deploy() error model; shifts are exactly linear in f',
        estimator_for_m_minus1='ring: (1/2)<-T> - D_Kohn (D computed separately by flux curvature, not same-sample); '
                               'open chain: (1/2)<-T> exactly',
        threshold_rules={'bias_only': 'tau_k = 1.96 * 0.02 * |estimated expectation| ((1/2)<-T> for m_-1)',
                         'full_rule': 'tau_k = 1.96 sqrt(Var(O_k^loc)/5e4 + (0.02 x)^2), x the estimated expectation'},
        f_star_def='f at which |Delta_k| = tau_k (fires above); None if the moment is blind',
        ring=ring, open_chain=obc,
        committed_ring_values={'half_negT': 5.024889909209261, 'kohn_stiffness_D': 4.700763076570524,
                               'm_minus1': 0.32412604304278386, 'm1': 35.14516230804411,
                               'source': 'data/2026-08-28_inverse_moment_falsifier.json deployment'},
        plan_expectation='tau_-1 = 0.197 (61% of m_-1 = 0.324); f* ~ 0.20 for m_-1 vs ~0.05 for m1; ring = 94% '
                         'cancellation (5.025 - 4.701); advantage expected to hold cleanly on open chains')
    record_theory_numerics('R8_inverse_moment_sensitivity', R8, 'run_inverse_moment_falsifier.py', rt, seed=seed,
                           command='cd src && python run_inverse_moment_falsifier.py --r8')


if __name__ == '__main__':
    if '--r8' in sys.argv:
        r8_main()
    else:
        main()
