# -*- coding: utf-8 -*-
"""Export the real (m1,m2) point for NEW fig_momentcone. Measure = doped Hubbard L=6 U=8 current-operator
Lehmann density, normalized (m0=1) and rescaled to [-1,1]. A valid moment sequence in this rescaled frame
must satisfy the Hankel-positivity geometry m1^2 <= m2 <= 1; the real point lies inside. (The m1>=0
shifted/Stieltjes cut is a physical [0,inf)-frame condition and does NOT apply after rescaling -- the real
point has m1<0 -- so it is intentionally not drawn.)"""
import os, sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)   # hubbard_ed.py is vendored in src/
import hubbard_ed as H
from run_sumrule_falsifier import current_operator
OUT = os.path.normpath(os.path.join(HERE, '..', 'paper', 'figs'))
L, U = 6, 8.0


def main():
    c = H.build_operators(L); cd = [ci.getH() for ci in c]
    Ham = H.hubbard_hamiltonian(L, t=1.0, U=U, pbc=True, c=c)
    Ntot = sum((cd[m] @ c[m] for m in range(2 * L)), sp.csr_matrix((2 ** (2 * L),) * 2, dtype=complex))
    J = current_operator(L, c)
    for mu_ in [0.0, -1.0, -2.0]:
        cand = eigsh(Ham - mu_ * Ntot, k=1, which='SA')[1][:, 0]
        if float(np.real(np.vdot(cand, Ntot @ cand))) <= L - 1.5:
            psi0 = cand; break
    E0 = float(np.real(np.vdot(psi0, Ham @ psi0)))
    Jpsi = J @ psi0
    Ev, Vv = np.linalg.eigh(Ham.toarray())
    om = Ev - E0; w = np.abs(Vv.conj().T @ Jpsi) ** 2
    keep = om > 1e-6; om, w = om[keep], w[keep]
    wn = w / w.sum()
    a, b = om.min(), om.max()
    xs = (2 * om - (a + b)) / (b - a)
    m1 = float(np.sum(wn * xs)); m2 = float(np.sum(wn * xs ** 2))
    with open(os.path.join(OUT, 'momentcone_point.dat'), 'w') as f:
        f.write("m1 m2\n%.5f %.5f\n" % (m1, m2))
    print("normalized measure on [-1,1]: m1=%.4f  m2=%.4f  (m1^2=%.4f, must have m2>=m1^2: %s)"
          % (m1, m2, m1 ** 2, m2 >= m1 ** 2))
    print("wrote momentcone_point.dat ->", OUT)


if __name__ == '__main__':
    main()
