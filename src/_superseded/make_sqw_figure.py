# -*- coding: utf-8 -*-
"""make_sqw_figure.py — the SECOND rich showpiece: the collective dynamical structure factors of the Hubbard
ring, S(q,omega) [charge, Mott-gapped] and S^zz(q,omega) [spin, gapless spinon continuum] — spin-charge
separation, echoing the companion paper. Beautiful warm heatmaps in the paper style; the screen operates on both."""
import os, sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)   # quarantined in src/_superseded/; hubbard_ed.py etc. live in src/
sys.path.insert(0, SRC); sys.path.insert(0, HERE)
import hubbard_ed as H
import paper_style as PS; PS.apply()
PAP = os.path.normpath(os.path.join(HERE, '..', 'paper', 'figs')); os.makedirs(PAP, exist_ok=True)

L, U = 6, 8.0
NW = 380


def structure_factors():
    c = H.build_operators(L); cd = [ci.getH() for ci in c]
    Ham = H.hubbard_hamiltonian(L, t=1.0, U=U, pbc=True, c=c)
    Nop = sum((cd[m] @ c[m] for m in range(2 * L)), sp.csr_matrix((2 ** (2 * L),) * 2, dtype=complex))
    Sz = sum((0.5 * (cd[2 * s] @ c[2 * s] - cd[2 * s + 1] @ c[2 * s + 1]) for s in range(L)),
             sp.csr_matrix((2 ** (2 * L),) * 2, dtype=complex))
    I = sp.identity(2 ** (2 * L), dtype=complex, format='csr')
    Hpen = Ham + 30 * ((Nop - L * I) @ (Nop - L * I)) + 30 * (Sz @ Sz)
    _, V0 = eigsh(Hpen, k=1, which='SA'); psi0 = V0[:, 0]
    E0 = float(np.real(np.vdot(psi0, Ham @ psi0)))
    Ev, Vv = np.linalg.eigh(Ham.toarray())
    # site density and spin-z operators
    ndens = [cd[2 * j] @ c[2 * j] + cd[2 * j + 1] @ c[2 * j + 1] for j in range(L)]
    szsite = [0.5 * (cd[2 * j] @ c[2 * j] - cd[2 * j + 1] @ c[2 * j + 1]) for j in range(L)]
    qs = list(range(1, L))                       # q = 2 pi m / L, skip q=0 (trivial)
    wgrid = np.linspace(0, 12, NW); eta_c, eta_s = 0.30, 0.16
    Sc = np.zeros((len(qs), NW)); Ss = np.zeros((len(qs), NW))
    for iq, m in enumerate(qs):
        q = 2 * np.pi * m / L
        rho_q = sum((np.exp(-1j * q * j) * ndens[j] for j in range(L)),
                    sp.csr_matrix((2 ** (2 * L),) * 2, dtype=complex))
        sz_q = sum((np.exp(-1j * q * j) * szsite[j] for j in range(L)),
                   sp.csr_matrix((2 ** (2 * L),) * 2, dtype=complex))
        for op, S, eta in ((rho_q, Sc, eta_c), (sz_q, Ss, eta_s)):
            ov = Vv.conj().T @ (op @ psi0)
            wt = np.abs(ov) ** 2; om = Ev - E0
            for e, weight in zip(om, wt):
                if e > 1e-6 and weight > 1e-9:
                    S[iq] += weight * (eta / np.pi) / ((wgrid - e) ** 2 + eta ** 2)
    return np.array(qs) / L * 2, wgrid, Sc, Ss     # q in units of pi


def main():
    qq, ww, Sc, Ss = structure_factors()
    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.5), gridspec_kw={'wspace': 0.30})
    for col, (Sd, ttl, cap) in enumerate([
            (Sc, r'charge  $S(q,\omega)$', 'Mott charge gap'),
            (Ss, r'spin  $S^{zz}(q,\omega)$', 'gapless spinon continuum')]):
        vmax = np.percentile(Sd, 99.5)
        im = ax[col].imshow(Sd.T, origin='lower', aspect='auto', cmap=PS.WARMMAP, vmin=0, vmax=vmax,
                            extent=[qq.min(), qq.max(), ww.min(), ww.max()])
        ax[col].set_title(ttl, fontsize=12.5, color=PS.INK)
        ax[col].set_xlabel(r'momentum transfer  $q/\pi$')
        ax[col].set_ylabel(r'frequency  $\omega/t$') if col == 0 else ax[col].set_yticklabels([])
        cb = fig.colorbar(im, ax=ax[col], fraction=0.046, pad=0.02); cb.ax.tick_params(labelsize=8)
        ax[col].text(0.5, 11.0, cap, color='white', fontsize=9.5, style='italic', ha='left')
    # mark the charge gap and the gapless spin onset
    ax[0].annotate('', xy=(qq[len(qq) // 2], 0), xytext=(qq[len(qq) // 2], 4.0),
                   arrowprops=dict(arrowstyle='<->', color='white', lw=1.2))
    ax[0].text(qq[len(qq) // 2] + 0.06, 1.8, r'$\Delta_c$', color='white', fontsize=11)
    ax[1].plot([qq.min(), qq.max()], [0, 0], color='white', lw=1.0, ls=(0, (4, 3)), alpha=0.7)
    fig.suptitle(r'Spin--charge separation: collective response of the Hubbard model, on which the screen also operates',
                 y=1.01, fontsize=11.5)
    fig.tight_layout()
    fig.savefig(os.path.join(PAP, 'fig_sqw.pdf')); fig.savefig(os.path.join(PAP, 'fig_sqw.png'), dpi=200)
    plt.close(fig); print('saved fig_sqw ->', PAP)


if __name__ == '__main__':
    main()
